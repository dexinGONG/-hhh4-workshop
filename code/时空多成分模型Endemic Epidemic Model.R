# 时空多成分模型（endemic-epidemic / hhh4）
# 运行入口：解压文件包，双击"时空多成分模型.Rproj"，运行本code即可
# 目录结构：data/（csv）与 map/（measles_map.gpkg）与本脚本同级

# 00 环境与依赖 ------------------------------------------------
# 首次运行：若安装包很慢/失败，可先在 RStudio 设置镜像（Tools -> Global Options -> Packages -> CRAN mirror）
# 也可手动：options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, surveillance, spdep, sf, sp)

# 01 数据载入 ------------------------------------------------
cases_wide <- read_csv("data/cases_weekly.csv")        # week + 各地区病例（宽表）
pop_df     <- read_csv("data/pop.csv")                # unit + population
W_wide     <- read_csv("data/W_adj.csv")              # unit + 邻接矩阵（宽表，0/1）
cov_long   <- read_csv("data/covariates_weekly.csv")  # week + unit + Sprop（易感人群:1-疫苗接种比率）

units <- setdiff(names(cases_wide), "week")           # 地区代码（来自病例宽表列名）

# 邻接矩阵 W（地区×地区）：行名必须是 unit，且顺序与 units 对齐
W_mat <- as.matrix(W_wide[, units])
rownames(W_mat) <- W_wide$unit
W <- W_mat[units, units]

# 协变量 Sprop（周×地区矩阵）
# pivot_wider：把 long 表（week, unit, Sprop）展开为 wide 表（week + 每个 unit 一列）
Sprop_wide <- cov_long %>%
  select(week, unit, Sprop) %>%
  pivot_wider(names_from = unit, values_from = Sprop) %>%
  arrange(week)

Sprop <- as.matrix(Sprop_wide[, units])
colnames(Sprop) <- units

# Sprop 保险处理：避免 log(0) / NA 造成拟合报错或数值异常
Sprop[is.na(Sprop)] <- min(Sprop, na.rm = TRUE)
Sprop <- pmax(Sprop, 1e-6)   # 避免 log(0)

# （协变量来源说明）
# Sprop 来自疫苗接种率，若有更多数据，比如气象、环境、GDP等，直接照葫芦画瓢即可


# 02 构建 sts -------------------------------------------------
Y <- as.matrix(cases_wide[, units])                   # 病例矩阵：周×地区

# population（周×地区矩阵）：人口不随时间变 → 复制成与 Y 同维度
pop <- pop_df$population
names(pop) <- pop_df$unit
pop <- pop[units]
pop_mat <- matrix(rep(pop, each = nrow(Y)), nrow = nrow(Y), byrow = FALSE)

# 地图对象：用于空间可视化与邻接推导
map_sf <- st_read("map/measles_map.gpkg", quiet = TRUE)
measles_map <- tryCatch(as(map_sf, "Spatial"), error = function(e) as_Spatial(map_sf)) #（更稳定的转换写法）

# 让 map 的区域 ID 匹配 observed 列名（否则 sts / 空间图会报错）
row.names(measles_map@data) <- colnames(Y)
measles_map <- spChFIDs(measles_map, as.character(colnames(Y)))

# sts 对象：hhh4 的输入数据结构
stsObj <- sts(
  observed      = Y,
  start         = c(2001, 1),                         # 2001年第1周
  frequency     = 52,                                 # 周尺度
  population    = pop_mat,
  neighbourhood = W,                                  # 先放 0/1 邻接；power-law 时在这里替换即可
  map           = measles_map
)

# 快速检查
str(stsObj); dim(observed(stsObj)); head(epoch(stsObj))

plot(stsObj, type = observed ~ time) #时序图

plot(stsObj, type = observed ~ unit) #空间图

plot(stsObj) #按unit分面的时序图

# 03 邻接阶数矩阵 ---------------------------------------------
adjmat      <- poly2adjmat(stsObj@map)                # 0/1 邻接（由地图推导）
nbOrder_mat <- nbOrder(adjmat)                        # 1阶/2阶/3阶…邻接阶数

# power-law 的 maxlag：阶数太大更慢/更易不收敛；或按需调小/调大
maxlag_use <- max(nbOrder_mat)   # 最大邻接阶数（由地图决定）

# 04 基础模型：无协变量 ---------------------------------------
ctrl_base <- list(
  end = list(
    f = addSeason2formula(~ 1, period = 52),          # endemic：加季节项（52周周期）
    offset = population(stsObj)                       # offset：人口（对数偏移）
  ),
  ar = list(
    f = ~1                                            # AR：本地延续（最简形式）
  ),
  ne = list(
    f = ~1,
    weights = (nbOrder_mat == 1)                      # NE：仅1阶邻居（0/1权重）
  ),
  family = "NegBin1"                                  # 负二项：允许过度离散
)

fit_base <- hhh4(stsObj = stsObj, control = ctrl_base)
summary(fit_base, idx2Exp = TRUE, amplitudeShift = TRUE, maxEV = TRUE)

plot(fit_base, type = "season", components = "end", main = "")
plot(fit_base, type = "fitted", total = TRUE, hide0s = TRUE, legend = FALSE)

# fitted：挑高发地区看拟合
districts2plot <- which(colSums(observed(stsObj)) > 50)
plot(fit_base, type = "fitted", units = districts2plot, hide0s = TRUE, legend = 1)

# 提取 fitted 分量 + 三成分占比（end / AR / NE）
fc <- plot(fit_base, type = "fitted", total = TRUE, hide0s = TRUE, legend = FALSE)
head(fc$Overall)
colSums(fc$Overall)[3:5] / sum(fc$Overall[, 1])

# 过度离散参数区间（负二项的关键参数）
confint(fit_base, parm = "overdisp")

# Poisson 对照：说明为什么常用 NegBin（AIC 通常更优）
AIC(fit_base, update(fit_base, family = "Poisson"))

# 05 模型拓展：加协变量 ---------------------------------------
ctrl_cov <- list(
  end = list(
    f = addSeason2formula(~ 1 + log(Sprop), period = 52), # Sprop 进入 endemic（取 log）
    offset = population(stsObj)
  ),
  ar = list(f = ~1),
  ne = list(f = ~1, weights = (nbOrder_mat == 1)),
  family = "NegBin1"
)

fit_cov <- hhh4(stsObj = stsObj, control = ctrl_cov)
summary(fit_cov, idx2Exp = TRUE, amplitudeShift = TRUE, maxEV = TRUE)

# 也可以用 update() 来“在已有模型上改一小部分”（更省代码），例如：
# fit_cov2 <- update(fit_base, end = list(f = update(formula(fit_base)$end, ~ . + log(Sprop))))

# 06 更改空间权重：power-law -----------------------------------
stsObj_order <- stsObj
neighbourhood(stsObj_order) <- nbOrder_mat            # W_powerlaw 需要“邻接阶数矩阵”

ctrl_powerlaw <- list(
  end = list(
    f = addSeason2formula(~ 1 + log(Sprop), period = 52),
    offset = population(stsObj_order)
  ),
  ar = list(f = ~1),
  ne = list(
    f = ~1,
    weights = W_powerlaw(maxlag = maxlag_use)         # 幂律衰减权重（含参数 d）
  ),
  family = "NegBin1",
  optimizer = list(control = list(iter.max = 500, eval.max = 500))
)

fit_powerlaw <- hhh4(stsObj = stsObj_order, control = ctrl_powerlaw)
summary(fit_powerlaw, idx2Exp = TRUE, amplitudeShift = TRUE, maxEV = TRUE)

# 幂律指数 d：越大→衰减越快（更本地化）
cbind("Estimate" = coef(fit_powerlaw), confint(fit_powerlaw))["neweights.d", ]

plot(
  fit_powerlaw,
  type = "maps",
  which = c("epi.own", "epi.neighbours", "endemic"),
  prop = TRUE,
  labels = list(cex = 0.6)
)

# 07 模型拓展：加随机效应 -------------------------------------
ctrl_ri <- list(
  end = list(
    f = addSeason2formula(~ 1 + log(Sprop) + ri() - 1, period = 52), # endemic：地区随机效应
    offset = population(stsObj_order)
  ),
  ar = list(
    f = ~ 1 + ri() - 1                                              # AR：地区随机效应
  ),
  ne = list(
    f = ~ 1 + ri() - 1,                                             # NE：地区随机效应
    weights = W_powerlaw(maxlag = maxlag_use)
  ),
  family = "NegBin1",
  optimizer = list(control = list(iter.max = 800, eval.max = 800)) 
  # 可能出现 non-convergence：可把 iter.max / eval.max 再调大；或尝试更换 optimizer
)

fit_ri <- hhh4(stsObj = stsObj_order, control = ctrl_ri)
summary(fit_ri, idx2Exp = TRUE, amplitudeShift = TRUE, maxEV = TRUE)

plot(fit_ri, type = "ri", component = "end", exp = TRUE, labels = list(cex = 0.6), main = "")
plot(fit_ri, type = "ri", component = "ar",  exp = TRUE, labels = list(cex = 0.6), main = "")
plot(fit_ri, type = "ri", component = "ne",  exp = TRUE, labels = list(cex = 0.6), main = "")

# 08 模型对比：AIC --------------------------------------------
AIC(fit_base, fit_cov, fit_powerlaw, fit_ri)
