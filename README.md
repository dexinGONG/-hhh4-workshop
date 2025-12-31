

# hhh4-workshop

**时空多成分模型（endemic–epidemic / hhh4）培训代码与示例**

本仓库提供一套**可直接运行、完整复现**的 R 代码，用于教学演示
**endemic–epidemic（hhh4）时空统计模型**在传染病监测与预警中的应用。

本项目以**教学与实践为导向**，重点放在模型结构理解、数据组织方式以及结果解读，而非数学推导。

---

## 一、项目简介

endemic–epidemic（hhh4）模型是一类常用于**多地区传染病监测数据**的统计模型，可将病例数分解为三类传播来源：

* **地方性成分（Endemic）**：地区的基线发病风险
* **时间自回归成分（Autoregressive）**：本地疫情的延续传播
* **空间传播成分（Neighbourhood）**：来自周边地区的输入传播

该模型已被广泛用于流感、麻疹、手足口病等传染病的时空分析与预警研究。

---

## 二、教学目标

完成本 workshop 后，学员应能够：

* 理解 hhh4 模型的**核心思想与三成分结构**
* 掌握将监测数据整理为 `sts` 对象的方法
* 构建并拟合基础的 hhh4 时空模型
* 逐步扩展模型：

  * 加入协变量
  * 更换空间权重（如 power-law）
  * 引入地区随机效应
* 从**疾控业务视角**解读模型结果，用于监测与预警支持

---

## 三、仓库结构

```text
hhh4-workshop/
├── hhh4_workshop_master.R    # 主脚本：运行此文件即可跑通全部流程
├── data/
│   ├── cases_weekly.csv      # 周病例数（week × region）
│   ├── pop.csv               # 各地区人口
│   ├── W_adj.csv             # 地区邻接矩阵（0/1）
│   └── covariates_weekly.csv # 示例协变量（如易感人群比例）
├── map/
│   └── measles_map.gpkg      # 地区空间矢量数据
└── README.md
```

---

## 四、如何运行（运行入口）

### 1️⃣ 设置工作目录

将 R 的工作目录设为本仓库所在文件夹：

```r
setwd("path/to/hhh4-workshop")
```

### 2️⃣ 一键运行完整流程
双击打开"时空多成分模型.Rproj"文件

脚本将自动完成：

1. 安装并加载所需 R 包
2. 读取示例数据
3. 构建 `sts` 对象
4. 拟合一系列 hhh4 模型
5. 输出关键图形与结果

---

## 五、数据说明

本仓库中的数据为**教学示例数据**，用于演示模型流程。

* **cases_weekly.csv**
  各地区每周报告病例数

* **pop.csv**
  各地区人口规模（用于 offset）

* **W_adj.csv**
  地区邻接关系矩阵（0/1）

* **covariates_weekly.csv**
  示例时间变化协变量（如易感人群比例）

* **measles_map.gpkg**
  地区空间边界，用于空间可视化与邻接构建

---

## 六、建模流程说明

主脚本按以下**递进逻辑**组织（与实际研究流程一致）：

1. **基础模型**
   Endemic + AR + NE 三成分

2. **加入协变量**
   协变量进入 endemic 成分

3. **更换空间权重**
   由 0/1 邻接 → power-law 衰减

4. **加入随机效应**
   描述地区间异质性

5. **模型对比**
   使用 AIC 比较模型拟合效果

---

## 七、适用对象

本仓库适合以下人群使用：

* 疾控中心从事传染病监测与预警的工作人员
* FETP 学员与授课教师
* 流行病学与公共卫生研究人员
* 对时空传染病建模感兴趣的 R 用户

建议具备**基础 R 使用经验**。

---

## 八、复现与教学说明

* 所有分析均可通过**单一脚本完全复现**
* 代码注释突出**关键参数与建模决策**
* 更适合教学与培训，而非追求计算效率

---

## 九、许可证（License）

本项目采用 **MIT License**，允许自由使用、修改和再发布（需保留署名）。

---

## 十、致谢

hhh4 模型由 R 包 **`surveillance`** 实现。
在科研或实际应用中使用该模型时，请引用相应的方法学文献。

---

---

# hhh4-workshop

**Teaching materials for endemic–epidemic (hhh4) spatiotemporal models**

This repository provides reproducible R code for a hands-on workshop on
**endemic–epidemic (hhh4) models** for infectious disease surveillance.

The focus is on **practical implementation, interpretation, and teaching**, rather than mathematical derivation.

---

## Overview

The endemic–epidemic (hhh4) model decomposes disease incidence into:

* **Endemic component** (baseline risk)
* **Autoregressive component** (local temporal persistence)
* **Neighbourhood component** (spatial transmission)

It is widely used for spatiotemporal analysis of infectious disease surveillance data.

---

## Learning Objectives

Participants will learn how to:

* Understand the structure of hhh4 models
* Prepare surveillance data in `sts` format
* Fit and extend endemic–epidemic models
* Interpret results for public health surveillance and early warning

---

## How to Run
Double click .Rproj file.

Run the script from the repository root directory.

---

## License

MIT License.

