# 基于方差的隐藏混杂敏感性分析

中文 | [English](README.md)

本仓库是一个可复现的 R/Python 项目，用于研究隐藏混杂下的 **variance-based sensitivity analysis, VBM**。主应用来自 Huang & Pimentel (Biometrika, 2025) 中的 NHANES 鱼类/贝类摄入与血汞水平案例。

核心问题是：

> 如果存在未观测混杂因素，这个隐藏混杂需要多强，才能推翻“高鱼贝摄入与更高血汞水平相关”的结论？

本仓库复现了主要的 VBM/MSM 分析，并加入四个 extension，用来 stress-test 不同的 hidden-confounding story。

---

## 项目包含什么

主程序会完成：

1. 读取 NHANES fish-consumption 数据；
2. 使用 inverse probability weighting 估计 ATT；
3. 计算 VBM sensitivity threshold 和 covariate benchmark；
4. 将 VBM 与 marginal sensitivity model, MSM benchmark 进行比较；
5. 自动生成基础图和 enhanced figures。

四个 extension 脚本如下：

| 脚本 | 作用 |
|---|---|
| `extension_01_hidden_strength_vbm_msm.R` | 比较不同 hidden-confounder 强度下 VBM 与 MSM 的表现。 |
| `extension_02_vbm_ps_misspecification.R` | 研究 observed propensity-score model 错设时 VBM 的表现。 |
| `extension_03_good_overlap_instability.R` | 在 very good overlap 情形下 stress-test VBM denominator instability。 |
| `extension_04_rgm_vbm_msm.R` | 探索 RGM / TV-L1 sensitivity model，并与 VBM 和 MSM 比较。 |

---

## 仓库结构

```text
.
├── main.R
├── setup.R
├── extension_01_hidden_strength_vbm_msm.R
├── extension_02_vbm_ps_misspecification.R
├── extension_03_good_overlap_instability.R
├── extension_04_rgm_vbm_msm.R
├── functions/
├── rgm_model/
├── scripts/
└── output/
```

本项目有意不再使用 `extensions/` 源码文件夹。四个 extension 脚本直接放在仓库根目录，方便读者快速找到并单独运行。

---

## 快速开始

克隆仓库后，在 R 中打开项目根目录，运行：

```r
source("setup.R")   # 可选：安装缺失的 R 包
source("main.R")    # 主分析 + enhanced plots + 四个 extensions
```

默认情况下，`main.R` 会依次运行：

1. 主程序 NHANES 分析；
2. `scripts/enhanced_plots.py`，生成主程序 enhanced figures；
3. 四个 extension 脚本；
4. 每个 extension 自己对应的 enhanced plotting script。

输出会重新生成到：

```text
output/tables/
output/figures/
output/figures_enhanced/
output/extension_results/
```

生成结果默认被 Git 忽略。重新运行 `source("main.R")` 即可复现全部表格和图片。

---

## 单独运行某个 extension

```r
source("extension_01_hidden_strength_vbm_msm.R")
source("extension_02_vbm_ps_misspecification.R")
source("extension_03_good_overlap_instability.R")
source("extension_04_rgm_vbm_msm.R")
```

---

## 配置项

配置文件位于：

```text
functions/config.R
```

主要默认开关包括：

```r
run_main_enhanced_plots = TRUE
run_extensions = TRUE
run_extension_01_hidden_strength_vbm_msm = TRUE
run_extension_02_vbm_ps_misspecification = TRUE
run_extension_03_good_overlap_instability = TRUE
run_extension_04_rgm_vbm_msm = TRUE
```

如果只想运行较短流程，可以把对应开关设为 `FALSE`。

---

## 依赖

R 包：

- `CrossScreening`
- `survey`
- `boot`
- `ggplot2`
- `dplyr`
- `tidyr`
- `readr`
- `scales`

enhanced plots 需要的 Python 包：

- `pandas`
- `numpy`
- `matplotlib`

安装 Python 包：

```bash
pip install pandas numpy matplotlib
```

---

## 核心思想

Sensitivity analysis 的作用不是消除所有隐藏假设，而是让隐藏假设变得可见。本项目比较了三种描述隐藏混杂的方式：

- **MSM / L-infinity**：防止最极端个体的权重扭曲。
- **VBM / L2-R2**：衡量系统性的 residual weight variation。
- **RGM / L1-TV**：衡量隐藏混杂最多能移动多少 probability mass。

NHANES 复现部分验证了主结果与论文接近；四个 extension 进一步说明，VBM 的解释依赖于 propensity-score model 的设定、权重方差的稳定性，以及我们选择哪一种 norm 来刻画 hidden perturbation。




---

## 参考文献

Huang & Pimentel (2025). Variance-based sensitivity analysis. *Biometrika*.
