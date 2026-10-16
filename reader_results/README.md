# ISR and VP results / 新旧版本结果包

This directory is self-contained. Read the CSV files directly; Julia, MATLAB, the original repository and the external disk are not required. Rebuilding the tables and plots requires only Python, NumPy and Matplotlib:

```sh
python -m pip install -r requirements.txt
python make_results.py
```

本目录可独立分发。CSV 可直接用 Excel、Python 或 Julia 读取；已有 PNG 和矢量 PDF 可直接使用。重画图不需要 Git、外接硬盘或本机 Codex 工具。

**了解新旧差别：** [更新逐项说明](CHANGES.md) 对照列出公式、算法、光滑、中心值及误差处理的变化，说明验证结果、实际影响和未纳入的不确定度。

## Start here / 文件索引

GitHub 主页使用四张主对比图：[ISR](figures/github_isr.png)、[dressed 截面](figures/github_sigma_dressed.png)、[Born 截面](figures/github_sigma_born.png)、[VP 与误差带](figures/github_vp.png)；均有同名 PDF。四图统一使用名义中心值，下方面板为新版相对旧版的百分比变化。蓝虚线为旧版、橙实线为新版；连接线仅辅助阅读。VP 阴影为名义值上下一个输入传播标准差，不是均值的标准误或包含全部模型误差的置信带。图中用到的名义值和标准差也汇总在 [github_series.csv](data/github_series.csv)。

| Version | ISR data | VP data | Factors and errors |
|---|---|---|---|
| Paper reproduction / 文章复现 | [CSV](data/paper_reproduction_isr.csv) | [CSV](data/paper_reproduction_vp.csv) | [PNG](figures/paper_reproduction.png), [PDF](figures/paper_reproduction.pdf) |
| Numerical updates / 数值更新 | [CSV](data/numerical_updates_isr.csv) | [CSV](data/numerical_updates_vp.csv) | [PNG](figures/numerical_updates.png), [PDF](figures/numerical_updates.pdf) |
| Same-definition comparison / 同定义对比 | [CSV](data/comparison_isr.csv) | [CSV](data/comparison_vp.csv) | ISR: [PNG](figures/comparison_isr.png), [PDF](figures/comparison_isr.pdf); VP: [PNG](figures/comparison_vp.png), [PDF](figures/comparison_vp.pdf) |

All factor tables contain the same 209 measured energy points, in increasing energy order. The artificial threshold point is excluded. Energies are in GeV; factors and their standard deviations are dimensionless. CSV files are UTF-8, comma-separated, with a header row and full floating-point precision.

## Definitions / 定义与使用

### Dressed cross sections in nb / dressed 截面数据

| 数据 | 文件 |
|---|---|
| 旧版 dressed 截面 | [paper_reproduction_sigma_dressed_nb.csv](data/paper_reproduction_sigma_dressed_nb.csv) |
| 新版 dressed 截面 | [numerical_updates_sigma_dressed_nb.csv](data/numerical_updates_sigma_dressed_nb.csv) |
| 新旧并列表 | [comparison_sigma_dressed_nb.csv](data/comparison_sigma_dressed_nb.csv) |
| 新版 dressed 截面协方差 | [numerical_updates_sigma_dressed_covariance_nb2.csv](data/numerical_updates_sigma_dressed_covariance_nb2.csv) |

每张单版本表包含能量、名义截面、总输入误差抽样均值、传播标准差和换算用参考缪子截面。截面及标准差单位为 nb，协方差单位为 $\mathrm{nb}^2$；全部保留 209 个测量点。

$$\sigma_{\rm dressed}=R_{\rm dressed}\,\sigma_{\mu\mu}^{\rm ref},\qquad
\operatorname{Cov}(\sigma_i,\sigma_j)=\sigma_{\mu\mu,i}^{\rm ref}\sigma_{\mu\mu,j}^{\rm ref}\operatorname{Cov}(R_i,R_j).$$

换算严格沿用各分支约定：旧版 $\sigma_{\mu\mu}^{\rm ref}=86.85/s\;\mathrm{nb}$；新版使用代码中的 $\alpha=1/137.035999074$、$m_\mu=0.105658\,\mathrm{GeV}$ 和换算常数计算带缪子质量修正的参考截面，具体数值已逐点附在 `sigma_mumu_reference_nb` 列。换算常数的不确定度未额外传播。

旧版表来自保存并经过 ISR 重放验证的 dressed $R_b$ 样本，定义为历史ISR 修正因子乘回原始输入后的谱，**尚未进行 VP 阶段第二次光滑**；均值与标准差按历史规则包含全部 10000 项。新版表同样来自原始实验截面乘以新版ISR 修正因子，误差统计排除首个名义样本，使用 9999 项。旧版该表的均值是总输入误差抽样均值，不是旧 ISR 因子表的统计误差抽样均值。误差由逐样本 dressed 结果计算，未将输入误差与 ISR 因子误差当作独立误差相加。

The ISR column is the ISR correction factor:

$$C_{\rm ISR}=\frac{1}{1+\delta}=\frac{\sigma_{\rm dressed}^{\rm solver}}{\sigma_{\rm observed,smoothed}},\qquad \sigma_{\rm dressed}^{\rm reported}=\sigma_{\rm obs,raw}C_{\rm ISR}.$$

The VP factor is:

$$F_{\rm VP}=|1-\Pi|^{-2},\qquad R_{\rm Born}=R_{\rm dressed}/F_{\rm VP}.$$

光滑输入用于求修正因子；最终实验结果采用光滑前的原始截面乘以修正因子。当需要新版 Born 结果时，直接使用配套的 [dressed/Born 表](data/numerical_updates_cross_sections.csv)，其中保存的是无量纲的 $R$ 值（不是单位为 nb 的截面）。

| Column | Meaning / 含义 |
|---|---|
| `energy_GeV` | Center-of-mass energy / 质心能量 |
| `factor_nominal` | Result from nominal input / 名义输入计算值 |
| `factor_stat_toy_mean` | Mean of historical statistical-input toys; old ISR only / 旧 ISR 统计误差抽样均值 |
| `factor_total_toy_mean` | Mean of total-input-error toys / 总输入误差抽样均值 |
| `factor_input_sd` | Sample standard deviation from total-input-error toys / 总输入误差传播标准差 |

To reproduce the paper's reported centers, use `factor_stat_toy_mean` for old ISR and `factor_total_toy_mean` for old VP. For the updated nominal calculation use `factor_nominal`. Means and nominal values are distinct estimators; do not mix them in a method comparison.

旧版本统计包含 10000 项，其中第一项为名义输入；新版均值及标准差使用其余 9999 个随机样本。旧 VP 名义值采用第 3 次迭代，与历史样本算法一致。所有误差条表示传播标准差，既不是均值的蒙特卡洛标准误，也不保证等于非高斯分布的 68% 区间。

The two comparison figures use **nominal versus nominal**, and total-input SD versus total-input SD. Each figure contains factor values, relative nominal changes, SD values and absolute SD changes. Absolute SD differences avoid dividing by the effectively zero ISR uncertainty at the first measured point. The separate version figures use the published mean convention for the old version and nominal centers for the new version; their center convention is written in the legends.

旧、新结果共用历史输入抽样，不能当作两组独立测量。误差范围沿用现有输入模型，不包含光滑方案、pQCD 模型/参数或新推断实验相关性带来的额外不确定度。低能区新版线性反演保留有符号的随机估计值，不对负样本截零。

## Born cross sections in nb / Born 截面

- [旧版 Born 截面](data/paper_reproduction_sigma_born_nb.csv)
- [新版 Born 截面](data/numerical_updates_sigma_born_nb.csv)
- [新旧并列表](data/comparison_sigma_born_nb.csv)
- [新版 Born 截面协方差，单位为平方 nb](data/numerical_updates_sigma_born_covariance_nb2.csv)
- [对比图 PNG](figures/github_sigma_born.png) · [PDF](figures/github_sigma_born.pdf)

每表列出 209 个能点的名义值、总输入误差抽样均值和标准差。两版最终截面都使用原始实验输入：

$$\sigma_{\rm dressed}^{(k)}=\sigma_{\rm obs,raw}^{(k)}C_{\rm ISR}^{(k)},\qquad
\sigma_{\rm Born}^{(k)}=\frac{\sigma_{\rm obs,raw}^{(k)}C_{\rm ISR}^{(k)}}{F_{\rm VP}^{(k)}}.$$

光滑只参与修正因子的求解。旧版使用全部 10000 项，更新版统计排除首个名义输入，使用 9999 个随机样本。不能用均值的比值代替比值的均值，也不能独立相加各项误差。旧 Born 表是本次从匹配的历史 dressed 与 VP 样本推导的表，不宣称已逐点复现一张独立的论文 Born 表。新版 paired Julia 重算记录见 [validation](source/new_applied_validation.toml)。ISR、VP 因子已与原更新结果核对保持一致。

## Covariance / 协方差

Updated covariance matrices are supplied for [ISR](data/numerical_updates_isr_covariance.csv), [VP](data/numerical_updates_vp_covariance.csv), [dressed R](data/numerical_updates_dressed_covariance.csv) and [Born R](data/numerical_updates_born_covariance.csv). Each CSV has one leading row-energy column, and 209 covariance columns whose headers give their energies. Thus, after the header, the numeric array is 209 by 210; discard its first column to obtain the matrix. Row and column ordering matches the factor tables.

$$\operatorname{Cov}_{ij}=\frac{1}{N-1}\sum_{k=1}^{N}(X_{ik}-\bar X_i)(X_{jk}-\bar X_j),\qquad N=9999.$$

The diagonal equals the squared reported SD. These are within-observable covariance matrices, not ISR–VP cross-covariances. For Born fits, use the supplied Born covariance, which already propagates the paired ISR and VP samples. Do not combine their separate errors as if independent. Historical covariance matrices are not supplied; the old SD columns do not imply independent energy points.

[VP error components](data/numerical_updates_error_components.csv) separately list bottom-statistical-only, background-only and all-input VP SDs. Born component errors from the previous solver-spectrum convention are not presented as errors of the raw-data Born cross sections. The background ensemble does not allow a further reliable statistical/systematic split.

## Reading example / 读取示例

```python
import numpy as np

vp = np.genfromtxt('data/numerical_updates_vp.csv', delimiter=',', names=True)
energy = vp['energy_GeV']
factor = vp['factor_nominal']
error = vp['factor_input_sd']
raw = np.loadtxt('data/numerical_updates_vp_covariance.csv', delimiter=',', skiprows=1)
covariance = raw[:, 1:]
assert np.allclose(raw[:, 0], energy)
assert np.allclose(np.diag(covariance), error**2)
```

## Provenance and checks / 来源与检查

Sources are frozen from `paper-reproduction` at `5cabb61` and `numerical-updates` at `27d9bd1`; complete commit IDs, source paths and SHA-256 hashes are in [source/provenance.json](source/provenance.json). Original factor tables are retained in `source/`. The later raw-data cross-section calculation is separately recorded in `source/new_applied_validation.toml` and source hashes; it supersedes the solver-spectrum dressed/Born columns for reader-facing results.

The builder checks source hashes, prior numerical acceptance, energy alignment, finite values, matrix symmetry and covariance diagonals. Plot geometry, visual review and display acceptance results are archived in `verification/`. The ISR/VP factors are retained; dressed/Born cross sections and covariance were recomputed with paired Julia ensembles using the unsmoothed original experimental cross sections. `SHA256SUMS` covers this distribution.

名义 ISR 最大相对变化为 0.617727%，位于 $10.6063\,\mathrm{GeV}$；名义 VP 最大相对变化为 0.046090%，位于 $10.6659\,\mathrm{GeV}$。完整逐点变化保存在两张 comparison 表中。
