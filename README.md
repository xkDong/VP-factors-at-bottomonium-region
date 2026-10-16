# PaperISR.jl：论文复现与数值更新

**直接使用结果：** [读者结果包](reader_results/README.md) 整理了新旧版本全部 209 个能点的 ISR、VP 因子及误差 CSV、PNG/PDF 图、新版协方差矩阵和读取示例。该目录可独立分发，读取结果或重画图无需外接硬盘和旧仓库。

## 主要结果：新旧版本对比

蓝色虚线为文章复现版本，橙色实线为数值更新版本。四张图均比较相同定义的**名义输入结果**；下方小图为 $100(X_{\rm new}/X_{\rm old}-1)\%$。曲线连接 209 个测量能点以方便阅读，不代表额外的连续谱拟合。论文采用的抽样均值另存于 CSV，详见 [更新说明](reader_results/CHANGES.md)。

### ISR 修正因子

展示逆 ISR 因子 $C_{\rm ISR}=1/(1+\delta)$。两版名义值大部分重合，主要差别在低能点；最大相对变化约为 0.618%。

![新旧版本逆 ISR 因子对比](reader_results/figures/github_isr.png)

[逐点数据](reader_results/data/comparison_isr.csv) · [矢量 PDF](reader_results/figures/github_isr.pdf)

### Dressed 截面

截面单位为 nb。两版均用光滑前的原始实验截面乘以各自逆 ISR 因子；光滑谱仅用于求解修正因子，不替代最终报告的实验数据。

![新旧版本 dressed 截面对比](reader_results/figures/github_sigma_dressed.png)

[逐点数据](reader_results/data/comparison_sigma_dressed_nb.csv) · [矢量 PDF](reader_results/figures/github_sigma_dressed.pdf)

### Born 截面

两版均按 $\sigma_{\rm Born}=\sigma_{\rm obs}^{\rm raw}C_{\rm ISR}/F_{\rm VP}$ 计算，单位为 nb。均值和误差来自逐样本比值，保留原始输入、ISR 与 VP 的相关性；旧版 Born 是由已验证的历史样本推导的结果，尚未与原论文独立的 Born 表逐点核对。

![新旧版本 Born 截面对比](reader_results/figures/github_sigma_born.png)

[逐点数据](reader_results/data/comparison_sigma_born_nb.csv) · [矢量 PDF](reader_results/figures/github_sigma_born.pdf)

### VP 因子与误差带

阴影分别表示各版名义中心值上下一个**输入误差传播标准差**。旧版标准差包含 10000 项，新版排除首个名义样本，使用 9999 个随机样本。误差带不包含光滑、pQCD 等额外模型不确定度，也不是均值的蒙特卡洛标准误。两版输入相关，误差带重叠程度不能直接用于独立测量的显著性判断。

![新旧版本 VP 因子及输入误差带](reader_results/figures/github_vp.png)

[逐点数据](reader_results/data/comparison_vp.csv) · [矢量 PDF](reader_results/figures/github_vp.pdf)

## 程序与分支

当前为 numerical-updates 分支。paper-reproduction 分支固定保存逐点对齐论文附录的历史重放，初始提交为 5cabb61。

更新实现修复明确的公式和单位错误，采用带残差检查的 ISR 线性反演、VP 收敛求解、一次按能量阈值配置的平滑以及逐样本协方差传播。新接口位于 src/updated.jl，历史接口不变。

    julia --project=. test/runtests.jl
    julia --project=. scripts/audit_updates.jl
    julia --project=. scripts/check_updated_methods.jl
    julia --project=. scripts/validate_updates.jl
    julia --project=. scripts/export_observed_cross_sections.jl

更新结果在 results/updates/；最终采用范围以 acceptance.toml 为准。负的反演估计值未截零，平滑及 pQCD 等模型误差不包含在输入误差协方差内。详见 [更新检查报告](docs/numerical_updates.md)。

其中 `results/updates/` 的 dressed/Born 列保留原先求解器光滑谱的审计记录。最终乘回原始实验数据的截面及协方差在 `results/observed_cross_sections/`，读者表和主页图采用后者；由 `scripts/export_observed_cross_sections.jl` 重新计算。

## 历史复现入口

独立 Julia 包，按大仓库 MATLAB 实现重算 [arXiv:2002.09838](https://arxiv.org/abs/2002.09838) 最终附录的逆 ISR、VP 因子及其误差。

**状态：ISR 和 VP 历史样本已完整重算，附录全部 $209$ 点验证通过，测试 $25/25$ 通过。** 现已包含完整 ISR → dressed $R_b$ → VP 重放；另已提供原始实验数据乘回修正因子的 dressed/Born 汇总和配对误差检查；共振拟合尚未整体验收。旧项目未被替换。

## 结构

- `src/`：解析积分、普通 LOWESS、逐能点迭代、随机抽样。
- `scripts/reproduce.jl`：ISR 历史样本重放和新种子计算。
- `scripts/reproduce_vp.jl`：默认重新计算 ISR 后接续 VP，或仅重放保存的 VP 输入。
- `scripts/vp_convention_check.jl`：原程序与论文强子虚部约定的单项差异检查。
- `scripts/export_legacy.py`：一次性从原 MATLAB 文件提取输入。
- `data/`：截面、误差、历史随机截面及来源记录。
- `reference/`：仅供验证的附录及原程序结果。
- `test/`：数值核与历史名义结果检查。
- `docs/`：方法、输入和未决事项。
- `results/`：本包新计算的结果及验证报告。

## 运行

仅依赖 Julia 标准库，开发测试版本为 Julia 1.12.6。在本包目录运行：

```sh
julia --project=. test/runtests.jl
julia --project=. --threads=4 scripts/reproduce.jl historical
julia --project=. --threads=4 scripts/reproduce.jl fresh 10000 2020
julia --project=. --threads=4 scripts/reproduce_vp.jl
julia --project=. scripts/reproduce_vp.jl saved
```

输入已随包提供，运行无需 MATLAB、Python、外置硬盘或旧项目。历史模式从原始随机截面重新平滑和反演；保存的 ISR 因子仅用于检验，不参与求解。新种子模式重新抽样，不应逐位匹配历史有限样本均值。

输出分别在 `results/historical/` 和 `results/fresh/`；最终验证报告为 `validation.toml`。未通过检查时不会写出新的成功报告。

## 对齐约定

使用历史网格，阈值为 $10.55852\,\mathrm{GeV}$，截面换算为 $\sigma=86.85R/s$（nb）。从第 $23$ 点，即 $10.6659\,\mathrm{GeV}$，开始对截面做普通 LOWESS，比例 $0.045$，无稳健重加权，仅初始化平滑。以 $s=E^2$ 插值，解析积分，按原程序逐能点求解，保留 $30$ 状态、$29$ 次更新。

附录中心值取统计抽样的逆因子均值；误差取总误差抽样的逆因子样本标准差。名义输入输出不等于抽样均值。报告另列 $30$ 与 $60$ 状态的敏感性，不将历史设置自动称为收敛极限。

## 一次性导入

维护者连接含 VP_B 的外盘后，用装有 NumPy/SciPy 的 Python 运行 `python scripts/export_legacy.py`。脚本只读原文件，导出实际计算输入并记录 SHA-256。导入后应将 `data/` 一并分发。

二进制格式为 little-endian Float64，每个样本的全部能点连续排列；形状写入 `data/provenance.json`。Julia 读取为“能点 × 样本”。每个数组约 $16.8\,\mathrm{MB}$。

附录中心值最大绝对差为 $8.39\times10^{-12}$，误差最大绝对差为 $2.57\times10^{-13}$。ISR 测试和新种子两组各 $10000$ 次已运行。VP 全链对比最大绝对差：中心值 $7.33\times10^{-15}$，误差 $7.49\times10^{-18}$。详细证据见 [验证记录](docs/reproduction.md)。

VP 复现约定、输入链和原程序与论文公式的差异见 [VP 验证说明](docs/vp_reproduction.md)。该历史模式保留原程序的强子虚部系数，不等于逐式实现论文公式。

尚未选择开源许可证；本地整理不等于已上传 GitHub。
