# PaperISR.jl：ISR 与 VP 复现

独立 Julia 包，按大仓库 MATLAB 实现重算 [arXiv:2002.09838](https://arxiv.org/abs/2002.09838) 最终附录的逆 ISR、VP 因子及其误差。

**状态：ISR 和 VP 历史样本已完整重算，附录全部 $209$ 点验证通过，测试 $25/25$ 通过。** 现已包含完整 ISR → dressed $R_b$ → VP 重放；尚未对 Born 截面汇总和共振拟合做整体验收。旧项目未被替换。

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
