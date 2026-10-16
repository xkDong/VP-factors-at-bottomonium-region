> 本文件保留前期审计记录，其中“当前 Julia”指旧 BottomoniumReproduction.jl 实现，不是本包。新包验收见 reproduction.md。

# 附录 ISR 来源确认

只读检查外盘 `VP_B/R_b/matlab/codes1`，脚本 `scripts/audit_matlab_isr.py` 可复查。未修改原始文件或当前 Julia 默认计算。

对 CPC 附录全部 $209$ 点，能量完全相同：

- 中心值来自 `ISR_factor_sta_err.mat` 的 `mean(1 ./ onedeltasamp, 1)`，最大绝对差 $8.46\times10^{-14}$。
- 误差来自 `ISR_factor_total_err.mat` 的ISR 修正因子样本标准差（样本标准差，分母 $N-1$），最大绝对差 $8.32\times10^{-14}$。
- 两组均为 $10000$ 个样本。中心值不是无扰动谱的结果，也不是总误差样本的均值；平均ISR 修正因子不能替换为平均因子的倒数。

独立 Python 实现了 MATLAB `g1.m` 解析积分、`ISRcalculation.m` 的普通 LOWESS，以及 `iteration.m` 的逐能量迭代。与总误差文件第一行（无随机扰动样本）的ISR 修正因子相比，最大绝对差 $6.61\times10^{-13}$。已设断言检验。不需要重新抽样即可验证中心值算法；未重新运行整套 MATLAB 或重新生成全部蒙特卡洛样本。

## 与当前 Julia 的关键不同

| 项目 | 恢复的 MATLAB | 当前 Julia |
| --- | --- | --- |
| 网格 | $209$ 个观测能量加一个阈值点 | $212$ 个观测，重复能量构造曲线时平均 |
| 阈值 | $10.55852\,\mathrm{GeV}$ | $10.5585\,\mathrm{GeV}$ |
| 转为截面 | $86.85R/s$ | 按程序物理常数计算 $\sigma_{\mu\mu}$ |
| 平滑对象 | 截面 | $R_b$ |
| 平滑起点 | MATLAB 第 $23$ 点，即 $10.6659\,\mathrm{GeV}$ | $E\geq10.66\,\mathrm{GeV}$，包含 $10.6613\,\mathrm{GeV}$ |
| 平滑方法 | 普通 LOWESS，比例 $0.045$，本数据为 $9$ 点邻域 | 三轮稳健重加权 LOWESS |
| 平滑频率 | 仅初始化 | 仅初始化 |
| 插值变量 | $s=E^2$ | $E$ |
| 求解顺序 | 逐能量求解；较低区间用第 $30$ 层已求结果，邻接区间用上一迭代层 | 全谱同步逐轮更新 |
| 迭代 | $30$ 个状态层，循环 $2:30$ | $20$ 次更新 |
| 积分 | `g1` 解析原函数 | 数值求积矩阵 |

直接调用链是 `ISR*.m -> ISRcalculation.m -> iteration.m -> g1.m`。`F.m` 和 `iterationint.m` 不在这条实际调用链中；不能用它们中的公式判断本次已验证结果。

在 $10.6613\,\mathrm{GeV}$，保存的无扰动结果为 $1.00309415$，统计样本均值（附录）为 $0.99614876$，误差为 $0.08918357$，当前 Julia 为 $1.19291348$。该点恰好处于两实现不同的平滑边界，但尚未进行逐项消融，不能声称单独一个设置解释全部差异。

逐点证据：`results/matlab_isr_audit/appendix_match.csv`、`central_replay.csv` 和 `summary.json`。
