# BES 数据替换工作计划与约定

> Historical input-only BES replacement, superseded by the current erratum results in `results/paper_method_corrected/`. See the repository README for current figures and data.

基线为 `paper-reproduction` 的 `5cabb61`。目标是在原数值实现下，将内部粲区输入换成 CPC 原文引用的 BES 包容 R 测量，重新运行 ISR、VP、dressed/Born 截面及抽样误差。

1. 使用 PDG 数据表中的 BES 2002 测量覆盖 $3.70$–$4.80\,\mathrm{GeV}$，用 BES 2000 的 $5.00\,\mathrm{GeV}$ 端点补齐。重复的 $4.60\,\mathrm{GeV}$ 点使用 2002 年结果。共 $77$ 点。
2. 用原始统计误差及 PDG 系统误差（百分数换算；2000 年各系统分量平方和）取代内部数据及其人为 $4\%$ 系统误差。
3. 原背景精确满足 $R_{i,t}=R_i+\delta R_{i,\rm stat}z_{t,\rm stat}+\delta R_{i,\rm sys}z_{t,\rm sys}$，重建残差低于 $2\times10^{-15}$。重用两条随机流，使每个误差分量跨能点完全相关；不把两种分量强行合为一个随机量。该模型是历史复现假设，不是 BES 实验完整协方差。
4. 新数据直接是包容 R，无需另加 $2.15$ 或开粲重子截面。pQCD 拼接起点从 $4.65$ 移至 $5.05\,\mathrm{GeV}$，避免与新增测量重叠。
5. ISR、底区原始样本、低能背景、两阶段平滑、历史虚部和常数全部保留。重新验证 ISR 不变，再计算新 VP。最终截面使用未光滑实验样本逐样本乘 ISR、除 VP。
6. 运行原测试、原结果重放、BES 输入检查、直接积分对照、VP 迭代敏感性检查；输出数据、误差协方差及对比图。

## 不确定性和限制

- 原文同时引用两篇 BES 测量，未指定重叠点合并规则；本分支明确采用较新测量，不作未知相关性的平均。
- 主表使用 PDG 精度，可能与论文四舍五入后的系统误差末位稍有差别。
- 原文“完全相关”按已恢复的历史双随机量模型解释；不等于实验本身的协方差。BES 2002 原文说明共同误差约 $3.3\%$，更精细误差模型留待独立更新。
- 原实现的强子虚部缺项、pion 换算常数、固定迭代和光滑起点问题不在本次输入替换中修复。此分支是数据来源一致性检查，不能称为全面修正的最终结果。
- 旧分支的 `results/vp` 等目录作为历史基线保存；当前结果唯一位置为 `results/besii`。
