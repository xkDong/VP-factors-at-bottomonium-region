# VP 数值复现

## 结果与输入链

已定位附录来源为大仓库 `matlab/vp_iterate/vp_smo_tot_err_0320.mat`。VP 中心值和误差均取该文件总误差样本的均值及样本标准差；与 ISR 的汇总方式不同。

默认入口 `scripts/reproduce_vp.jl` 使用历史 ISR 随机截面，重新计算逆 ISR 因子，然后按

$$
R_{b,k}^{\mathrm{dressed}}(s)=\frac{s}{86.85}\sigma_{k}^{\mathrm{observed}}(s)(1+\delta_k)^{-1}
$$

生成 VP 输入。与原 VP 输入的最大绝对差为 $3.59\times10^{-12}$。低能区和粲区使用 VP 保存的随机输入，未从附录反推数据。数据来源与哈希见 `data/vp/provenance.json`。

| 全链检查 | 最大绝对差 |
| --- | ---: |
| 全部 $10000\times210$ 个 VP 因子与 MATLAB | $4.89\times10^{-15}$ |
| 附录 $209$ 点均值 | $7.33\times10^{-15}$ |
| 附录 $209$ 点标准差 | $7.49\times10^{-18}$ |
| 名义输入五轮结果与 MATLAB | $4.45\times10^{-15}$ |

直接使用保存的 dressed 输入时，VP 样本差异仅 $8.89\times10^{-16}$。矩阵实现与直接梯形求和实现均已核对，前者还对全部历史样本逐点检验。结果见 [逐点对比](../results/vp/appendix_comparison.csv) 和 [验证报告](../results/vp/validation.toml)。

## 原程序约定

- VP 先对 dressed $R_b$ 再做普通 LOWESS，从第 $22$ 点，即 $10.6613\,\mathrm{GeV}$ 开始，比例 $0.045$。ISR 阶段的平滑起点和对象不同，不能合并两次平滑。
- 原始网格按 Newenergy 五等分，完整积分网格有 $5592$ 点，底区子网格有 $1046$ 点。积分以 $s=E^2$ 为变量。
- 轻子贡献、窄共振、pQCD、高能尾部和单位换算按原代码重写。窄共振电子宽度保留原程序中的强子分支修正。
- pQCD 的底阈值采用 RpQCD.m 中的 $5.62\,\mathrm{GeV}$；运行耦合中的参数仍为 $5.6\,\mathrm{GeV}$，不混用。
- 连续谱采用减去局部截面后的梯形积分；奇点处数组元素直接置零。这不是解析主值积分的同一离散实现。
- 每轮仅更新底夸克部分，固定轻夸克背景。历史随机样本取第 $3$ 轮；名义检查保留 $5$ 轮。名义输入第 $3$ 至第 $5$ 轮最大绝对变化为 $4.80\times10^{-8}$。
- 当前大仓库的 vp_iterate.m 已改成 smoothing spline；附录对应的历史保存结果采用普通 LOWESS。全样本平滑结果已核对，不以脚本修改日期猜测来源。

## 原代码与论文公式的实际差异

这是 **source-faithful 历史数值复现**，不等于论文逐式计算或精度升级版本。

论文的强子虚部写为

$$
\operatorname{Im}\Pi_h(s)=-\frac{s\,\sigma(s)}{4\pi\alpha}.
$$

原 `PIH.m` 配合 `vacuum.m` 的符号转换，实际得到

$$
\operatorname{Im}\Pi_h^{\rm legacy}(s)=-\frac{s\,\sigma(s)}{4\pi},
$$

即缺少 $1/\alpha$。代码中的截面由 sigmu.m 生成，单位换算也已显式核对，未在别处补回这个因子。复现模式保留它，才能对齐附录。论文原文摘录保存在 `reference/vp/paper_vp_excerpt.txt`。

另有原 pipi.m 的换算常数 $380379000$ 与 sigmu.m 使用的 $389379000$ 不同；当前历史模式同样保留原值，未静默纠正。

单独补回强子虚部的 $1/\alpha$，其余设置不变，名义 VP 最大相对变化约 $0.04424\%$，出现在 $10.5788\,\mathrm{GeV}$。这只是单项敏感性诊断，并非完整修正后结果，也未重新计算其抽样误差。见 [诊断数据](../results/vp/imaginary_convention.csv)。

## 使用

```sh
julia --project=. --threads=4 scripts/reproduce_vp.jl
julia --project=. scripts/reproduce_vp.jl saved
julia --project=. scripts/vp_convention_check.jl
```

分别输出到 `results/vp/`、`results/vp_saved/`。默认完整链和保存输入路径均不需要连接外盘。保持包名 PaperISR，以兼容已有 ISR 调用；新增 API 为 VPKernel、vp_replay、vp_literal。
