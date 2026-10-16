# 原程序与 Julia 对应

| 原程序 | Julia | 约定 |
| --- | --- | --- |
| 保存文件的 datacombis | data 输入表 | 使用实际截面及误差 |
| Newsigma.m | sample_inputs | 相关系统误差、负值重抽样、首样本设为名义值 |
| ISRcalculation.m | smoothing_matrix | 普通局部线性 LOWESS |
| g1.m | primitive、segment | MATLAB dilog 对应 $\operatorname{Li}_2(1-x)$ |
| iteration.m | unfold、unfold_literal | 相邻区间用上一状态，较低区间用最终状态 |
| 保存样本的汇总 | ensemble_summary | 统计均值与总误差样本标准差 |

设 $q_k=(1+\delta_k)^{-1}$，附录使用

$$
\bar q_{\rm stat}=\frac1N\sum_kq_{k,\rm stat},\qquad
s_{q,\rm total}=\sqrt{\frac{\sum_k(q_{k,\rm total}-\bar q_{\rm total})^2}{N-1}}.
$$

缓存系数利用积分对端点截面的线性性。直接原函数求值路径用于对照；历史模式另逐样本、逐能点核对保存的 MATLAB 因子。
