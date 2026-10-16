# 输入

- stat_input.csv、total_input.csv：对应运行实际使用的能量、名义截面和四类误差。
- stat_draws.f64、total_draws.f64：历史随机截面输入，不是 ISR 输出。
- reference 中 forward_factors：原 MATLAB 正 ISR 因子，仅供检验。
- reference/cpc.csv：最终发表版附录，仅用于比较。
- reference/central_replay.csv：前期审计的名义结果，仅用于比较。

历史计算输入已从原保存文件导入；数据来源与哈希记录在 data/provenance.json。没有从附录反推输入。
