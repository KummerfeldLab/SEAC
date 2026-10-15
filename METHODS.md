# SEAC 方法总结

## BLADE：本包实现的范围

SEAC 重现 `VisHD-main/SEAC.py::add_blade_column()` 中用于 Visium HD 的局部
BLADE 掩膜，而不是重做已发表 BLADE 工作中的完整切片级伪影检测流程。

完整的 `tissue_positions.parquet` 中，capture array 矩形边界上的位置以及
`in_tissue == 0` 的位置均作为种子。将所有种子的 `(array_row, array_col)`
坐标建立为一棵 KD-tree；再对每个位置查询最近种子的精确欧氏距离。距离小于
或等于 1 的位置标为 `blade` 并在 QC 后移除。

因此，这里不是“每个 spot 分别找邻居”的规则，也不是 `k = 8` 的 kNN 模式。
对角位置距离为 sqrt(2)，不会仅因与种子对角相邻而被移除。

## MaxMin：通用基因集空间共定位统计量

用户提供至少两个具名、彼此不重叠的基因集。对每一个基因集，SEAC 找到 QC 后
至少含一个该基因集 UMI 的 bin，并用这些 bin 构造 KD-tree。每个待评估 bin 到
该基因集最近阳性 bin 的欧氏距离向上取整。一个 bin 的 MaxMin 是所有基因集
距离中的最大值：数值越小，表示所有指定基因集在该位置附近越紧密地共同出现。

该统计量不依赖 SEAC 基因表。`seacGeneSets` 仅作为可运行的示例；任意符合输入
规则的基因集均可用于 `maxMinScore()` 和 `maxMinPermutationTest()`。

基因集既可以使用 feature row name（例如 Ensembl ID），也可以使用 10x HDF5
中保存在 `rowData(x)$Symbol` 的基因符号。同一符号对应多个 feature row 时会全部
计入；若 ID/符号别名导致两个基因集最终映射到同一 feature，则报错以避免重复
计算同一 UMI。

## UMI-level permutation test

置换保留每个 bin 的总 UMI 数，并将属于各输入基因集的 UMI 类别与 `Other` 类别
随机重新分配。`within` 模式保留每个样本内各基因集的 UMI 总数；`pooled` 模式
在多个样本间汇总后保留全局类别总数。每次置换重新计算 MaxMin，但只累计小于或
等于观测统计量的次数，因此不会将完整 null matrix 放入内存。

输出包括 Monte Carlo p 值、BH FDR、原始管线兼容的 1% lower-tail hotspot 标记，
以及可恢复的分块 checkpoint。

## 输出与可重复性

`runBladeQC()` 写入 HDF5-backed `SpatialExperiment`、完整 BLADE audit mask 和
manifest。`maxMinPermutationTest()` 将紧凑结果写为 HDF5；不会复制原始 counts。
每个 permutation 有确定的独立随机种子，因此在不改变参数或输入签名时，改变
chunk 大小或并行调度不会改变结果。

Visium HD 位置表由轻量的 `nanoparquet` 读取，HDF5 counts 由 `TENxIO` 以
`TENxMatrix` 磁盘后端导入，不需要安装 Apache Arrow C++ runtime。置换 checkpoint
采用输入内容哈希和临时文件原子改名；若最近一次写入被中断，可回退到上一份完整
checkpoint。
