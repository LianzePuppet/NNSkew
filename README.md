# NNSkew

**NNSkew (dinucleotide skew)** 是一个计算二核苷酸偏倚度的命令行工具。

## 定义
给定二核苷酸 XY：
(NNSkew) = (#XY - #YX) / (#XY + #YX)

例如：
- GC skew = ( #GC - #CG ) / ( #GC + #CG )
- AT skew = ( #AT - #TA ) / ( #AT + #TA )

## 安装
```bash
pip install -e .
```

依赖：
- Python >= 3.9
- bedtools
- bedGraphToBigWig (UCSC tools)
- 可选：bedClip (UCSC tools)

## 用法
```bash
# 以 50bp 窗口计算 GC skew
nnskew -f hg38.fa -o hg38.gcskew.50bp.bw -b 50 --nn GC

# 以 100bp 窗口计算 AT skew
nnskew -f hg38.fa -o hg38.atskew.100bp.bw -b 100 --nn AT
```
