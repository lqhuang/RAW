# RAW Wrapper
Get tired of trivial GUI operations? Try these scripts to make life easier.

[English](./README.md), 中文简介]

# 使用指南

适用于在 `SSRF` 采集的实验数据


实验目录结构:

    ONE-ROOT-DIR
    |- Data
    |- Processed
    |- Averaged
    |- Subtracted
    |- GNOM
    |- config.yml
    |- setup.yml


`ProcessedFilePath`: `default=Processed`
`AveragedFilePath`: `default=Averaged`
`SubtractedFilePath`: `default=Subtracted`
`GnomFilePath`: `default=GNOM`

在计算过程中，如果这些目录不存在会自动创建。但会报错`NotADirectoryError`如果目标文件夹不是一个文件夹。

`skip_frames`: `default=1`

在 SSRF 的实验采集中，若使用自动曝光方式，将产生 `ionchamber`（后缀`.Iochamber`），记录文件后缀 `.log` 以及相对应的 TIFF 图像文件后缀 `.tif`。
