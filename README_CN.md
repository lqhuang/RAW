# RAW Wrapper

Get tired of trivial GUI operations? Try these scripts to make life easier.

[English](./README.md), 中文简介

## 环境依赖

可运行的最小依赖 (PY2/PY3):

    pip install numpy scipy pillow six pyyaml

通过 c++ 扩展(需编译)加快散射图像到 1D 曲线的计算，仅在 Python 2 + Linux/macOS 环境下可用:

    pip install weave

支持更多更全的散射图像格式:

    pip install pyfai fabio

快速安装:

    pip install -r requirements.txt
    # or
    conda install --file requirements.txt

conda env:

    conda create -n raw --file requirement.txt
    source activate raw

## 使用指南

适用于在 `SSRF` 采集的实验数据

### raw_cfg_printer 输出 raw_cfg 配置

美观地输出 `raw_settings.cfg` 的各项配置内容

    python raw_cfg_printer.py /path/to/raw_settings.cfg

重定向到 `raw_cfg.log` 文件中

    python raw_cfg_printer.py /path/to/raw_settings.cfg > raw_cfg.log

### raw_script 自动脚本

查看帮助

    python raw_script.py --help

首先需要确定出单次实验的根目录，并在根目录中创造 `config.yml` 配置文件(不放在根目录也可以，但需在配置文件中指明根目录的位置)，同时指定脚本所需要的参数，示例：

    # exp_root_path: 'blabla/SSRF-MagR-201805/Analysis/EXP41'
    raw_cfg_path: 'blabla/SSRF-MagR-201805/cfg/20180503-recentering-remasking.cfg'

    # Relative path to root directory
    # ProcessedFilePath: 'Processed'
    # AveragedFilePath: 'Averaged'
    # SubtractedFilePath: 'Subtracted'
    # GnomFilePath: 'GNOM'

    window_size: 1
    buffer_skip_frames: 1
    skip_frames: 1

执行脚本

    python raw_script.py /path/to/config.yml

查看执行记录

    less /path/of/run-root-directory/summary_timestamp.log

可以从命令行中补充特定的参数，并将**替代**配置文件中的参数，但不会写入到配置文件中，例如:

    python raw_script.py /path/to/config.yml --raw_cfg_path=/path/to/raw_settings.cfg --window_size=1

若需要将命令行补充的参数写入到配置文件中，可以通过 `--write_to_file` 的参数实现(由于 `pyyaml` 的实现原因，重新写入后会丢失掉文件中的所有注释型内容):

    python raw_script.py --write_to_file /path/to/config.yml --raw_cfg_path=/path/to/raw_settings.cfg --window_size=1

暂时只允许以 `--key=arg` 的方式进行解析。

如果没有在 `.yml/.yaml` 中配置或命令行中指明，部分参数将默认补全，以下为列表:

    exp_root_path     : 使用 .yml 所在的目录作为参数
    skip_frames       : 1  # 平均时单组内跳过 1 帧
    buffer_skip_frames: 1  # buffer 组平均时跳过 1 帧
    SourceFilePath    : 'Data'  # 相对于根目录的位置或文件名
    ProcessedFilePath : 'Processed'
    AveragedFilePath  : 'Averaged'
    SubtractedFilePath: 'Subtracted'
    GnomFilePath      : 'GNOM'
    img_ext           : '.tif'
    dat_ext           : '.dat'
    ionchamber_ext    : '.Iochamber'
    record_ext        : '.log'
    overwrite         : False  # 是否重新计算 img -> 1D profile 的过程，覆盖已存在的 1D profiles.

其他参数:

    raw_cfg_path      : /path/to/raw_settings.cfg (必要，缺少将报错)
    window_size       : 1  # 窗口大小，目前只有 `1` 和 `unset` 的两种。为 `1` 时不针对 `ionchamber records` 文件进行平均，而是把每次曝光记录单独处理。不设置(也就是不设为 `1` 时将根据 `ionchamber` 记录自动决定窗口大小，自动平均)。Buffer 不受这个选项控制。
    base_intensity    : 提供一个参考的 `ionchamber_intensity`，将所有的 Profile 缩放到这个光强上以提高可对比性。若放空的情况下将默认使用该次实验中的 buffer ionchamber 作为参考。(eg: 7e-10)
    buffer_ionchamber : 在部分情况下 `SourceFilePath` 中没有任何的 `buffer` 文件，脚本将自动从 `AveragedFilePath` 中寻找已有的 buffer profile，但由于没有 ionchamber 文件来说明该 buffer 的强度，因此从参数中提供该 buffer 的强度。

这些参数均可以修改，但不会自动添加默认值。

### 实验目录结构与具体操作过程 (Step by Step)

    ONE-ROOT-DIR
    |- Data
    |- Processed
    |- Averaged
    |- Subtracted
    |- GNOM
    |- config.yml
    |- setup.yml

在计算过程中，如果这些目录不存在会自动创建。但如果目标文件夹不是一个文件夹，会报错 `NotADirectoryError`。

在 `Data` 文件夹中将保存所有的原始数据。在 SSRF 的实验采集中，若使用自动曝光方式，将产生 `ionchamber` 文件 (`.Iochamber`)，记录文件 (`.log`) 以及相对应的 TIFF 图像文件 (`.tif`)。

例如一次连续曝光 5 次进行数据采集，将生成 5 张散射图像。

    beads_14_1_00001.log
    beads_14_1_1.Iochamber
    beads_14_1_00001.tif
    beads_14_1_00002.tif
    beads_14_1_00003.tif
    beads_14_1_00004.tif
    beads_14_1_00005.tif

### 模式一 (scale = 'ionchamber')：通过 log record 文件 Average, Scale and Subtract

`ionchamber` 文件记录了本次曝光整个过程中光强的变化，`.log` 文件，示例:

    t req.  t meas.     endTime               =========== File ============
    1.0000  1.0000  2018-05-04T13:50:50.752  /public/Data/2018/201805/20180504/EXP27/Data/beads_14_1_00001.tif
    1.0000  1.0000  2018-05-04T13:50:51.755  /public/Data/2018/201805/20180504/EXP27/Data/beads_14_1_00002.tif
    1.0000  1.0000  2018-05-04T13:50:52.758  /public/Data/2018/201805/20180504/EXP27/Data/beads_14_1_00003.tif
    1.0000  1.0000  2018-05-04T13:50:53.761  /public/Data/2018/201805/20180504/EXP27/Data/beads_14_1_00004.tif
    1.0000  1.0000  2018-05-04T13:50:53.761  /public/Data/2018/201805/20180504/EXP27/Data/beads_14_1_00005.tif
    --------- End ----------------------------

因此，自动脚本将找出在 `Data` 文件夹下所有的 `.log` 文件，并以 `.log` 文件中所记录的组别和对应的 `ionchamber` 文件，进行 `average` 操作。

                                          /beads_14_1_00001.dat
                                         | beads_14_1_00002.dat
    A_beads_14_1_00002 is averaged from -| beads_14_1_00003.dat
                                         | beads_14_1_00004.dat
                                          \beads_14_1_00005.dat
    Scale profile: A_beads_14_1_00002 *= 1.0558532832387264

在这个模式下，鉴于 `ionchamber` 是连续的，没有单独的记录，是先将组内的数据先平均，再对平均数据进行 `scale` 操作，最后减去 buffer，伪代码：

    scale_factor = current_ionchamber / first_group_ionchamber
    averaged_profile = radial_average(group_profiles) * scale_factor
    subtractd_profile = averaged_profile - buffer_profile

### 模式二 (scale = 'background') ：通过统计背景强度方式进行尾部对齐

通过算术平均的方式针对特定的 `q` 范围进行计算，得到缩放系数 (scaling factor)。

在 `config.yml` 中设定 `scale: background`，将通过散射背景强度进行对齐 (align)，指定 `scale: background` 时最好也指定 `window_size`, `scale_qmin` 和 `scale_qmin` 的参数。未设定这些参数将默认设为 `winddow_size: 5`, `scale_qmin: 0.23`, `scale_qmax: 0.26`。

## 暂时处理不了的情形 (TODO)

1. 一个组内如果有一帧异常，不想加到 average 的范围内，暂时没有很好的办法去处理。(可能的解决方法：加一个 `masked_profile` 的参数，每次做 average 的时候去检查包不包含不想要的内容)
2. 可以添加一个 `Option` 来输出中间过程。(用图画出平均前的 Profiles 和平均后的 Profile，比较方便 debug)
