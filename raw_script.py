# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import

import os
import re
import sys
import glob
import argparse
from collections import Iterable
from datetime import datetime
from io import open  # python 2/3 compatibility

import yaml
from six import iteritems, string_types

from RAW import RAWSimulator
from fileio import get_mean_ionchamber, load_record


def remove_processed(data_list, processed_path):
    """Remove processed image data from given list"""
    processed_files = sorted(glob.glob1(processed_path, '*.dat'))
    processed_data = [os.path.splitext(fname)[0] for fname in processed_files]
    for filepath in reversed(data_list):
        fname = os.path.splitext(os.path.split(filepath)[1])[0]
        if fname in processed_data:
            data_list.remove(filepath)
    return data_list


def convert_ext(file_list, new_ext, new_dir=None):
    """System-based file extension converter"""
    new_file_list = list()
    for filepath in file_list:
        old_dir, fname = os.path.split(filepath)
        new_fname = '{}{}'.format(os.path.splitext(fname)[0], new_ext)
        if new_dir is None:
            new_dir = old_dir
        new_file_list.append(os.path.join(new_dir, new_fname))
    return new_file_list


def check_arguments():
    pass


def now():
    return u'[{}]'.format(str(datetime.now()))


def is_buffer(filename):
    return True if 'buffer' in filename.lower() else False


def split_buffer_sample(filenames):
    """Split file list into buffer file list and sample file list"""
    # Incredibily ingenious:) Just for fun. Readability is always more important.
    buffer_files, sample_files = [], []
    for fname in filenames:
        # (False, True)
        (sample_files, buffer_files)[is_buffer(fname)].append(fname)
    # equal to
    # buffer_files = [fname for fname in filenames if is_buffer(fname)]
    # sample_files = [fname for fname in filenames if not is_buffer(fname)]
    return buffer_files, sample_files


def strip_dir(filename):
    strip_func = lambda x: os.path.split(x)[-1]
    if isinstance(filename, string_types):
        return strip_func(filename)
    elif isinstance(filename, Iterable):
        return [strip_func(each) for each in filename]


def get_max_length(filenames):
    striped = map(strip_dir, filenames)
    length = map(len, striped)
    return max(length)


def average_for_one_record():
    pass


def run_RAW(exp_config):
    # TODO: complete this scripts

    raw_cfg_path = exp_config['raw_cfg_path']
    exp_root_path = exp_config['exp_root_path']

    source_data_path = os.path.join(exp_root_path, 'Data')
    num_skip = exp_config.get('skip_frames', 1)
    buffer_num_skip = exp_config.get('buffer_skip_frames', 1)
    num_frames_per_group = exp_config.get('window_size', 5)

    alignment = exp_config.get('scale', 'statistics')
    # alignment = 'statistics'
    # alignment = None
    scale_qmin = exp_config.get('scale_qmin', 0.23)
    sclae_qmax = exp_config.get('scale_qmax', 0.26)

    buffer_scaling_factor = exp_config.get('buffer_scaling_factor', 1.0)

    ProcessedFilePath = exp_config.get('ProcessedFilePath', 'Processed')
    AveragedFilePath = exp_config.get('AveragedFilePath', 'Averaged')
    SubtractedFilePath = exp_config.get('SubtractedFilePath', 'Subtracted')
    GnomFilePath = exp_config.get('GnomFilePath', 'GNOM')
    raw_settings = {
        'ProcessedFilePath': os.path.join(exp_root_path, ProcessedFilePath),
        'AveragedFilePath': os.path.join(exp_root_path, AveragedFilePath),
        'SubtractedFilePath': os.path.join(exp_root_path, SubtractedFilePath),
        'GnomFilePath': os.path.join(exp_root_path, GnomFilePath),
        'AutoSaveOnImageFiles': True,
        'AutoSaveOnSub': True,
        'AutoSaveOnAvgFiles': True,
        'AutoSaveOnGnom': False,
        'DatHeaderOnTop': True,
    }

    raw_simulator = RAWSimulator(
        raw_cfg_path,
        # exp_config['log_file'],
        do_analysis=False,
    )
    raw_simulator.set_raw_settings(**raw_settings)

    img_ext = 'tif'
    file_pattern = os.path.join(source_data_path, '*.' + img_ext)
    source_data_list = sorted(glob.glob(file_pattern))

    if not exp_config.get('overwrite', False):
        source_data_list = remove_processed(source_data_list,
                                            raw_settings['ProcessedFilePath'])

    if source_data_list:
        source_frames = raw_simulator.loadSASMs(source_data_list)
    else:
        file_pattern = os.path.join(raw_settings['ProcessedFilePath'], '*.dat')
        processed_files_list = sorted(glob.glob(file_pattern))
        source_frames = raw_simulator.loadSASMs(processed_files_list)

    buffer_frames = []
    sample_frames = []

    for each_sasm in source_frames:
        filename = each_sasm.getParameter('filename')
        # hard encoding, recognize 'buffer' string
        # or use intial_buffer_frame and final_buffer_frame to locate buffer?
        if 'buffer' in filename:
            buffer_frames.append(each_sasm)
        else:
            sample_frames.append(each_sasm)

    if buffer_frames and len(buffer_frames) > 1:
        average_buffer_sasm = raw_simulator.averageSASMs(
            buffer_frames[buffer_num_skip:])
        raw_simulator.saveSASM(
            average_buffer_sasm,
            '.dat',
            save_path=raw_settings['AveragedFilePath'],
        )
    elif buffer_frames and len(
            buffer_frames) == 1 and buffer_frames[0].startswith('A_'):
        # no buffer in `Processed` directory but A_buffer exist in `Processed`.
        # It's a backward compatibility for previous experiment data.
        # These raw data only saved average A_buffer_file in 'Valid_Frames'
        # or 'Frames' directory which has been renamed as `Processsed`
        # directory and save all buffer frames
        average_buffer_sasm = raw_simulator.loadSASMs(buffer_frames)[0]
    else:
        avg_buffer_pattern = os.path.join(raw_settings['AveragedFilePath'],
                                          '*buffer*.dat')
        avg_buffer_list = glob.glob(avg_buffer_pattern)
        if not avg_buffer_pattern:
            raise FileNotFoundError('No averaged buffer dat found.')
        elif len(avg_buffer_list) > 1:
            raise Warning(
                'Exist two or more buffer dats. The first one will be used.')
        else:
            average_buffer_sasm = raw_simulator.loadSASMs(avg_buffer_list)[0]

    raw_simulator.scaleSASMs([average_buffer_sasm], [buffer_scaling_factor])

    # TODO: save figures of middle process for debugging before alignment and after alignment
    if alignment == 'statistics':
        raw_simulator.alignSASMs(sample_frames[num_skip], sample_frames,
                                 (scale_qmin, sclae_qmax))
    else:
        # Not scale. Do nothing
        pass

    if num_frames_per_group > 1:
        sample_frames_by_group = [
            sample_frames[i:i + num_frames_per_group]
            for i in range(0, len(sample_frames), num_frames_per_group)
        ]
        average_sasm_list = [
            raw_simulator.averageSASMs(per_group[num_skip:])
            for per_group in sample_frames_by_group
        ]
    else:
        average_sasm_list = sample_frames

    subtracted_sasm_list = raw_simulator.subtractSASMs(average_buffer_sasm,
                                                       average_sasm_list)

    # for each_sasm in subtracted_sasm_list:
    #     raw_simulator.analyse(each_sasm)


def run_automatic(exp_config, log_file=sys.stdout):
    # setup configuration
    raw_cfg_path = exp_config['raw_cfg_path']
    exp_root_path = exp_config['exp_root_path']

    source_data_path = os.path.join(exp_root_path, 'Data')
    num_skip = exp_config.get('skip_frames', 1)
    buffer_num_skip = exp_config.get('buffer_skip_frames', 1)
    num_frames_per_group = exp_config.get('window_size', 5)

    ProcessedFilePath = exp_config.get('ProcessedFilePath', 'Processed')
    AveragedFilePath = exp_config.get('AveragedFilePath', 'Averaged')
    SubtractedFilePath = exp_config.get('SubtractedFilePath', 'Subtracted')
    GnomFilePath = exp_config.get('GnomFilePath', 'GNOM')
    raw_settings = {
        'ProcessedFilePath': os.path.join(exp_root_path, ProcessedFilePath),
        'AveragedFilePath': os.path.join(exp_root_path, AveragedFilePath),
        'SubtractedFilePath': os.path.join(exp_root_path, SubtractedFilePath),
        'GnomFilePath': os.path.join(exp_root_path, GnomFilePath),
        'AutoSaveOnImageFiles': True,
        'AutoSaveOnSub': True,
        'AutoSaveOnAvgFiles': True,
        'AutoSaveOnGnom': False,
        'DatHeaderOnTop': True,
    }

    intensity_base = 2.509468563829787e-10

    raw_simulator = RAWSimulator(
        raw_cfg_path,
        # exp_config['log_file'],
        # log_file=log_file,
        do_analysis=False,
    )
    raw_simulator.set_raw_settings(**raw_settings)

    img_ext = '.tif'
    dat_ext = '.dat'
    ionchamber_ext = '.Iochamber'
    record_ext = '.log'

    # generate file map
    img_pattern = os.path.join(source_data_path, '*%s' % img_ext)
    ionchamber_pattern = os.path.join(source_data_path, '*%s' % ionchamber_ext)
    record_pattern = os.path.join(source_data_path, '*%s' % record_ext)

    scatter_image_files = glob.glob(img_pattern)
    ionchamber_files = glob.glob(ionchamber_pattern)
    record_files = glob.glob(record_pattern)

    # Stage ?: convert scattering image to scattering profile

    overwirte_processed = exp_config.get('overwrite', False)
    num_total_images = len(scatter_image_files)
    num_processed = len(
        glob.glob(os.path.join(raw_settings['ProcessedFilePath'], '*.dat')))
    num_unprocessed = num_total_images - num_processed
    output = u"{} Find {} scattering images, processed {}, unprocessed {}. 'overwrite_processed={}`.".format(
        now(), num_total_images, num_processed, num_unprocessed,
        overwirte_processed)
    print(output, file=log_file)
    if not overwirte_processed:
        unprocessed_files = remove_processed(scatter_image_files,
                                             raw_settings['ProcessedFilePath'])
    else:
        unprocessed_files = scatter_image_files
    if unprocessed_files:
        raw_simulator.loadSASMs(unprocessed_files)
        for each in unprocessed_files:
            fname = os.path.splitext(strip_dir(each))[0]
            print(
                u'  {0}{ext1} -> {0}{ext2}'.format(
                    fname, ext1=img_ext, ext2=dat_ext),
                file=log_file)
    print(u'', file=log_file)

    # Averaged Strategy
    # 1. by-window-size (average all buffer, then average each group by window-size)
    # 2. by-log-record (extract log record, average each group in record)

    # processed_pattern = os.path.join(raw_settings['ProcessedFilePath'], '*%s' % dat_ext) # yapf: disable
    # processed_files = glob.glob(processed_pattern)
    # processed_sasms = raw_simulator.loadSASMs(processed_files)

    # Stage ?: Average processed profile

    buffer_ionchamber, sample_ionchamber = split_buffer_sample(
        ionchamber_files)
    buffer_records, sample_records = split_buffer_sample(record_files)

    # sort filename list with last integer (eg: ***_1.ext, ***_002.ext)
    regpattern = re.compile(r'\d+')
    key_func = lambda x: int(regpattern.findall(os.path.split(x)[-1])[-1])
    buffer_records.sort(key=key_func)
    sample_records.sort(key=key_func)
    buffer_ionchamber.sort(key=key_func)
    sample_ionchamber.sort(key=key_func)

    ## Buffer
    for i, (rec, ion) in enumerate(zip(buffer_records, buffer_ionchamber)):
        print(
            now(),
            strip_dir(rec.decode('utf-8')),
            strip_dir(ion.decode('utf-8')),
            file=log_file)
        _, _, _, filegroup = load_record(rec)

        mean_intensity = get_mean_ionchamber(ion)
        scale_factor = intensity_base / mean_intensity

        proc_dat = convert_ext(filegroup, dat_ext, new_dir=raw_settings['ProcessedFilePath'])
        loaded_sasms = raw_simulator.loadSASMs(proc_dat)
        averaged_sasm = raw_simulator.averageSASMs(loaded_sasms[num_skip:])
        print(u'  ' + averaged_sasm.getParameter('filename').decode('utf-8'), file=log_file)
        # print(u'  averaged from ', u'  \n'.join(strip_dir(filegroup)), file=log_file)
        print(u'', file=log_file)
        # 如果有很多个 buffer 的记录，排序第一个的 buffer 将会被使用
        if i == 0:
            averaged_buffer_sasm = averaged_sasm

    ## 没有 buffer 的话，直接去 Averaged 里面找带 buffer 的第一个 .dat 文件
    if not buffer_records:
        buffer_file = glob.glob(os.path.join(raw_settings['AveragedFilePath'], '*buffer*%s' % dat_ext))[0]
        averaged_buffer_sasm = raw_simulator.loadSASMs([buffer_file])[0]

        # 这个 buffer 的 intensity 要怎么处理
        scale_factor = intensity_base / exp_config.get('buffer_ionchamber', intensity_base)
        averaged_buffer_sasm.scale(scale_factor)

    ## Sample
    for rec, ion in zip(sample_records, sample_ionchamber):
        print(
            now(),
            strip_dir(rec.decode('utf-8')),
            strip_dir(ion.decode('utf-8')),
            file=log_file)
        _, _, _, filegroup = load_record(rec)

        mean_intensity = get_mean_ionchamber(ion)
        scale_factor = intensity_base / mean_intensity

        proc_dat = convert_ext(filegroup, dat_ext, new_dir=raw_settings['ProcessedFilePath'])
        loaded_sasms = raw_simulator.loadSASMs(proc_dat)
        if exp_config['window_size'] != 1:
            averaged_sasm = raw_simulator.averageSASMs(loaded_sasms[num_skip:])
            averaged_sasm.scale(scale_factor)
            print(u'  ' + averaged_sasm.getParameter('filename').decode('utf-8'), file=log_file)
            # print(u'  averaged from ', u'  \n'.join(strip_dir(filegroup)), file=log_file)

            subtracted_sasm = raw_simulator.subtractSASMs(averaged_buffer_sasm, [averaged_sasm])[0]
            print(u'  ' + subtracted_sasm.getParameter('filename').decode('utf-8'), file=log_file)
            print(u'', file=log_file)
        else:
            subtracted_sasms = raw_simulator.subtractSASMs(averaged_buffer_sasm, loaded_sasms)

    max_length = get_max_length(ionchamber_files)
    for ion_file in ionchamber_files:
        mean_intensity = get_mean_ionchamber(ion_file)
        scale_factor = intensity_base / mean_intensity
        output = u'{:<{l}}: {} {}'.format(
            strip_dir(ion_file), mean_intensity, scale_factor, l=max_length)
        print(output, file=log_file)

    # print(u'Summary:', file=log_file)
    # print(u'Number of source data', file=log_file)


def main():

    config_file = os.path.realpath(sys.argv[1])
    with open(config_file, 'r', encoding='utf-8') as fstream:
        exp_config = yaml.load(fstream)

    # overwrite with input argument

    # check essential arguments
    if exp_config.get('exp_root_path', None) is None:
        exp_config['exp_root_path'] = os.path.dirname(config_file)
        print(exp_config['exp_root_path'])
    if exp_config.get('raw_cfg_path', None) is None:
        raise ValueError('please specify a cfg file.')

    # processing_summary.txt
    summary_file = os.path.join(exp_config['exp_root_path'], 'summary.txt')
    with open(summary_file, 'w', encoding='utf-8') as summary:
        # print configuration from config file.
        print(
            u'{} Current configuration from file'.format(now()), file=summary)
        max_length = get_max_length(exp_config.keys())
        for key, val in iteritems(exp_config):
            print(
                u'  {:<{l}}: {}'.format(key, val, l=max_length), file=summary)
        print(u'', file=summary)
        print(
            u'{} Appended arguments from command line:'.format(now()),
            file=summary)
        print(u'', file=summary)
        print(
            u'{} Check and set default value for missing arguments:'.format(
                now()),
            file=summary)
        print(u'', file=summary)
        run_automatic(exp_config, summary)

    # run_RAW(exp_config)


if __name__ == '__main__':
    main()
