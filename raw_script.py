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

PY2 = sys.version_info < (3, 0)

def unicodify(text):
    if PY2 and isinstance(text, (bytes, str)):
        return text.decode('utf-8')
    else:
        return text


def remove_processed(data_list, processed_path, target_ext='.dat'):
    """Remove processed image data from given list"""
    processed_files = sorted(glob.glob1(processed_path, '*%s' % target_ext))
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


def parse_args(args):
    # args (list):
    args_dict = {}
    for item in args:
        try:
            key, arg = unicodify(item).strip('-').split('=')
            args_dict[key] = arg
        except ValueError:
            # ValueError: not enough values to unpack (expected 2, got 1)
            raise ValueError('Error: Unrecognized argument {}. Please retry with --key=val format.'.format(item))
    return args_dict


def now():
    return u'[{}]'.format(str(datetime.now()))


def now_under_score():
    return u'{}'.format(
        datetime.strftime(datetime.now(), r'%Y_%m_%d_%H_%M_%S'))


def is_buffer(filename):
    return True if 'buffer' in filename.lower() else False


def split_buffer_sample(filenames):
    """Split file list into buffer file list and sample file list"""
    # Incredibily ingenious:) Just for fun.
    # Readability is always more important.
    buffer_files, sample_files = [], []
    for fname in filenames:
        # (False, True)
        (sample_files, buffer_files)[is_buffer(fname)].append(fname)
    # equal to
    # buffer_files = [fname for fname in filenames if is_buffer(fname)]
    # sample_files = [fname for fname in filenames if not is_buffer(fname)]
    return buffer_files, sample_files


def strip_dir(filename):
    def strip_func(x):
        return os.path.split(x)[-1]

    if isinstance(filename, string_types):
        return strip_func(filename)
    elif isinstance(filename, Iterable):
        return [strip_func(each) for each in filename]


def get_max_length(filenames):
    striped = map(strip_dir, filenames)
    length = map(len, striped)
    return max(length)


def gen_average_form(target, filelist, prepadding=2):
    filelist = strip_dir(filelist)
    target_leng = len(target)
    const_leng = 19
    size = len(filelist)
    target_line = u' ' * prepadding + target + u' is averaged from -'

    if size < 2:
        form = target_line + filelist[0]
    else:
        brace = [u' /'] + [u'| '] * (size - 2) + [u' \\']
        blank = [u' ' * (prepadding + target_leng + const_leng)]

        prefix_header = blank * ((size - 1) // 2) + \
            [target_line] + blank * (size // 2)
        form = u'\n'.join(
            map(lambda x: u''.join(x), zip(prefix_header, brace, filelist)))
    # print(form)
    return form


def check_essential_arguments(exp_config, config_file):
    missing_args = []

    check_info = u'{} Check and set default value for missing arguments:\n'.format(now())
    if exp_config.get('exp_root_path', None) is None:
        missing_args.append('exp_root_path')
        exp_config['exp_root_path'] = os.path.abspath(os.path.dirname(config_file))
    if exp_config.get('raw_cfg_path', None) is None:
        check_info += u'  Error: No raw_cfg_path found.'
        missing_args.append('raw_cfg_path')
        raise ValueError('please specify a cfg file.')

    default_options = {
        'skip_frames': 1,
        'buffer_skip_frames': 1,
        'SourceFilePath': 'Data',
        'ProcessedFilePath': 'Processed',
        'AveragedFilePath': 'Averaged',
        'SubtractedFilePath': 'Subtracted',
        'GnomFilePath': 'GNOM',
        'img_ext': '.tif',
        'dat_ext': '.dat',
        'ionchamber_ext': '.Iochamber',
        'record_ext': '.log',
        'overwrite': False,
        'scale': 'ionchamber',
    }
    existent_args = exp_config.keys()

    for key, val in iteritems(default_options):
        if key not in existent_args:
            missing_args.append(key)
            exp_config[key] = val

    max_length = get_max_length(missing_args)
    for key in missing_args:
        val = exp_config[key]
        check_info += u'  {:<{l}}: {}\n'.format(key, val, l=max_length)

    # check scale_qmin, scale_qmax
    if ('scale' in exp_config and exp_config['scale'] == 'background') or \
       ('scale' in missing_args):
        if 'scale_qmin' not in exp_config:
            exp_config['scale_qmin'] = 0.23
            check_info += u'  {:<{l}}: {}\n'.format('scale_qmin', 0.23, l=max_length)
        if 'scale_qmax' not in exp_config:
            exp_config['scale_qmax'] = 0.26
            check_info += u'  {:<{l}}: {}\n'.format('scale_qmax', 0.26, l=max_length)
        if 'window_size' not in exp_config:
            exp_config['window_size'] = 5
            check_info += u'  {:<{l}}: {}\n'.format('window_size', 5, l=max_length)

    check_info += '\n'
    return check_info


def run_RAW(exp_config):
    """Old version"""
    # TODO: complete this scriptss
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

    num_skip = exp_config['skip_frames']
    buffer_num_skip = exp_config['buffer_skip_frames']

    SourceFilePath = exp_config['SourceFilePath']
    ProcessedFilePath = exp_config['ProcessedFilePath']
    AveragedFilePath = exp_config['AveragedFilePath']
    SubtractedFilePath = exp_config['SubtractedFilePath']
    GnomFilePath = exp_config['GnomFilePath']
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
    source_data_path = os.path.join(exp_root_path, SourceFilePath)

    img_ext = exp_config['img_ext']  # '.tif'
    dat_ext = exp_config['dat_ext']  # '.dat'
    ionchamber_ext = exp_config['ionchamber_ext']  # '.Iochamber'
    record_ext = exp_config['record_ext']  # '.log'

    overwirte_processed = exp_config['overwrite']

    base_intensity = exp_config.get('base_ionchamber', None)

    alignment = exp_config.get('scale', None)

    raw_simulator = RAWSimulator(
        raw_cfg_path,
        # exp_config['log_file'],
        # log_file=log_file,
        do_analysis=False,
    )
    raw_simulator.set_raw_settings(**raw_settings)

    # ================= generate file map =========================== #
    img_pattern = os.path.join(source_data_path, '*%s' % img_ext)
    ionchamber_pattern = os.path.join(source_data_path, '*%s' % ionchamber_ext)
    record_pattern = os.path.join(source_data_path, '*%s' % record_ext)

    scatter_image_files = glob.glob(img_pattern)
    ionchamber_files = glob.glob(ionchamber_pattern)
    record_files = glob.glob(record_pattern)

    buffer_ionchamber, sample_ionchamber = split_buffer_sample(
        ionchamber_files)
    buffer_records, sample_records = split_buffer_sample(record_files)

    # sort filename list with last integer (eg: ***_1.ext, ***_002.ext)
    regpattern = re.compile(r'\d+')

    def key_func(x):
        return int(regpattern.findall(os.path.split(x)[-1])[-1])

    ionchamber_files.sort()
    buffer_records.sort(key=key_func)
    sample_records.sort(key=key_func)
    buffer_ionchamber.sort(key=key_func)
    sample_ionchamber.sort(key=key_func)

    # ===== convert scattering image to scattering profile ========== #
    num_total_images = len(scatter_image_files)
    processed_files = sorted(
        glob.glob1(raw_settings['ProcessedFilePath'], '*%s' % dat_ext))
    num_processed = len(processed_files)
    num_unprocessed = num_total_images - num_processed
    output = u"{} Found {} scattering images, processed {}, unprocessed {}. 'overwrite_processed={}`.".format(
        now(), num_total_images, num_processed, num_unprocessed,
        overwirte_processed)
    print(output, file=log_file)
    if not overwirte_processed:
        unprocessed_files = remove_processed(
            scatter_image_files, raw_settings['ProcessedFilePath'], dat_ext)
    else:
        unprocessed_files = scatter_image_files
    if processed_files:
        print('  List of processed profiles:', file=log_file)
        print('\n'.join(map(lambda x: '  %s' % x, processed_files)),
            file=log_file)
    if unprocessed_files:
        raw_simulator.loadSASMs(unprocessed_files)
        for each in unprocessed_files:
            fname = os.path.splitext(strip_dir(each))[0]
            print(
                u'  {0}{ext1} -> {0}{ext2}'.format(
                    fname, ext1=img_ext, ext2=dat_ext),
                file=log_file)
    print(u'', file=log_file)

    # ================ Show ionchamber intensity ==================== #
    print(u"{} Found {} ionchamber files, ".format(now(), len(ionchamber_files)),
        file=log_file)
    if base_intensity is None:
        print(u'  No base ionchamber intensity has been set. '
            'The last ionchamber file will be used as base intensity.',
            file=log_file)
        if buffer_ionchamber:
            base_intensity = get_mean_ionchamber(buffer_ionchamber[-1])
        else:
            if exp_config.get('buffer_ionchamber', None) is None:
                print(u'  Error: `buffer_ionchamber` argument is missing.', file=log_file)
                raise ValueError('`buffer_ionchamber` argument is missing.')
            print(u'  Error: no ionchamber file found.', file=log_file)
            raise FileNotFoundError(
                'At least one buffer ionchamber file is required while no base intensity provided.'
            )

    print(u'  Current base intensity: {}'.format(base_intensity), file=log_file)
    print(u'  Display scaling factor for each record', file=log_file)
    max_length = get_max_length(ionchamber_files)
    print(u'  {:<{l}} {:>13} {:>14}'.format('Filename', 'Ion Intensity', 'Scaling Factor', l=max_length + 1),
        file=log_file)
    for ion_file in ionchamber_files:
        mean_intensity = get_mean_ionchamber(ion_file)
        scale_factor = base_intensity / mean_intensity
        output = u'  {:<{l}}: {:13.5e} {:14.5f}'.format(
            strip_dir(ion_file), mean_intensity, scale_factor, l=max_length)
        print(output, file=log_file)
    print('', file=log_file)

    # Averaged Strategy
    # 1. by-window-size (average all buffer, then average each group by window-size)
    # 2. by-log-record (extract log record, average each group in record)

    # processed_pattern = os.path.join(raw_settings['ProcessedFilePath'], '*%s' % dat_ext) # yapf: disable
    # processed_files = glob.glob(processed_pattern)
    # processed_sasms = raw_simulator.loadSASMs(processed_files)

    # Stage ?: Average processed profile

    # ================== Buffer ===================================== #
    for i, (rec, ion) in enumerate(zip(buffer_records, buffer_ionchamber)):
        print(
            now(),
            unicodify(strip_dir(rec)),
            unicodify(strip_dir(ion)),
            file=log_file)
        _, _, _, filegroup = load_record(rec)

        mean_intensity = get_mean_ionchamber(ion)
        scale_factor = base_intensity / mean_intensity

        proc_dat = convert_ext(
            filegroup, dat_ext, new_dir=raw_settings['ProcessedFilePath'])
        loaded_sasms = raw_simulator.loadSASMs(proc_dat)
        averaged_sasm = raw_simulator.averageSASMs(
            loaded_sasms[buffer_num_skip:])

        curr_filename = unicodify(averaged_sasm.getParameter('filename'))
        print(gen_average_form(curr_filename, proc_dat), file=log_file)
        # if multiple buffer records exist, the last record will be used as average buffer.
        if i == len(buffer_records) - 1:
            print(u'  *** this curve will be used as average buffer profile.',
                file=log_file)
            averaged_buffer_sasm = averaged_sasm

    # 没有 buffer 的话，直接去 Averaged 里面找带 buffer 的第一个 .dat 文件
    if not buffer_records:
        print(u'%s Warning: no buffer file in source data folder. Try to find first buffer from AveragedFilePath:' % now(),
            file=log_file)
        buffer_files = glob.glob(os.path.join(
            raw_settings['AveragedFilePath'], '*buffer*%s' % dat_ext))
        if not buffer_files:
            print(u'  Error: no buffer profile in AveragedFilePath', file=log_file)
            raise FileNotFoundError('No buffer profile in AveragedFilePath')
        averaged_buffer_sasm = raw_simulator.loadSASMs([buffer_files[0]])[0]
        curr_filename = unicodify(averaged_buffer_sasm.getParameter('filename'))
        print(u'  Found {} as buffer profile'.format(curr_filename), file=log_file)
        # 这个 buffer 的 intensity 要怎么处理?
        buffer_ion_intensity = exp_config.get('buffer_ionchamber', None)
        if buffer_ion_intensity is None:
            print(u'  Loading buffer profile from average profile but no buffer ionchamber set in config',
                file=log_file)
            raise ValueError('Please set a ionchamber intensity for buffer.')
        print(u'  Get buffer ionchamber from config -> buffer_ionchamber: {}'.format(buffer_ion_intensity),
            file=log_file)
        scale_factor = (base_intensity / buffer_ion_intensity)
        print(u' scale buffer profile with base_intensity: {} *= {}'.format(curr_filename, scale_factor),
            file=log_file)
        averaged_buffer_sasm.scale(scale_factor)

    averaged_buffer_fname = unicodify(
        averaged_buffer_sasm.getParameter('filename'))
    print(u'', file=log_file)

    if alignment == 'ionchamber' or alignment is None:
        # ================== Sample ===================================== #
        for rec, ion in zip(sample_records, sample_ionchamber):
            print(
                now(),
                unicodify(strip_dir(rec)),
                unicodify(strip_dir(ion)),
                file=log_file)
            _, _, _, filegroup = load_record(rec)

            mean_intensity = get_mean_ionchamber(ion)
            scale_factor = base_intensity / mean_intensity

            proc_dat = convert_ext(
                filegroup, dat_ext, new_dir=raw_settings['ProcessedFilePath'])
            loaded_sasms = raw_simulator.loadSASMs(proc_dat)

            if exp_config['window_size'] != 1:
                averaged_sasm = raw_simulator.averageSASMs(loaded_sasms[num_skip:])
                averaged_fname = unicodify(averaged_sasm.getParameter('filename'))
                print(gen_average_form(averaged_fname, proc_dat), file=log_file)

                averaged_sasm.scale(scale_factor)
                print(u'  Scale profile: {} *= {}'.format(averaged_fname, scale_factor),
                    file=log_file)

                subtracted_sasm = raw_simulator.subtractSASMs(
                    averaged_buffer_sasm, [averaged_sasm])[0]
                subtracted_fname = unicodify(
                    subtracted_sasm.getParameter('filename'))
                print(u'  Subtract: {} = {} - {}'.format(subtracted_fname, averaged_fname, averaged_fname),
                    file=log_file)
                print(u'', file=log_file)
            else:  # window_size == 1
                raw_simulator.scaleSASMs(
                    loaded_sasms, [scale_factor] * len(loaded_sasms))
                print(u'  Scale all loaded profiles by scaling_factor:', scale_factor,
                    file=log_file)
                print(u'  Then do', file=log_file)
                subtracted_sasms = raw_simulator.subtractSASMs(
                    averaged_buffer_sasm, loaded_sasms)
                for sub_sasm, sasm in zip(subtracted_sasms, loaded_sasms):
                    print(u'  Subtract: {} = {} - {}'.format(
                        unicodify(strip_dir(sub_sasm.getParameter('filename'))),
                        unicodify(strip_dir(sasm.getParameter('filename'))),
                        unicodify(strip_dir(averaged_buffer_fname)),
                    ), file=log_file)

    elif alignment == 'background':
        scale_qmin = exp_config['scale_qmin']
        scale_qmax = exp_config['scale_qmax']
        window_size = exp_config['window_size']

        file_pattern = os.path.join(raw_settings['ProcessedFilePath'], '*%s' % dat_ext)
        processed_files_list = sorted(glob.glob(file_pattern))
        source_frames = raw_simulator.loadSASMs(processed_files_list)

        sample_frames = []

        for each_sasm in source_frames:
            filename = each_sasm.getParameter('filename')
            if 'buffer' not in filename:
                sample_frames.append(each_sasm)

        raw_simulator.alignSASMs(sample_frames[num_skip], sample_frames, (scale_qmin, scale_qmax))

        if window_size > 1:
            sample_frames_by_group = [
                sample_frames[i:i + window_size]
                for i in range(0, len(sample_frames), window_size)
            ]
            average_sasm_list = [
                raw_simulator.averageSASMs(per_group[num_skip:])
                for per_group in sample_frames_by_group
            ]
        else:
            average_sasm_list = sample_frames

        subtracted_sasm_list = raw_simulator.subtractSASMs(
            averaged_buffer_sasm, average_sasm_list)

    # print(u'Summary:', file=log_file)
    # print(u'Number of source data', file=log_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help="Defined configuration.")
    parser.add_argument('--write_to_file', action='store_true', help='Write appended arguments to configuration file.')
    parser.add_argument('...', nargs=argparse.REMAINDER, help='Other supplementary arguments in --key=val format.')
    args = parser.parse_args()

    config_file = args.config_file
    write_to_file = args.write_to_file
    rest_args = getattr(args, '...')
    # Arguments in this two positions are different.
    # python raw_script.py .\config.yml --write_to_file
    # python raw_script.py --write_to_file .\config.yml
    if '--write_to_file' in rest_args:
        write_to_file = True
        rest_args.remove('--write_to_file')
    rest_config = parse_args(rest_args)
    if 'write_to_file' in rest_config:
        rest_config.pop('write_to_file')

    # load configuration from config file.
    with open(config_file, 'r', encoding='utf-8') as fstream:
        exp_config = yaml.load(fstream)

    # print configuration from config file.
    cf_info = u'{} Current configuration from file\n'.format(now())
    max_length = get_max_length(exp_config.keys())
    for key, val in iteritems(exp_config):
        cf_info += u'  {:<{l}}: {}\n'.format(key, val, l=max_length)
    cf_info += '\n'

    # overwrite with input argument
    appending_info = u'{} Appended arguments from command line:\n'.format(now())
    if rest_config:
        max_length = get_max_length(rest_config.keys())
        for key, val in iteritems(rest_config):
            appending_info += u'  {:<{l}}: {}\n'.format(key, val, l=max_length)
    appending_info += u'\n'
    exp_config.update(rest_config)
    if rest_config and write_to_file:
        with open(config_file, 'w', encoding='utf-8') as fstream:
            yaml.dump(exp_config, fstream)

    # check essential arguments
    # Though `exp_config` isn't in return, this will modify its content/value.
    check_info = check_essential_arguments(exp_config, config_file)

    # processing_summary.log
    summary_file = os.path.join(
        exp_config['exp_root_path'], 'summary_{}.log'.format(now_under_score()))
    with open(summary_file, 'w', encoding='utf-8') as summary:
        summary = open(summary_file, 'w', encoding='utf-8')

        # write log
        summary.write(cf_info)
        summary.write(appending_info)
        summary.write(check_info)

        # run scripts
        run_automatic(exp_config, summary)


if __name__ == '__main__':
    main()
