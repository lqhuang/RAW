from __future__ import print_function, division, absolute_import

import os
import sys
import glob
import argparse

from io import open  # python 2/3 compatibility
import yaml

from RAW import RAWSimulator


def remove_processed(data_list, processed_path):
    """Remove processed image data from given list"""
    processed_files = sorted(glob.glob1(processed_path, '*.dat'))
    processed_data = [os.path.splitext(fname)[0] for fname in processed_files]
    for filepath in reversed(data_list):
        fname = os.path.splitext(os.path.split(filepath)[1])[0]
        if fname in processed_data:
            data_list.remove(filepath)
    return data_list


def run_RAW(exp_config):
    # TODO: complete this scripts

    # exp_config = {
    #     'raw_cfg_path': sys.argv[1],
    #     'log_file': os.path.join(exp_root_path, 'log.txt'),
    #     'source_data_path': os.path.join(exp_root_path, 'Data'),
    #     'overwrite': False,
    #     'SKIP_FRAMES': 1,
    #     'WINDOW_SIZE': 5,
    #     'SCALE_QMIN': 0.23,
    #     'SCALE_QMAX': 0.26,
    # }

    raw_cfg_path = exp_config.get('raw_cfg_path', None)
    if raw_cfg_path is None:
        raise ValueError('please specify a raw cfg file.')

    exp_root_path = exp_config.get('exp_root_path', None)

    source_data_path = os.path.join(exp_root_path, 'Data')
    num_skip = exp_config.get('skip_frames', 1)
    buffer_num_skip = exp_config.get('buffer_skip_frames', 1)
    num_frames_per_group = exp_config.get('window_size', 5)

    alignment = exp_config.get('scale', 'statistics')
    # alignment = 'statistics'
    # alignment = 'ionchamber'
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

    # Init directories
    if not os.path.exists(raw_settings['ProcessedFilePath']):
        os.makedirs(raw_settings['ProcessedFilePath'])
    if not os.path.exists(raw_settings['AveragedFilePath']):
        os.makedirs(raw_settings['AveragedFilePath'])
    if not os.path.exists(raw_settings['SubtractedFilePath']):
        os.makedirs(raw_settings['SubtractedFilePath'])
    # if not os.path.exists(raw_settings['GnomFilePath']):
    #     os.makedirs(raw_settings['GnomFilePath'])

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
    elif alignment == 'ionchamber':
        ionchamber_pattern = os.path.join(source_data_path, '*.[Ii]o*chamber')
        ionchamber_list = sorted(glob.glob(ionchamber_pattern))
        for ion_name in reversed(ionchamber_list):
            if 'buffer' in os.path.basename(ion_name):
                ionchamber_list.remove(ion_name)
        ionchamber_factors = list()
        for ion_name in ionchamber_list:
            # TODO: how to decide number of points to skip in ionchamber?
            # sometimes, bad intensity log is negative, but sometimes not.
            # TODO: match ionchamber file to id list of frames
            ion_intensity = dat.load_ionchamber(ion_name, skip=2)
            factor = sum(ion_intensity) / len(ion_intensity)
            ionchamber_factors.append(factor)
        scaling_factors = [
            each / ionchamber_factors[0] for each in ionchamber_factors
        ]
        print('\n'.join(str(each) for each in scaling_factors))
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

    if alignment == 'ionchamber':
        raw_simulator.scaleSASMs(average_sasm_list, scaling_factors)

    subtracted_sasm_list = raw_simulator.subtractSASMs(average_buffer_sasm,
                                                       average_sasm_list)

    # for each_sasm in subtracted_sasm_list:
    #     raw_simulator.analyse(each_sasm)


def main():

    config_file = os.path.realpath(sys.argv[1])
    with open(config_file, 'r', encoding='utf-8') as fstream:
        exp_config = yaml.load(fstream)

    if exp_config.get('exp_root_path', None) is None:
        exp_config['exp_root_path'] = os.path.dirname(config_file)
        print(exp_config['exp_root_path'])
    run_RAW(exp_config)


if __name__ == '__main__':
    main()
