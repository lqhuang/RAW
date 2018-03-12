from __future__ import print_function, division

import os
import sys
import glob

from RAWWrapper import RAWSimulator


def remove_processed(data_list, processed_path):
    """Remove processed image data from given list"""
    processed_files = sorted(glob.glob1(processed_path, '*.dat'))
    processed_data = [os.path.splitext(fname)[0] for fname in processed_files]
    for filepath in reversed(data_list):
        fname = os.path.splitext(os.path.split(filepath)[1])[0]
        if fname in processed_data:
            data_list.remove(filepath)
    return data_list


def main():
    # TODO: complete this scripts
    exp_root_path = sys.argv[0]

    exp_config = {
        'raw_cfg_path': sys.argv[1],
        'log_file': os.path.join(exp_root_path, 'log.txt'),
        'source_data_path': os.path.join(exp_root_path, 'Data'),
        'overwrite': False,
        'SKIP_FRAMES': 1,
        'WINDOW_SIZE': 5,
        'SCALE_QMIN': 0.23,
        'SCALE_QMAX': 0.26,
    }

    num_skip = exp_config['SKIP_FRAMES']
    num_frames_per_group = exp_config['WINDOW_SIZE']
    scale_qmin = exp_config['SCALE_QMIN']
    sclae_qmax = exp_config['SCALE_QMAX']

    raw_settings = {
        'ProcessedFilePath': os.path.join(exp_root_path, 'Processed'),
        'AveragedFilePath': os.path.join(exp_root_path, 'Averaged'),
        'SubtractedFilePath': os.path.join(exp_root_path, 'Subtracted'),
        'GnomFilePath': os.path.join(exp_root_path, 'GNOM'),
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
    if not os.path.exists(raw_settings['GnomFilePath']):
        os.makedirs(raw_settings['GnomFilePath'])

    raw_simulator = RAWSimulator(
        exp_config['raw_cfg_path'],
        # exp_config['log_file'],
        do_analysis=True,
    )
    raw_simulator.set_raw_settings(**raw_settings)

    img_ext = 'tif'
    file_pattern = os.path.join(exp_config['source_data_path'], '*.' + img_ext)
    source_data_list = sorted(glob.glob(file_pattern))

    if not exp_config['overwrite']:
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

    average_buffer_sasm = raw_simulator.averageSASMs(buffer_frames[num_skip:])

    # TODO: save figures of middle process for debugging
    # before alignment and after alignment
    raw_simulator.alignSASMs(sample_frames[0], sample_frames,
                             (scale_qmin, sclae_qmax))

    sample_frames_by_group = [
        sample_frames[i:i + num_frames_per_group]
        for i in range(0, len(sample_frames), num_frames_per_group)
    ]
    average_sasm_list = [
        raw_simulator.averageSASMs(per_group[num_skip:])
        for per_group in sample_frames_by_group
    ]

    subtracted_sasm_list = raw_simulator.subtractSASMs(average_buffer_sasm,
                                                       average_sasm_list)

    for each_sasm in subtracted_sasm_list:
        raw_simulator.analyse(each_sasm)


if __name__ == '__main__':
    main()
