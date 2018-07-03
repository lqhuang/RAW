from __future__ import print_function, division, absolute_import

from io import open


def load_ionchamber(filepath, skip=0):
    intensity_list = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                intensity_list.append(float(line.split()[-1]))
    return intensity_list[skip:]


def get_mean_ionchamber(filepath, skip=0):
    ionchamber = load_ionchamber(filepath, skip=skip)
    return sum(ionchamber) / len(ionchamber)


def load_record(filepath):
    """read log file"""
    t_req = list()
    t_meas = list()
    end_time = list()
    file_list = list()
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if i > 0 and 'End' not in line:
                t0, t1, time, fp = line.split()
                t_req.append(t0)
                t_meas.append(t1)
                end_time.append(time)
                file_list.append(fp)
    return t_req, t_meas, end_time, file_list
