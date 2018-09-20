from __future__ import print_function, division, absolute_import

import sys
from RAW import RAWSettings
from raw_script import get_max_length


def raw_cfg_printer(raw_cfg_path):
    raw_settings = RAWSettings.RawGuiSettings()
    success = RAWSettings.loadSettings(raw_settings, raw_cfg_path)
    if success:
        raw_settings.set('CurrentCfg', raw_cfg_path)
    else:
        raise IOError('Load failed, config file might be corrupted.')

    params = raw_settings.getAllParams()
    max_length = get_max_length(params.keys())

    for key in params.keys():
        print('{:<{l}}: {}'.format(key, raw_settings.get(key), l=max_length))


if __name__ == '__main__':
    cfg_fpath = sys.argv[1]
    raw_cfg_printer(cfg_fpath)
