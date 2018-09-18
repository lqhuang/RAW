# RAW Wrapper

Get tired of trivial GUI operations? Try these scripts to make life easier:)

English, [中文简介](./README_CN.md)

## Requirements

Minimal required packages (PY2/PY3):

    pip install numpy scipy pillow six pyyaml

Speed up profile calculation (image to curve) by using native c++ extension with python 2 interpreter in Linux/macOS system:

    pip install weave

Support extra formats of scattering image:

    pip install pyfai fabio

Quick install:

    pip install -r requirements.txt
    # or
    conda install --file requirements.txt

`conda env`:

    conda create -n raw --file requirement.txt
    source activate raw

## How to use?

Still writing ...... (Please see chinese version for more information.)

Run script:

    python raw_script.py /path/to/config.yml

Check log:

    less /path/of/run-root-directory/summary_timestamp.log

Append arguments from command line and **replace** those arguments from config file:

    python raw_script.py /path/to/config.yml --raw_cfg_path=/path/to/raw_settings.cfg --window_size=1

If we want to overwrite these parameters in configuration file, Add option `--write_to_file` to command line.:

    python raw_script.py --write_to_file /path/to/config.yml --raw_cfg_path=/path/to/raw_settings.cfg --window_size=1

Notice 1: Because of `pyyaml` packages, rewrite to file will lost all comments.

Notcie 2: Only parse `--key=arg` format for now.
