from __future__ import print_function, division

import os
import sys
import platform
import glob
import subprocess

RAW_DIR = os.path.dirname(os.path.abspath(__file__))
if RAW_DIR not in sys.path:
    sys.path.append(RAW_DIR)
import SASExceptions


def findATSASDirectory():
    """Find ATSAS direcotry"""
    # from RAWOptions
    opsys = platform.system()

    if opsys == 'Darwin':
        default_path = '/Applications/ATSAS/bin'
    elif opsys == 'Windows':
        default_path = 'C:\\atsas\\bin'
    elif opsys == 'Linux':
        default_path = '~/atsas'
        default_path = os.path.expanduser(default_path)

        if os.path.exists(default_path):
            dirs = glob.glob(default_path + '/*')

            for item in dirs:
                if item.split('/')[-1].lower().startswith('atsas'):
                    default_path = item
                    break

            default_path = os.path.join(default_path, 'bin')

    is_atsas = os.path.exists(default_path)

    if is_atsas:
        return default_path

    if opsys != 'Windows':
        which = subprocess.Popen(['which', 'dammif'], stdout=subprocess.PIPE)
        output = which.communicate()

        atsas_path = output[0].strip()

        if atsas_path != '':
            return os.path.dirname(atsas_path)

    try:
        path = os.environ['PATH']
    except Exception:
        path = None

    if path != None:
        if opsys == 'Windows':
            split_path = path.split(';')
        else:
            split_path = path.split(':')

        for item in split_path:
            if item.lower().find('atsas') > -1 and item.lower().find(
                    'bin') > -1:
                return item

    try:
        atsas_path = os.environ['ATSAS']
    except Exception:
        atsas_path = None

    if atsas_path != None:
        if atsas_path.lower().find('atsas') > -1:
            atsas_path = atsas_path.rstrip('\\')
            atsas_path = atsas_path.rstrip('/')
            if atsas_path.endswith('bin'):
                return atsas_path
            else:
                if os.path.exists(os.path.join(atsas_path, 'bin')):
                    return os.path.join(atsas_path, 'bin')

    # return ''
    raise SASExceptions.NoATSASError('No ATSAS found')


class ErrorPrinter():
    def __init__(self, raw_settings, stdout):
        self._raw_settings = raw_settings
        self._stdout = stdout

    def showDataFormatError(self,
                            filename,
                            include_ascii=True,
                            include_sec=False):
        img_fmt = self._raw_settings.get('ImageFormat')

        if include_ascii:
            ascii_fmt = ' or any of the supported ASCII formats'
        else:
            ascii_fmt = ''

        if include_sec:
            sec = ' or the RAW SEC format'
        else:
            sec = ''

        print('Error loading file:', file=self._stdout)
        print(
            'The selected file: {}\n'
            'could not be recognized as a {} image format {} {}. '
            'This can be caused by failing to load the correct configuration file.\n'.
            format(filename, img_fmt, ascii_fmt, sec),
            'You can change the image format under Advanced Options in '
            'the Options menu.',
            file=self._stdout)

    def showSubtractionError(self, sasm, sub_sasm):
        filename1 = sasm.getParameter('filename')
        q1_min, q1_max = sasm.getQrange()
        points1 = len(sasm.i[q1_min:q1_max])
        filename2 = sub_sasm.getParameter('filename')
        q2_min, q2_max = sub_sasm.getQrange()
        points2 = len(sub_sasm.i[q2_min:q2_max])
        print(
            'Subtraction Error:\n',
            '  {} has {} data points.\n'.format(filename1, points1),
            '  {} has {} data points.\n'.format(filename2, points2),
            '  Subtraction is not possible. Data files must have equal number of points.',
            sep='\n',
            file=self._stdout,
        )

    def showAverageError(self, err_no):
        print('Average Error:', end=' ', file=self._stdout)
        if err_no == 1:
            print(
                'The selected items must have the same total number of points to be averaged.',
                file=self._stdout)
        elif err_no == 2:
            print(
                'Please select at least two items to be averaged.',
                file=self._stdout)
        elif err_no == 3:
            print(
                'The selected items must have the same q vectors to be averaged.',
                file=self._stdout)

    def showPleaseSelectItemsError(self, typename):
        print('No items selected Error:', end=' ', file=self._stdout)
        if typename == 'average':
            print(
                'Please select the items you want to average.',
                file=self._stdout)
        elif typename == 'subtract':
            print(
                'Please select the items you want the marked (star) item subtracted from.',
                file=self._stdout)
        elif typename == 'superimpose':
            print(
                'Please select the items you want the marked (star) item superimposed from.',
                file=self._stdout)
        elif typename == 'align':
            print(
                'Please select the items you want the marked (star) item aligned from.',
                file=self._stdout)
        elif typename == 'scale':
            print(
                'Please select the items you want the marked (star) item scaled from.',
                file=self._stdout)

    def showPleaseMarkItemError(self, typename):
        print('No item marked Error:', end=' ', file=self._stdout)
        if typename == 'subtract':
            print(
                'Please mark (star) the item you are using for subtraction',
                file=self._stdout)
        elif typename == 'merge':
            print(
                'Please mark (star) the item you are using as the main curve for merging',
                file=self._stdout)
        elif typename == 'superimpose':
            print(
                'Please mark (star) the item you want to superimpose to.',
                file=self._stdout)
        elif typename == 'interpolate':
            print(
                'Please mark (star) the item you are using as the main curve for interpolation',
                file=self._stdout)
        elif typename == 'align':
            print(
                'Please mark (star) the item you are using as the main curve for alignment',
                file=self._stdout)

    def showSaveError(self, err_type):
        if err_type == 'header':
            print(
                'Invalid Header Values: ',
                'Header values could not be saved, file was saved without them.',
                file=self._stdout)

    def showQvectorsNotEqualWarning(self, sasm, sub_sasm):
        sub_filename = sub_sasm.getParameter('filename')
        filename = sasm.getParameter('filename')

        print(
            'Q vectors do not match Warning:\n'
            '  Q vectors for {} and {} are not the same.\n'
            '  Continuing subtraction will attempt to find matching q regions in\n'
            '  or create matching q regions by binning.'.format(
                filename, sub_filename),
            file=self._stdout)
