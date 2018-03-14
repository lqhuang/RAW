from __future__ import print_function

import os
import sys
import copy
import glob

import numpy as np

RAW_DIR = os.path.dirname(os.path.abspath(__file__))
if RAW_DIR not in sys.path:
    sys.path.append(RAW_DIR)
import RAWSettings
import SASM
import SASFileIO
import SASImage
import SASExceptions

from RAWUtils import findATSASDirectory, ErrorPrinter
from RAWAnalysisWrapper import RAWAnalysisSimulator


class RAWSimulator():
    """ RAW operator """

    def __init__(self, raw_cfg_path, log_file=None, do_analysis=False):
        # load configuration
        print(raw_cfg_path)
        self._raw_settings = RAWSettings.RawGuiSettings()
        success = RAWSettings.loadSettings(self._raw_settings, raw_cfg_path)
        if success:
            self._raw_settings.set('CurrentCfg', raw_cfg_path)
        else:
            raise IOError('Load failed, config file might be corrupted.')

        # set ATSAS package
        find_atsas = self._raw_settings.get('autoFindATSAS')
        if find_atsas:
            atsas_dir = findATSASDirectory()
            self._raw_settings.set('ATSASDir', atsas_dir)

        # print to log file
        if log_file is None:
            self._stdout = sys.stdout
        else:
            self._stdout = open(log_file, mode='w')

        self.error_printer = ErrorPrinter(self._raw_settings, self._stdout)
        # create mask
        self._createMasks()

        if do_analysis:
            self.analysis_simulator = RAWAnalysisSimulator(
                self._raw_settings, self._stdout)

    # def __exit__():
    # close file
    #    self._stdout.close()

    def _createMasks(self, overwrite_cached=False):
        """Create mask from mask objects."""
        print('Please wait while creating masks...', file=self._stdout)
        mask_dict = self._raw_settings.get('Masks')
        img_dim = self._raw_settings.get('MaskDimension')

        cfg_path = self._raw_settings.get('CurrentCfg')
        cfg_root = os.path.split(cfg_path)[0]
        cached_masks = glob.glob(os.path.join(cfg_root, '*Mask*.npy'))
        cached_keys = [
            os.path.splitext(os.path.split(mask_path)[1])[0]
            for mask_path in cached_masks
        ]

        for each_key in mask_dict.keys():
            # each_key: 'TransparentBSMask', 'BeamStopMask', 'ReadOutNoiseMask'
            # mask_dict[key] = [mask_matrix, mask_object]
            masks = mask_dict[each_key][1]
            if masks != None:
                if each_key in cached_keys:
                    mask_img = np.load(
                        cached_masks[cached_keys.index(each_key)])
                else:
                    mask_img = SASImage.createMaskMatrix(img_dim, masks)
                mask_param = mask_dict[each_key]
                mask_param[0] = mask_img
                mask_param[1] = masks

                if overwrite_cached:
                    np.save(os.path.join(cfg_root, each_key), mask_img)
                elif each_key not in cached_keys:
                    np.save(os.path.join(cfg_root, each_key), mask_img)

    def get_raw_settings(self):
        return self._raw_settings

    def get_raw_settings_value(self, key):
        return self._raw_settings.get(key)

    def get_stdout(self):
        return self._stdout

    def set_raw_settings(self, **kwargs):
        for key, val in kwargs.items():
            self._raw_settings.set(key, val)

    def analyse(self, sasm):
        """
        sasm is SASM object instead of a list object
        """
        self.analysis_simulator.analyse(sasm)

    def alignSASMs(self, marked_sasm, selected_sasms, qrange):
        """Align selected sasms in-place with marked sasm.

        Parameters
        ----------
        marked_sasm : SASM object
            Reference sasm for alignment.
        selected_sasms : SASM object
            Scale selected sasm to align with reference.
        qrange : (qmin, qmax), float
            Q range for alignment and scaling. 

        Returns
        -------
        None
        """
        print('Please wait while aligning and plotting...', file=self._stdout)

        if marked_sasm is None:
            self.error_printer.showPleaseMarkItemError('subtract')
            return None
        elif len(selected_sasms) == 0:
            self.error_printer.showPleaseSelectItemsError('subtract')
            return None

        qmin, qmax = qrange

        ref_sasm = copy.deepcopy(marked_sasm)
        ref_q = ref_sasm.q
        ref_i = ref_sasm.i
        if qmax < 0:
            ref_indices = ref_q >= qmin
        else:
            ref_indices = np.logical_and(ref_q >= qmin, ref_q < qmax)

        # TODO: Improve this is rough scaling method that assumes curves are parallel.
        stat_func = np.mean
        for each_sasm in selected_sasms:
            curve_q = each_sasm.q
            curve_i = each_sasm.i

            if qmax < 0:
                curve_indices = curve_q >= qmin
            else:
                curve_indices = np.logical_and(curve_q >= qmin, curve_q < qmax)

            scaling_factor = stat_func(ref_i[ref_indices]) / stat_func(
                curve_i[curve_indices])

            each_sasm.scale(scaling_factor)

    def calibrateSASM(self, sasm):
        pass
        # sd_distance = self._raw_settings.get('SampleDistance')
        # pixel_size = self._raw_settings.get('DetectorPixelSize')
        # wavelength = self._raw_settings.get('WaveLength')
        # sasm.calibrateQ(sd_distance, pixel_size, wavelength)

    def loadSASMs(self, filename_list):
        """Load image or dat files."""
        print('Please wait while loading files...', file=self._stdout)
        sasm_list = []
        do_auto_save = self._raw_settings.get('AutoSaveOnImageFiles')

        try:
            for each_filename in filename_list:
                # file_ext = os.path.splitext(each_filename)[1]
                sasm, img = SASFileIO.loadFile(each_filename,
                                               self._raw_settings)

                if img is not None:
                    # qrange = sasm.getQrange()
                    start_point = self._raw_settings.get('StartPoint')
                    end_point = self._raw_settings.get('EndPoint')

                    if isinstance(sasm, list):
                        qrange = (start_point,
                                  len(sasm[0].getBinnedQ()) - end_point)
                        for each_sasm in sasm:
                            each_sasm.setQrange(qrange)
                    else:
                        qrange = (start_point,
                                  len(sasm.getBinnedQ()) - end_point)
                        sasm.setQrange(qrange)

                    if do_auto_save:
                        save_path = self._raw_settings.get('ProcessedFilePath')
                        try:
                            if isinstance(sasm, list):
                                for each in sasm:
                                    self.saveSASM(each, '.dat', save_path)
                            else:
                                self.saveSASM(sasm, '.dat', save_path)
                        except IOError as error:
                            self._raw_settings.set('AutoSaveOnImageFiles',
                                                   False)
                            do_auto_save = False

                if isinstance(sasm, list):
                    sasm_list.extend(sasm)
                else:
                    sasm_list.append(sasm)

        except (SASExceptions.UnrecognizedDataFormat,
                SASExceptions.WrongImageFormat) as error:
            self.error_printer.showDataFormatError(
                os.path.split(each_filename)[1])
            raise error
        except SASExceptions.HeaderLoadError as error:
            print(
                str(error),
                'Error Loading Headerfile:'
                'Please check that the header file is in the directory with the data.',
                file=self._stdout)
            raise error
        except SASExceptions.MaskSizeError as error:
            print(
                str(error),
                'Saved mask does not fit loaded image',
                file=self._stdout)
            raise error
        except SASExceptions.HeaderMaskLoadError as error:
            print(
                str(error),
                'Mask information was not found in header',
                file=self._stdout)
            raise error

        return sasm_list

    def loadIFTMs(self, filename_list):
        """Load GNOM .ift/.out files."""
        print('Please wait while loading files...', file=self._stdout)
        iftm_list = []

        try:
            for each_filename in filename_list:
                file_ext = os.path.splitext(each_filename)[1]

                if file_ext == '.ift' or file_ext == '.out':
                    iftm, _ = SASFileIO.loadFile(each_filename,
                                                 self._raw_settings)

                    if file_ext == '.ift':
                        item_colour = 'blue'
                    else:
                        item_colour = 'black'

                    if isinstance(iftm, list):
                        iftm_list.extend(iftm)
                    else:
                        iftm_list.append(iftm)

        except (SASExceptions.UnrecognizedDataFormat,
                SASExceptions.WrongImageFormat) as error:
            self.error_printer.showDataFormatError(
                os.path.split(each_filename)[1])
            raise error
        except SASExceptions.HeaderLoadError as error:
            print(
                str(error),
                'Error Loading Headerfile:'
                'Please check that the header file is in the directory with the data.',
                file=self._stdout)
            raise error
        except SASExceptions.MaskSizeError as error:
            print(
                str(error),
                'Saved mask does not fit loaded image',
                file=self._stdout)
            raise error
        except SASExceptions.HeaderMaskLoadError as error:
            print(
                str(error),
                'Mask information was not found in header',
                file=self._stdout)
            raise error

        return iftm_list

    def superimposeSASMs(self, marked_sasm, selected_sasms):
        """Superimpose seleceted sasms with marked sasm."""
        superimposed_sasms = copy.deepcopy(selected_sasms)
        if marked_sasm in selected_sasms:
            superimposed_sasms.remove(marked_sasm)

        if marked_sasm is None:
            self.error_printer.showPleaseMarkItemError('superimpose')
            return None

        if len(superimposed_sasms) == 0:
            self.error_printer.showPleaseSelectItemsError('superimpose')
            return None

        SASM.superimpose(marked_sasm, superimposed_sasms)
        return superimposed_sasms

    def subtractSASMs(self, marked_sasm, selected_sasms, yes_to_all=False):
        """Subtracts the marked sasm from other selected sasms in
        the manipulation list.

        Parameters
        ----------
        marked_sasm : SASM object
            Minuend in subtraction.
        selected_sasms: list of SASM object
            Subtrahend in subtraction.
        yes_to_all : bool, optional
            Whether to force subtraction for all sasms, if q vectors are not equal.
            Default is False.

        Returns
        -------
        subtracted_list : list of SASM object
        """
        print(
            'Please wait while subtracting and plotting...', file=self._stdout)
        do_auto_save = self._raw_settings.get('AutoSaveOnSub')

        if marked_sasm in selected_sasms:
            selected_sasms.remove(marked_sasm)

        if marked_sasm is None:
            self.error_printer.showPleaseMarkItemError('subtract')
            return None
        elif len(selected_sasms) == 0:
            self.error_printer.showPleaseSelectItemsError('subtract')
            return None

        sub_sasm = marked_sasm
        sub_qmin, sub_qmax = sub_sasm.getQrange()
        subtracted_list = []

        for sasm in selected_sasms:
            print(sasm.getParameter('filename'))
            qmin, qmax = sasm.getQrange()

            if np.all(
                    np.round(sasm.q[qmin:qmax], 5) == np.round(
                        sub_sasm.q[sub_qmin:sub_qmax], 5)) is False:
                self.error_printer.showQvectorsNotEqualWarning(sasm, sub_sasm)
                if not yes_to_all:
                    raise SASExceptions.DataNotCompatible

            try:
                subtracted_sasm = SASM.subtract(
                    sasm, sub_sasm, forced=yes_to_all)
                self.insertSasmFilenamePrefix(subtracted_sasm, 'S_')
                subtracted_list.append(subtracted_sasm)

                if do_auto_save:
                    save_path = self._raw_settings.get('SubtractedFilePath')
                    try:
                        self.saveSASM(subtracted_sasm, '.dat', save_path)
                    except IOError as error:
                        self._raw_settings.set('AutoSaveOnSub', False)
                        do_auto_save = False
                        print(
                            str(error) +
                            '\n\nAutosave of subtracted images has been disabled. If you are using a config file from a different computer please go into Advanced Options/Autosave to change the save folders, or save you config file to avoid this message next time.',
                            'Autosave Error')

            except SASExceptions.DataNotCompatible:
                self.error_printer.showSubtractionError(sasm, sub_sasm)
                return None

        return subtracted_list

    def averageSASMs(self, selected_sasms, yes_to_all=False):
        """ average selected sasms in the manipulation list.
        
        Parameters
        ----------
        selected_sasms: list of SASM object
            Sasm list for averaging
        yes_to_all : bool, optional
            Whether to force subtraction for all sasms, if q vectors are not equal.
            Default is False.

        Returns
        -------
        average_sasm : SASM object
            Return average sasm.
        """
        print('Please wait while averaging and plotting...', file=self._stdout)
        do_auto_save = self._raw_settings.get('AutoSaveOnAvgFiles')

        if len(selected_sasms) < 2:
            self.error_printer.showAverageError(2)
            return None

        try:
            avg_sasm = SASM.average(selected_sasms, forced=yes_to_all)
        except SASExceptions.DataNotCompatible:
            self.error_printer.showAverageError(3)
            return None

        self.insertSasmFilenamePrefix(avg_sasm, 'A_')

        if do_auto_save:
            save_path = self._raw_settings.get('AveragedFilePath')
            try:
                self.saveSASM(avg_sasm, '.dat', save_path)
            except IOError as error:
                self._raw_settings.set('AutoSaveOnAvgFiles', False)
                print(
                    str(error) +
                    '\n\nAutosave of averaged images has been disabled. If you are using a config file from a different computer please go into Advanced Options/Autosave to change the save folders, or save you config file to avoid this message next time.',
                    'Autosave Error')

        return avg_sasm

    def weightedAverageSASMs(self, selected_sasms, yes_to_all=False):
        """Average selected sasms with weights in the manipulation list.
       
        Parameters
        ----------
        selected_sasms: list of SASM object
            Sasm list for averaging
        yes_to_all : bool, optional
            Whether to force subtraction for all sasms, if q vectors are not equal.
            Default is False.

        Returns
        -------
        average_sasm : SASM object
            Return average sasm.
        """
        print('Please wait while averaging and plotting...', file=self._stdout)
        do_auto_save = self._raw_settings.get('AutoSaveOnAvgFiles')

        if len(selected_sasms) < 2:
            self.error_printer.showAverageError(2)
            return None

        weightByError = self._raw_settings.get('weightByError')
        weightCounter = self._raw_settings.get('weightCounter')

        if not weightByError and weightCounter == '':
            print(
                'Weighted Average Error:\n',
                '  An appropriate counter to weight the data is not selected and error weighting is not enabled. Weighted average aborted.',
                file=self._stdout)
            return

        if not weightByError:
            has_header = []
            for each_sasm in selected_sasms:
                header_keys = []
                if each_sasm.getAllParameters().has_key('counters'):
                    file_hdr = each_sasm.getParameter('counters')
                    header_keys = header_keys + file_hdr.keys()
                if each_sasm.getAllParameters().has_key('imageHeader'):
                    img_hdr = each_sasm.getParameter('imageHeader')
                    header_keys = header_keys + img_hdr.keys()

                if weightCounter in header_keys:
                    has_header.append(True)
                else:
                    has_header.append(False)
            if not np.all(has_header):
                print(
                    'Weighted Average Error:',
                    'Not all selected items had the counter value selected as the weight. Weighted average aborted.',
                    file=self._stdout)
                return None

        try:
            avg_sasm = SASM.weightedAverage(
                selected_sasms,
                weightByError,
                weightCounter,
                forced=yes_to_all)
        except SASExceptions.DataNotCompatible:
            self.error_printer.showAverageError(3)
            return

        self.insertSasmFilenamePrefix(avg_sasm, 'A_')

        if do_auto_save:
            save_path = self._raw_settings.get('AveragedFilePath')
            try:
                self.saveSASM(avg_sasm, '.dat', save_path)
            except IOError as error:
                self._raw_settings.set('AutoSaveOnAvgFiles', False)
                print(
                    str(error) +
                    '\n\nAutosave of averaged images has been disabled. If you are using a config file from a different computer please go into Advanced Options/Autosave to change the save folders, or save you config file to avoid this message next time.',
                    'Autosave Error')

        return avg_sasm

    def rebinSASMs(self, selected_sasms, rebin_factor, log_rebin):
        rebinned_list = []
        for sasm in selected_sasms:
            points = np.floor(len(sasm.q) / rebin_factor)
            if log_rebin:
                rebin_sasm = SASM.logBinning(sasm, points)
            else:
                rebin_sasm = SASM.rebin(sasm, rebin_factor)
            self.insertSasmFilenamePrefix(rebin_sasm, 'R_')
            rebinned_list.append(rebin_sasm)
        return rebinned_list

    def mergeSASMs(self, marked_sasm, selected_sasms):
        """Merge selected sasms with marked sasm."""
        if marked_sasm in selected_sasms:
            selected_sasms.remove(marked_sasm)

        if marked_sasm is None:
            self.error_printer.showPleaseMarkItemError('merge')
            return None

        merged_sasm = SASM.merge(marked_sasm, selected_sasms)
        filename = marked_sasm.getParameter('filename')
        merged_sasm.setParameter('filename', filename)
        self.insertSasmFilenamePrefix(merged_sasm, 'M_')

        return merged_sasm

    def interpolateItems(self, marked_sasm, selected_sasms):
        """Interpolate sasms with marked sasm."""
        if marked_sasm in selected_sasms:
            selected_sasms.remove(marked_sasm)

        if marked_sasm is None:
            self.error_printer.showPleaseMarkItemError('interpolate')
            return None

        interpolated_list = []
        for sasm in selected_sasms:
            interpolate_sasm = SASM.interpolateToFit(marked_sasm, sasm)
            filename = sasm.getParameter('filename')
            interpolate_sasm.setParameter('filename', filename)
            self.insertSasmFilenamePrefix(interpolate_sasm, 'I_')
            interpolated_list.append(interpolate_sasm)

        return interpolated_list

    def insertSasmFilenamePrefix(self, sasm, prefix='', extension=''):
        filename = sasm.getParameter('filename')
        new_filename, _ = os.path.splitext(filename)
        sasm.setParameter('filename', prefix + new_filename + extension)

    def saveSASM(self, sasm, filetype='dat', save_path=''):
        """Save SASM object to file."""
        newext = filetype

        filename = sasm.getParameter('filename')
        check_filename, _ = os.path.splitext(filename)
        check_filename = check_filename + newext

        filepath = os.path.join(save_path, check_filename)
        # file_exists = os.path.isfile(filepath)
        filepath = save_path

        try:
            SASFileIO.saveMeasurement(
                sasm, filepath, self._raw_settings, filetype=newext)
        except SASExceptions.HeaderSaveError:
            self.error_printer.showSaveError('header')
