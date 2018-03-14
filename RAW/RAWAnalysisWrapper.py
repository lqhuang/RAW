from __future__ import print_function, division

import os
import sys
import re
import copy
import platform
import subprocess

from scipy import polyval, polyfit
import numpy as np

RAW_DIR = os.path.dirname(os.path.abspath(__file__))
if RAW_DIR not in sys.path:
    sys.path.append(RAW_DIR)
import SASM
import SASFileIO
import SASExceptions
import SASCalc
from RAWUtils import ErrorPrinter


class GuinierAnalyzer():
    """Wrapper for GuinierControlPanel"""

    def __init__(self, raw_settings, stdout):
        self.raw_settings = raw_settings
        self._stdout = stdout

        self.curr_sasm = None
        self.old_analysis = None

        self.spinctrlIDs = {
            'qstart': 0,
            'qend': 0,
        }

        self.infodata = {
            'I0': (0, 0),
            'Rg': (0, 0),
            'qRg_max': 0,
            'qRg_min': 0,
            'rsq': 0,
        }

    def analyse(self, sasm):
        """
        sasm: (not list) this will modify input sasm object
        """
        self.curr_sasm = sasm
        self._runAutoRg()
        self._calcFit()
        self._saveInfo()

    def _runAutoRg(self):
        rg, rger, i0, i0er, idx_min, idx_max = SASCalc.autoRg(self.curr_sasm)

        self.spinctrlIDs['rg'] = rg
        self.spinctrlIDs['rger'] = rger
        self.spinctrlIDs['i0'] = i0
        self.spinctrlIDs['i0er'] = i0er

        if rg == -1:
            print(
                'AutoRG Failed:',
                'AutoRG could not find a suitable interval to calculate Rg.',
                file=self._stdout)
        else:
            try:
                self.curr_sasm.q[int(idx_min)]
                self.curr_sasm.q[int(idx_max)]
                self.spinctrlIDs['qstart'] = int(idx_min)
                self.spinctrlIDs['qend'] = int(idx_max)
            except IndexError:
                print(
                    'AutoRG Failed:',
                    'AutoRG did not produce a useable result. Please report this to the developers.',
                    file=self._stdout)
                raise IndexError()

    def _calcFit(self):
        """ calculate fit and statistics """
        qstart = self.spinctrlIDs['qstart']
        qend = self.spinctrlIDs['qend']

        q_roi = self.curr_sasm.q[qstart:qend]
        i_roi = self.curr_sasm.i[qstart:qend]
        err_roi = self.curr_sasm.err[qstart:qend]

        x = np.power(q_roi, 2)
        y = np.log(i_roi)
        err = y * np.absolute(err_roi / i_roi)

        #Remove NaN and Inf values:
        x = x[np.where(np.isnan(y) == False)]
        err = err[np.where(np.isnan(y) == False)]
        y = y[np.where(np.isnan(y) == False)]

        x = x[np.where(np.isinf(y) == False)]
        err = err[np.where(np.isinf(y) == False)]
        y = y[np.where(np.isinf(y) == False)]

        #Get 1.st order fit:
        ar, br = polyfit(x, y, 1)

        #This uses error weighted points to calculate the Rg. Probably the correct way to do it, but different
        #from how it has always been done.
        # f = lambda x, a, b: a+b*x
        # opt, cov = scipy.optimize.curve_fit(f, x, y, sigma = err, absolute_sigma = True)
        # ar = opt[1]
        # br = opt[0]

        #Obtain fit values:
        y_fit = polyval([ar, br], x)

        #Get fit statistics:
        error = y - y_fit
        SS_tot = np.sum(np.power(y - np.mean(y), 2))
        SS_err = np.sum(np.power(error, 2))
        rsq = 1 - SS_err / SS_tot

        I0 = br
        Rg = np.sqrt(-3 * ar)

        if np.isnan(Rg):
            Rg = 0

        ######## CALCULATE ERROR ON PARAMETERS ###############

        N = len(error)
        stde = SS_err / (N - 2)
        std_slope = stde * np.sqrt((1 / N) + (
            np.power(np.mean(x), 2) / np.sum(np.power(x - np.mean(x), 2))))
        std_interc = stde * np.sqrt(1 / np.sum(np.power(x - np.mean(x), 2)))

        ######################################################

        if np.isnan(std_slope):
            std_slope = -1
        if np.isnan(std_interc):
            std_interc = -1

        newInfo = {
            'I0': (np.exp(I0), std_interc),
            'Rg': (Rg, std_slope),
            'qRg_max': Rg * np.sqrt(x[-1]),
            'qRg_min': Rg * np.sqrt(x[0]),
            'rsq': rsq
        }

        return x, y_fit, br, error, newInfo

    def _saveInfo(self):
        x_fit, y_fit, I0, error, newInfo = self._calcFit()

        for key, value in newInfo.items():
            self.infodata[key] = value

        info_dict = copy.deepcopy(self.infodata)

        qstart_val = self.spinctrlIDs['qstart']
        qend_val = self.spinctrlIDs['qend']

        info_dict['qStart'] = qstart_val
        info_dict['qEnd'] = qend_val

        analysis_dict = self.curr_sasm.getParameter('analysis')
        analysis_dict['guinier'] = info_dict


class GNOMAnalyzer():
    """Wrapper for GNOMControlPanel """

    def __init__(self, raw_settings, stdout=None):
        self.raw_settings = raw_settings
        if stdout is None:
            self._stdout = sys.stdout
        else:
            self._stdout = stdout

        self.gnom_settings = {
            'expert': self.raw_settings.get('gnomExpertFile'),
            'rmin_zero': self.raw_settings.get('gnomForceRminZero'),
            'rmax_zero': self.raw_settings.get('gnomForceRmaxZero'),
            'npts': self.raw_settings.get('gnomNPoints'),
            'alpha': self.raw_settings.get('gnomInitialAlpha'),
            'angular': self.raw_settings.get('gnomAngularScale'),
            'system': self.raw_settings.get('gnomSystem'),
            'form': self.raw_settings.get('gnomFormFactor'),
            'radius56': self.raw_settings.get('gnomRadius56'),
            'rmin': self.raw_settings.get('gnomRmin'),
            'fwhm': self.raw_settings.get('gnomFWHM'),
            'ah': self.raw_settings.get('gnomAH'),
            'lh': self.raw_settings.get('gnomLH'),
            'aw': self.raw_settings.get('gnomAW'),
            'lw': self.raw_settings.get('gnomLW'),
            'spot': self.raw_settings.get('gnomSpot'),
            'expt': self.raw_settings.get('gnomExpt'),
        }

        # self.out_list = {}
        self.curr_sasm = None
        self.curr_iftm = None

        self.spinctrlIDs = {
            'qstart': 0,
            'qend': 0,
            'dmax': 0,
        }

        self.infodata = {
            'guinierI0': (0, 0),
            'guinierRg': (0, 0),
            'gnomI0': (0, 0),
            'gnomRg': (0, 0),
            'TE': 0,  # 'Total Estimate'
            'gnomQuality': 0,
            'chisq': 0,
        }

        self._getGnomVersion()

    def _getGnomVersion(self):
        """Checks if we have gnom4 or gnom5"""
        atsasDir = self.raw_settings.get('ATSASDir')

        opsys = platform.system()

        if opsys == 'Windows':
            dammifDir = os.path.join(atsasDir, 'dammif.exe')
        else:
            dammifDir = os.path.join(atsasDir, 'dammif')

        if os.path.exists(dammifDir):
            process = subprocess.Popen(
                '%s -v' % (dammifDir),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True
            )  #gnom4 doesn't do a proper -v!!! So use something else
            output, error = process.communicate()
            output = output.strip().decode('utf-8')
            error = error.strip().decode('utf-8')

            dammif_re = 'ATSAS\s*\d+[.]\d+[.]\d*'
            version_match = re.search(dammif_re, output)
            version = version_match.group().split()[-1]

            if int(version.split('.')[0]) > 2 or (
                    int(version.split('.')[0]) == 2
                    and int(version.split('.')[1]) >= 8):
                self.new_gnom = True
            else:
                self.new_gnom = False

    def analyse(self, sasm):
        self.curr_sasm = sasm

        analysis_dict = sasm.getParameter('analysis')
        if 'guinier' in analysis_dict:
            self.spinctrlIDs['qstart'] = analysis_dict['guinier']['qStart']
            self.spinctrlIDs['qend'] = analysis_dict['guinier']['qEnd']

        self.curr_iftm = self._initGNOM(sasm)
        self._saveInfo()

    def _initGNOM(self, sasm):
        analysis_dict = sasm.getParameter('analysis')
        if 'GNOM' in analysis_dict:
            iftm = self._initGnomValues(sasm)
            assert False
        else:
            path = self.raw_settings.get('GnomFilePath')  # TODO: temp files?
            cwd = os.getcwd()
            savename = 't_dat.dat'

            while os.path.isfile(os.path.join(path, savename)):
                savename = 't' + savename

            save_sasm = SASM.SASM(
                copy.deepcopy(sasm.i), copy.deepcopy(sasm.q),
                copy.deepcopy(sasm.err),
                copy.deepcopy(sasm.getAllParameters()))
            save_sasm.setParameter('filename', savename)
            save_sasm.setQrange(sasm.getQrange())

            try:
                SASFileIO.saveMeasurement(
                    save_sasm, path, self.raw_settings, filetype='.dat')
            except SASExceptions.HeaderSaveError as error:
                printer = ErrorPrinter(self.raw_settings, self._stdout)
                printer.showSaveError('header')

            os.chdir(path)
            try:
                init_iftm = SASCalc.runDatgnom(savename, sasm,
                                               self.raw_settings)
            except SASExceptions.NoATSASError as error:
                print(
                    'Error running GNOM/DATGNOM:',
                    str(error),
                    file=self._stdout)
                self.cleanupGNOM(path, savename=savename)
                os.chdir(cwd)
                return None
            os.chdir(cwd)

            if init_iftm is None:
                outname = 't_datgnom.out'
                while os.path.isfile(outname):
                    outname = 't' + outname

                if 'guinier' in analysis_dict:
                    rg = float(analysis_dict['guinier']['Rg'][0])  # TODO: [0]?
                    dmax = int(rg * 3.)  #Mostly arbitrary guess at Dmax
                    print("????:", rg)
                else:
                    print(
                        'No DMAX found warning:',
                        'No Guinier analysis found, arbirary value 80 will be set to DMAX.',
                        file=self._stdout)
                    dmax = 80  #Completely arbitrary default setting for Dmax

                os.chdir(path)
                try:
                    init_iftm = SASCalc.runGnom(
                        savename,
                        outname,
                        dmax,
                        self.gnom_settings,
                        new_gnom=self.new_gnom,
                        raw_settings=self.raw_settings)
                except SASExceptions.NoATSASError as error:
                    print(
                        'Error running GNOM/DATGNOM',
                        str(error),
                        file=self._stdout)
                    self.cleanupGNOM(path, savename=savename, outname=outname)
                    os.chdir(cwd)
                    return None
                os.chdir(cwd)

                self.cleanupGNOM(path, outname=outname)

            self.cleanupGNOM(path, savename=savename)

            iftm = self._initDatgnomValues(sasm, init_iftm)

        # plotPanel.plotPr(iftm)
        return iftm

    def _initGnomValues(self, sasm):
        dmax = sasm.getParameter('analysis')['GNOM']['Dmax']
        iftm = self._calcGNOM(dmax)
        return iftm

    def _initDatgnomValues(self, sasm, iftm):
        dmax = int(round(iftm.getParameter('dmax')))
        if dmax != iftm.getParameter('dmax'):
            iftm = self._calcGNOM(dmax)
        return iftm

    def _calcGNOM(self, dmax):
        start = int(self.spinctrlIDs['qstart'])
        end = int(self.spinctrlIDs['qend'])
        self.gnom_settings['npts'] = 0

        path = self.raw_settings.get('GnomFilePath')  # TODO: temp path
        cwd = os.getcwd()
        savename = 't_dat.dat'

        while os.path.isfile(os.path.join(path, savename)):
            savename = 't' + savename

        outname = 't_out.out'
        while os.path.isfile(os.path.join(path, outname)):
            outname = 't' + outname

        save_sasm = SASM.SASM(
            copy.deepcopy(self.curr_sasm.i), copy.deepcopy(self.curr_sasm.q),
            copy.deepcopy(self.curr_sasm.err),
            copy.deepcopy(self.curr_sasm.getAllParameters()))
        save_sasm.setParameter('filename', savename)
        save_sasm.setQrange((start, end))

        try:
            SASFileIO.saveMeasurement(
                save_sasm, path, self.raw_settings, filetype='.dat')
        except SASExceptions.HeaderSaveError as error:
            printer = ErrorPrinter(self.raw_settings, self._stdout)
            printer.showSaveError('header')

        os.chdir(path)
        try:
            iftm = SASCalc.runGnom(
                savename,
                outname,
                dmax,
                self.gnom_settings,
                new_gnom=self.new_gnom,
                raw_settings=self.raw_settings)
        except SASExceptions.NoATSASError as error:
            print('Error running GNOM/DATGNOM:', str(error), file=self._stdout)
            self.cleanupGNOM(path, savename, outname)
            os.chdir(cwd)
            return None

        os.chdir(cwd)
        # self.cleanupGNOM(path, savename, outname)

        return iftm

    def cleanupGNOM(self, path, savename='', outname=''):
        savefile = os.path.join(path, savename)
        outfile = os.path.join(path, outname)

        if savename != '':
            if os.path.isfile(savefile):
                try:
                    os.remove(savefile)
                except Exception as error:
                    print(
                        error,
                        'GNOM cleanup failed to remove the .dat file!',
                        file=self._stdout)

        if outname != '':
            if os.path.isfile(outfile):
                try:
                    os.remove(outfile)
                except Exception as error:
                    print(
                        error,
                        'GNOM cleanup failed to remove the .out file!',
                        file=self._stdout)

    def _saveInfo(self):
        gnom_results = {}

        dmax = int(round(self.curr_iftm.getParameter('dmax')))
        start_idx = self.spinctrlIDs['qstart']
        end_idx = self.spinctrlIDs['qend']

        gnom_results['Dmax'] = dmax
        gnom_results['Total_Estimate'] = self.curr_iftm.getParameter('TE')
        gnom_results['Real_Space_Rg'] = self.curr_iftm.getParameter('rg')
        gnom_results['Real_Space_I0'] = self.curr_iftm.getParameter('i0')
        gnom_results['qStart'] = self.curr_sasm.q[start_idx]
        gnom_results['qEnd'] = self.curr_sasm.q[end_idx]
        # gnom_results['GNOM_ChiSquared'] = self.curr_iftm['chisq']
        # gnom_results['GNOM_Quality_Assessment'] = self.curr_iftm['gnomQuality']

        analysis_dict = self.curr_sasm.getParameter('analysis')
        analysis_dict['GNOM'] = gnom_results

        iftm = self.curr_iftm
        iftm.setParameter(
            'filename',
            os.path.splitext(self.curr_sasm.getParameter('filename'))[0] +
            '.out')

        if self.raw_settings.get('AutoSaveOnGnom'):
            if os.path.isdir(self.raw_settings.get('GnomFilePath')):
                self.saveIFTM(iftm, self.raw_settings.get('GnomFilePath'))
            else:
                self.raw_settings.set('GnomFilePath', False)
                print(
                    'Autosave Error:',
                    'The folder:\n' + self.raw_settings.get('GNOMFilePath') +
                    '\ncould not be found. Autosave of GNOM files has been disabled. If you are using a config file from a different computer please go into Advanced Options/Autosave to change the save folders, or save you config file to avoid this message next time.',
                    file=self._stdout)

    def saveIFTM(self, iftm, save_path):
        """Save IFTM object to file."""
        if iftm.getParameter('algorithm') == 'GNOM':
            newext = '.out'
        else:
            newext = '.ift'

        filename = iftm.getParameter('filename')
        check_filename, ext = os.path.splitext(filename)
        check_filename = check_filename + newext

        filepath = os.path.join(save_path, check_filename)
        # file_exists = os.path.isfile(filepath)
        filepath = save_path

        try:
            SASFileIO.saveMeasurement(
                iftm, filepath, self.raw_settings, filetype=newext)
        except SASExceptions.HeaderSaveError:
            printer = ErrorPrinter(self.raw_settings, self._stdout)
            printer.showSaveError('header')


class MolecularWeightAnalyzer():
    """Wrapper for MolWeightFrame"""

    def __init__(self, raw_settings):
        self.raw_settings = raw_settings

        self.infodata = {
            'I0': ('I0 :', 0, 0),
            'Rg': ('Rg :', 0, 0),
        }

    def analyse(self, sasm):
        self.curr_sasm = sasm

        if 'molecularWeight' in self.curr_sasm.getParameter('analysis'):
            self.old_analysis = copy.deepcopy(
                self.curr_sasm.getParameter('analysis')['molecularWeight'])

    def updateGuinierInfo(self):
        pass


class BIFTAnalyzer():
    """Wrapper for BIFTControlPanel"""

    def __init__(self, raw_settings):
        self.raw_settings = raw_settings

        self.bift_settings = (self.raw_settings.get('PrPoints'),
                              self.raw_settings.get('maxAlpha'),
                              self.raw_settings.get('minAlpha'),
                              self.raw_settings.get('AlphaPoints'),
                              self.raw_settings.get('maxDmax'),
                              self.raw_settings.get('minDmax'),
                              self.raw_settings.get('DmaxPoints'))

        self.infodata = {
            'dmax': ('Dmax :', 0),
            'alpha': ('Alpha :', 0),
            'guinierI0': ('I0 :', 0),
            'guinierRg': ('Rg :', 0),
            'biftI0': ('I0 :', 0),
            'biftRg': ('Rg :', 0),
            'chisq': ('chi^2 (fit) :', 0),
        }

        self.iftm = None

    def analyse(self, sasm):
        self.curr_sasm = sasm
        if 'BIFT' in self.curr_sasm.getParameter('analysis'):
            self.old_analysis = copy.deepcopy(
                self.curr_sasm.getParameter('analysis')['BIFT'])


class RAWAnalysisSimulator():
    """RAW Data Analysis"""
    ANALYSIS = {
        'guinier': GuinierAnalyzer,
        'GNOM': GNOMAnalyzer,
        'molecularWeight': None,
        'BIFT': None,
    }

    def __init__(self, raw_settings, stdout=None):
        self.raw_settings = raw_settings

        if stdout is None:
            self._stdout = sys.stdout
        else:
            self._stdout = stdout

        self._analyzer = dict()
        for key, analyzer_cls in self.ANALYSIS.items():
            if analyzer_cls is not None:
                self._analyzer[key] = analyzer_cls(self.raw_settings,
                                                   self._stdout)

    def analyse(self, sasm):
        self._analyzer['guinier'].analyse(sasm)
        self._analyzer['GNOM'].analyse(sasm)
