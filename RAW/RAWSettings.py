'''
Created on Jul 16, 2010

@author: Soren Nielsen

#******************************************************************************
# This file is part of RAW.
#
#    RAW is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RAW is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RAW.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************************************************
'''
from __future__ import print_function, division

# import wx, cPickle, copy, os
try:
    import cPickle as pickle  # python 2
except ModuleNotFoundError:
    import pickle  # python 3
import copy, os
import RAWGlobals, SASFileIO


CURR_ID = []
def NewId():
    CURR_ID.append(len(CURR_ID))
    return CURR_ID[-1]


class RawGuiSettings:
    '''
    This object contains all the settings nessecary for the GUI.

    '''
    def __init__(self, settings = None):
        '''
        Accepts a dictionary argument for the parameters. Uses default is no settings are given.
        '''

        self._params = settings

        if settings == None:
            self._params = {
                            #'NormalizeConst'    : [1.0,   NewId(), 'float'],
                            #'NormalizeConstChk' : [False, NewId(),  'bool'],
                            #'NormalizeM2'       : [False, NewId(),  'bool'],
                            #'NormalizeTime'     : [False, NewId(),  'bool'],
                            #'NormalizeM1'       : [False, NewId(),  'bool'],

							'NormFlatfieldEnabled'	: [False,   NewId(),  'bool'],

                            'NormAbsWater'      	: [False,   NewId(),  'bool'],
                            'NormAbsWaterI0'    	: [0.01632, NewId(),  'float'],
                            'NormAbsWaterTemp'  	: ['25',    NewId(),  'choice'],
                            'NormAbsWaterConst' 	: [1.0,     NewId(),  'float'],

                            'NormalizeTrans'    : [False, NewId(),  'bool'],
                            'Calibrate'         : [False, NewId(),  'bool'],  # Calibrate AgBe
                            'CalibrateMan'      : [True, NewId(),  'bool'],  # Calibrate manual (wavelength / distance)
                            'AutoBgSubtract'    : [False, NewId(),  'bool'],
                            'CountNormalize'    : [1.0,   NewId(), 'float'],

                            'AutoBIFT'          : [False, NewId(), 'bool'],
                            'AutoAvg'           : [False, NewId(), 'bool'],
                            'AutoAvgRemovePlots': [False, NewId(), 'bool'],

                            'AutoAvgRegExp'     : ['', NewId(), 'text'],
                            'AutoAvgNameRegExp' : ['', NewId(), 'text'],
                            'AutoAvgNoOfFrames' : [1,  NewId(),  'int'],
                            'AutoBgSubRegExp'   : ['', NewId(), 'text'],

                            'UseHeaderForMask': [False, NewId(), 'bool'],
                            'DetectorFlipped90':[False, NewId(), 'bool'],

                            #CORRECTIONS
                            'DoSolidAngleCorrection' : [True, NewId(), 'bool'],


                            #CENTER / BINNING
                            'Binsize'    : [1,     NewId(), 'int'],
                            'Xcenter'    : [512.0, NewId(), 'float'],
                            'Ycenter'    : [512.0, NewId(), 'float'],
                            'QrangeLow'  : [25,    NewId(), 'int'],
                            'QrangeHigh' : [9999,  NewId(), 'int'],
                            'StartPoint' : [0,     NewId(), 'int'],
                            'EndPoint'   : [0,     NewId(), 'int'],
                            'ImageDim'   : [[1024,1024]],

                            #MASKING
                            'SampleFile'              : [None, NewId(), 'text'],
                            'BackgroundSASM'          : [None, NewId(), 'text'],

                            'DataSECM'                : [None, NewId(), 'text'],

                            'NormAbsWaterFile'        : [None, NewId(), 'text'],
                            'NormAbsWaterEmptyFile'   : [None, NewId(), 'text'],
							'NormFlatfieldFile'		  : [None, NewId(), 'text'],

                            'TransparentBSMask'       : [None],
                            'TransparentBSMaskParams' : [None],
                            'BeamStopMask'            : [None],
                            'BeamStopMaskParams'      : [None],
                            'ReadOutNoiseMask'        : [None],
                            'ReadOutNoiseMaskParams'  : [None],
                                                                                #mask, mask_patches
                            'Masks'                   : [{'BeamStopMask'     : [None, None],
                                                          'ReadOutNoiseMask' : [None, None],
                                                          'TransparentBSMask': [None, None],
                                                         }],

                            'MaskDimension'          : [1024,1024],

                            #Q-CALIBRATION
                            'WaveLength'          : [1.0,  NewId(), 'float'],
                            'SampleDistance'      : [1000, NewId(), 'float'],
                            'SampleThickness'     : [0.0, NewId(), 'float'],
                            'ReferenceQ'          : [0.0, NewId(), 'float'],
                            'ReferenceDistPixel'  : [0,   NewId(), 'int'],
                            'ReferenceDistMm'     : [0.0, NewId(), 'float'],
                            'DetectorPixelSize'   : [70.5, NewId(), 'float'],
                            'SmpDetectOffsetDist' : [0.0, NewId(), 'float'],


							#SANS Parameters
							'SampleThickness'		: [0.1,  NewId(), 'float'],
							'DarkCorrEnabled'		: [False,   NewId(),  'bool'],
							'DarkCorrFilename'		: [None, NewId(), 'text'],


                            #DEFAULT BIFT PARAMETERS
                            'maxDmax'     : [400.0,  NewId(), 'float'],
                            'minDmax'     : [10.0,   NewId(), 'float'],
                            'DmaxPoints'  : [10,     NewId(), 'int'],
                            'maxAlpha'    : [1e10,   NewId(), 'float'],
                            'minAlpha'    : [150.0,  NewId(), 'float'],
                            'AlphaPoints' : [16,     NewId(), 'int'],
                            'PrPoints'    : [50,     NewId(), 'int'],

                            #DEFAULT pyGNOM PARAMETERS
                            'pygnomMaxAlpha'    : [60,   NewId(), 'float'],
                            'pygnomMinAlpha'    : [0.01, NewId(), 'float'],
                            'pygnomAlphaPoints' : [100,  NewId(), 'int'],
                            'pygnomPrPoints'    : [50,   NewId(), 'int'],
                            'pygnomFixInitZero' : [True, NewId(), 'bool'],

                            'pyOSCILLweight'    : [3.0, NewId(), 'float'],
                            'pyVALCENweight'    : [1.0, NewId(), 'float'],
                            'pyPOSITVweight'    : [1.0, NewId(), 'float'],
                            'pySYSDEVweight'    : [3.0, NewId(), 'float'],
                            'pySTABILweight'    : [3.0, NewId(), 'float'],
                            'pyDISCRPweight'    : [1.0, NewId(), 'float'],

                            #DEFAULT IFT PARAMETERS:
                            'IFTAlgoList'        : [['BIFT', 'pyGNOM']],
                            'IFTAlgoChoice'      : [['BIFT']],

                            #ARTIFACT REMOVAL:
                            'ZingerRemovalRadAvg'    : [False, NewId(), 'bool'],
                            'ZingerRemovalRadAvgStd' : [4.0,     NewId(), 'float'],

                            'ZingerRemoval'     : [False, NewId(), 'bool'],
                            'ZingerRemoveSTD'   : [4,     NewId(), 'int'],
                            'ZingerRemoveWinLen': [10,    NewId(), 'int'],
                            'ZingerRemoveIdx'   : [10,    NewId(), 'int'],

                            'ZingerRemovalAvgStd'  : [8,     NewId(), 'int'],
                            'ZingerRemovalAvg'     : [False, NewId(), 'bool'],

                            #SAVE DIRECTORIES
                            'ProcessedFilePath'    : [None,  NewId(), 'text'],
                            'AveragedFilePath'     : [None,  NewId(), 'text'],
                            'SubtractedFilePath'   : [None,  NewId(), 'text'],
                            'BiftFilePath'         : [None,  NewId(), 'text'],
                            'GnomFilePath'         : [None,  NewId(), 'text'],
                            'AutoSaveOnImageFiles' : [False, NewId(), 'bool'],
                            'AutoSaveOnAvgFiles'   : [False, NewId(), 'bool'],
                            'AutoSaveOnSub'        : [False, NewId(), 'bool'],
                            'AutoSaveOnBift'       : [False, NewId(), 'bool'],
                            'AutoSaveOnGnom'       : [False, NewId(), 'bool'],

                            #IMAGE FORMATS
                            'ImageFormatList'      : [SASFileIO.all_image_types],
                            'ImageFormat'          : ['Pilatus', NewId(), 'choice'],

                            #HEADER FORMATS
                            'ImageHdrFormatList'   : [SASFileIO.all_header_types],
                            'ImageHdrFormat'       : ['None', NewId(), 'choice'],

                            'ImageHdrList'         : [None],
                            'FileHdrList'          : [None],

                            'UseHeaderForCalib'    : [False, NewId(), 'bool'],

                            # Header bind list with [(Description : parameter key, header_key)]
                            'HeaderBindList'       : [{'Beam X Center'            : ['Xcenter',           None, ''],
                                                       'Beam Y Center'            : ['Ycenter',           None, ''],
                                                       'Sample Detector Distance' : ['SampleDistance',    None, ''],
                                                       'Wavelength'               : ['WaveLength',        None, ''],
                                                       'Detector Pixel Size'      : ['DetectorPixelSize', None, '']}],
                                                       # 'Number of Frames'         : ['NumberOfFrames',    None, '']}],

                            'NormalizationList'    : [None, NewId(), 'text'],
                            'EnableNormalization'  : [True, NewId(), 'bool'],

                            'OnlineFilterList'     : [None, NewId(), 'text'],
                            'EnableOnlineFiltering': [False, NewId(), 'bool'],
                            'OnlineModeOnStartup'  : [False, NewId(), 'bool'],
	                        'OnlineStartupDir'     : [None, NewId(), 'text'],

                            'MWStandardMW'         : [0, NewId(), 'float'],
                            'MWStandardI0'         : [0, NewId(), 'float'],
                            'MWStandardConc'       : [0, NewId(), 'float'],
                            'MWStandardFile'       : ['', NewId(), 'text'],

                            #Initialize volume of correlation molecular mass values.
                            #Values from Rambo, R. P. & Tainer, J. A. (2013). Nature. 496, 477-481.
                            'MWVcType'             : ['Protein', NewId(), 'choice'],
                            'MWVcAProtein'         : [1.0, NewId(), 'float'], #The 'A' coefficient for proteins
                            'MWVcBProtein'         : [0.1231, NewId(), 'float'], #The 'B' coefficient for proteins
                            'MWVcARna'             : [0.808, NewId(), 'float'], #The 'A' coefficient for proteins
                            'MWVcBRna'             : [0.00934, NewId(), 'float'], #The 'B' coefficient for proteins

                            #Initialize porod volume molecularm ass values.
                            'MWVpRho'              : [0.83*10**(-3), NewId(), 'float'], #The density in kDa/A^3

                            #Initialize Absolute scattering calibration values.
                            #Default values from Mylonas & Svergun, J. App. Crys. 2007.
                            'MWAbsRhoMprot'         : [3.22*10**23, NewId(), 'float'], #e-/g, # electrons per dry mass of protein
                            'MWAbsRhoSolv'          : [3.34*10**23, NewId(), 'float'], #e-/cm^-3, # electrons per volume of aqueous solvent
                            'MWAbsNuBar'            : [0.7425, NewId(), 'float'], #cm^3/g, # partial specific volume of the protein
                            'MWAbsR0'               : [2.8179*10**-13, NewId(), 'float'], #cm, scattering lenght of an electron

                            'CurrentCfg'         : [None],
                            'CompatibleFormats'  : [['.rad', '.tiff', '.tif', '.img', '.csv', '.dat', '.txt', '.sfrm', '.dm3', '.edf',
                                                     '.xml', '.cbf', '.kccd', '.msk', '.spr', '.h5', '.mccd', '.mar3450', '.npy', '.pnm',
                                                      '.No', '.imx_0', '.dkx_0', '.dkx_1', '.png', '.mpa', '.ift', '.sub', '.fit', '.fir',
                                                      '.out', '.mar1200', '.mar2400', '.mar2300', '.mar3600', '.int', '.ccdraw'], None],


                            #SEC Settings:
                            'secCalcThreshold'      : [1.02, NewId(), 'float'],

                            #GUI Settings:
                            'csvIncludeData'      : [None],
                            'ManipItemCollapsed'  : [False, NewId(), 'bool'] ,
                            'CurrentFilePath'     : [None],


                            'DatHeaderOnTop'      : [False, NewId(), 'bool'],
                            'PromptConfigLoad'    : [True, NewId(), 'bool'],

                            #ATSAS settings:
                            'autoFindATSAS'       : [True, NewId(), 'bool'],
                            'ATSASDir'            : ['', NewId(), 'bool'],

                            #GNOM settings
                            'gnomExpertFile'        : ['', NewId(), 'text'],
                            'gnomForceRminZero'     : ['Y', NewId(), 'choice'],
                            'gnomForceRmaxZero'     : ['Y', NewId(), 'choice'],
                            'gnomNPoints'           : [171, NewId(), 'int'],
                            'gnomInitialAlpha'      : [0.0, NewId(), 'float'],
                            'gnomAngularScale'      : [1, NewId(), 'int'],
                            'gnomSystem'            : [0, NewId(), 'int'],
                            'gnomFormFactor'        : ['', NewId(), 'text'],
                            'gnomRadius56'          : [-1, NewId(), 'float'],
                            'gnomRmin'              : [-1, NewId(), 'float'],
                            'gnomFWHM'              : [-1, NewId(), 'float'],
                            'gnomAH'                : [-1, NewId(), 'float'],
                            'gnomLH'                : [-1, NewId(), 'float'],
                            'gnomAW'                : [-1, NewId(), 'float'],
                            'gnomLW'                : [-1, NewId(), 'float'],
                            'gnomSpot'              : ['', NewId(), 'text'],
                            'gnomExpt'              : [0, NewId(), 'int'],

                            #DAMMIF settings
                            'dammifMode'            : ['Fast', NewId(), 'choice'],
                            'dammifSymmetry'        : ['P1', NewId(), 'choice'],
                            'dammifAnisometry'      : ['Unknown', NewId(), 'choice'],
                            'dammifUnit'            : ['Unknown', NewId(), 'choice'],
                            'dammifChained'         : [False, NewId(), 'bool'],
                            'dammifConstant'        : ['', NewId(), 'text'],
                            'dammifOmitSolvent'     : [True, NewId(), 'bool'],
                            'dammifDummyRadius'     : [-1, NewId(), 'float'],
                            'dammifSH'              : [-1, NewId(), 'int'],
                            'dammifPropToFit'       : [-1, NewId(), 'float'],
                            'dammifKnots'           : [-1, NewId(), 'int'],
                            'dammifCurveWeight'     : ['e', NewId(), 'choice'],
                            'dammifRandomSeed'      : ['', NewId(), 'text'],
                            'dammifMaxSteps'        : [-1, NewId(), 'int'],
                            'dammifMaxIters'        : [-1, NewId(), 'int'],
                            'dammifMaxStepSuccess'  : [-1, NewId(), 'int'],
                            'dammifMinStepSuccess'  : [-1, NewId(), 'int'],
                            'dammifTFactor'         : [-1, NewId(), 'float'],
                            'dammifRgPen'           : [-1, NewId(), 'float'],
                            'dammifCenPen'          : [-1, NewId(), 'float'],
                            'dammifLoosePen'        : [-1, NewId(), 'float'],
                            'dammifAnisPen'         : [-1, NewId(), 'float'],
                            'dammifMaxBeadCount'    : [-1, NewId(), 'int'],
                            'dammifReconstruct'     : [15, NewId(), 'int'],
                            'dammifDamaver'         : [True, NewId(), 'bool'],
                            'dammifDamclust'        : [False, NewId(), 'bool'],
                            'dammifRefine'          : [True, NewId(), 'bool'],
                            'dammifProgram'         : ['DAMMIF', NewId(), 'choice'],

                            #DAMMIN settings that are not included in DAMMIF settings
                            'damminInitial'         : ['S', NewId(), 'choice'], #Initial DAM
                            'damminKnots'           : [20, NewId(), 'int'],
                            'damminConstant'        : [0, NewId(), 'float'],
                            'damminDiameter'        : [-1, NewId(), 'float'],
                            'damminPacking'         : [-1, NewId(), 'float'],
                            'damminCoordination'    : [-1, NewId(), 'float'],
                            'damminDisconPen'       : [-1, NewId(), 'float'],
                            'damminPeriphPen'       : [-1, NewId(), 'float'],
                            'damminCurveWeight'     : ['1', NewId(), 'choice'],
                            'damminAnealSched'      : [-1, NewId(), 'float'],

                            #Weighted Average Settings
                            'weightCounter'         : ['', NewId(), 'choice'],
                            'weightByError'         : [True, NewId(), 'bool'],
                            }

    def get(self, key):
        return self._params[key][0]

    def set(self, key, value):
        self._params[key][0] = value

    def getId(self, key):
        return self._params[key][1]

    def getType(self, key):
        return self._params[key][2]

    def getIdAndType(self, key):
        return (self._params[key][1], self._params[key][2])

    def getAllParams(self):
        return self._params



def fixBackwardsCompatibility(raw_settings):

    #Backwards compatibility for BindList:
    bind_list = raw_settings.get('HeaderBindList')
    for each_key in bind_list.keys():
        if len(bind_list[each_key]) == 2:
            bind_list[each_key] = [bind_list[each_key][0], bind_list[each_key][1], '']


def loadSettings(raw_settings, loadpath):

    file_obj = open(loadpath, 'rb')
    try:
        try:
            loaded_param = pickle.load(file_obj, encoding='latin1')  # python 3
        except TypeError:
            loaded_param = pickle.load(file_obj)  # python 2
    except (KeyError, EOFError, ImportError, IndexError, AttributeError, pickle.UnpicklingError) as e:
        print('Error type: %s, error: %s' %(type(e).__name__, e))
        file_obj.close()
        return False
    file_obj.close()

    keys = loaded_param.keys()
    all_params = raw_settings.getAllParams()

    for each_key in keys:
        if each_key in all_params:
            all_params[each_key][0] = copy.copy(loaded_param[each_key])
        else:
            print('WARNING: ' + str(each_key) + " not found in RAWSettings.")

    # main_frame = wx.FindWindowByName('MainFrame')
    # main_frame.queueTaskInWorkerThread('recreate_all_masks', None)

    postProcess(raw_settings)

    return True

def postProcess(raw_settings):
    fixBackwardsCompatibility(raw_settings)

    dir_check_list = [('AutoSaveOnImageFiles', 'ProcessedFilePath'), ('AutoSaveOnAvgFiles', 'AveragedFilePath'),
                    ('AutoSaveOnSub', 'SubtractedFilePath'), ('AutoSaveOnBift', 'BiftFilePath'),
                    ('AutoSaveOnGnom', 'GnomFilePath'), ('OnlineModeOnStartup', 'OnlineStartupDir')
                    ]

    file_check_list = [('NormFlatfieldEnabled', 'NormFlatfieldFile')]

    change_list = []

    message_dir = {'AutoSaveOnImageFiles'   : '- AutoSave processed image files',
                    'AutoSaveOnAvgFiles'    : '- AutoSave averaged files',
                    'AutoSaveOnSub'         : '- AutoSave subtracted files',
                    'AutoSaveOnBift'        : '- AutoSave BIFT files',
                    'AutoSaveOnGnom'        : '- AutoSave GNOM files',
                    'OnlineModeOnStartup'   : '- Start online mode when RAW starts',
                    'NormFlatfieldEnabled'  : '- Apply a flatfield correction'
                    }

    for item in dir_check_list:
        if raw_settings.get(item[0]):
            if not os.path.isdir(raw_settings.get(item[1])):
                raw_settings.set(item[0], False)
                change_list.append(item[0])

    for item in file_check_list:
        if raw_settings.get(item[0]):
            if not os.path.isfile(raw_settings.get(item[1])):
                raw_settings.set(item[0], False)
                change_list.append(item[0])

    text = ''
    for item in change_list:
        text = text + message_dir[item] +'\n'

    # if len(change_list) > 0:
    #     wx.CallAfter(wx.MessageBox, 'The following settings have been disabled because the appropriate directory/file could not be found:\n'+text+'\nIf you are using a config file from a different computer please go into Advanced Options to change the settings, or save you config file to avoid this message next time.', 'Load Settings Warning', style = wx.ICON_ERROR | wx.OK | wx.STAY_ON_TOP)

    # if raw_settings.get('autoFindATSAS'):
    #     main_frame = wx.FindWindowByName('MainFrame')
    #     main_frame.findAtsas()

    return

def saveSettings(raw_settings, savepath):
    param_dict = raw_settings.getAllParams()
    keys = param_dict.keys()

    exclude_keys = ['ImageFormatList', 'ImageHdrFormatList', 'BackgroundSASM', 'CurrentCfg', 'csvIncludeData', 'CompatibleFormats', 'DataSECM']

    save_dict = {}

    for each_key in keys:
        if each_key not in exclude_keys:
            save_dict[each_key] = param_dict[each_key][0]

    save_dict = copy.deepcopy(save_dict)

    #remove big mask arrays from the cfg file
    masks = save_dict['Masks']

    for key in masks.keys():
        masks[key][0] = None

    file_obj = open(savepath, 'wb')
    try:
        # cPickle.dump(save_dict, file_obj, cPickle.HIGHEST_PROTOCOL)
        pickle.dump(save_dict, file_obj, protocol=2)
    except Exception as e:
        print('<Error> type: %s, message: %s' %(type(e).__name__, e))
        file_obj.close()
        return False

    file_obj.close()

    ## Make a backup of the config file in case of crash:
    backup_file = os.path.join(RAWGlobals.RAWWorkDir, 'backup.cfg')

    FileObj = open(backup_file, 'wb')
    try:
        # pickle.dump(save_dict, FileObj, pickle.HIGHEST_PROTOCOL)
        pickle.dump(save_dict, FileObj, protocol=2)
    except Exception as e:
        print('Error type: %s, error: %s' %(type(e).__name__, e))
        FileObj.close()
        return False
    FileObj.close()

    dummy_settings = RawGuiSettings()

    test_load = loadSettings(dummy_settings, savepath)

    if not test_load:
        os.remove(savepath)

    return test_load



# Table from http://physchem.kfunigraz.ac.at/sm/Services.htm
water_scattering_table = {0 : 0.01692,
                        1 : 0.01686,
                        2 : 0.01680,
                        3 : 0.01675,
                        4 : 0.01670,
                        5 : 0.01665,
                        6 : 0.01661,
                        7 : 0.01657,
                        8 : 0.01653,
                        9 : 0.01650,
                        10 : 0.01647,
                        11 : 0.01645,
                        12 : 0.01642,
                        13 : 0.01640,
                        14 : 0.01638,
                        15 : 0.01637,
                        16 : 0.01635,
                        17 : 0.01634,
                        18 : 0.01633,
                        19 : 0.01633,
                        20 : 0.01632,
                        21 : 0.01632,
                        22 : 0.01632,
                        23 : 0.01632,
                        24 : 0.01632,
                        25 : 0.01633,
                        26 : 0.01634,
                        27 : 0.01635,
                        28 : 0.01636,
                        29 : 0.01637,
                        30 : 0.01638,
                        31 : 0.01640,
                        32 : 0.01641,
                        33 : 0.01643,
                        34 : 0.01645,
                        35 : 0.01647,
                        36 : 0.01650,
                        37 : 0.01652,
                        38 : 0.01655,
                        39 : 0.01658,
                        40 : 0.01660,
                        41 : 0.01663,
                        42 : 0.01666,
                        43 : 0.01670,
                        44 : 0.01673,
                        45 : 0.01677,
                        46 : 0.01680,
                        47 : 0.01684,
                        48 : 0.01688,
                        49 : 0.01692,
                        50 : 0.01696,
                        51 : 0.01700,
                        52 : 0.01704,
                        53 : 0.01709,
                        54 : 0.01713,
                        55 : 0.01718,
                        56 : 0.01723,
                        57 : 0.01728,
                        58 : 0.01732,
                        59 : 0.01738,
                        60 : 0.01743,
                        61 : 0.01748,
                        62 : 0.01753,
                        63 : 0.01759,
                        64 : 0.01764,
                        65 : 0.01770,
                        66 : 0.01776,
                        67 : 0.01781,
                        68 : 0.01787,
                        69 : 0.01793,
                        70 : 0.01800,
                        71 : 0.01806,
                        72 : 0.01812,
                        73 : 0.01818,
                        74 : 0.01825,
                        75 : 0.01831,
                        76 : 0.01838,
                        77 : 0.01845,
                        78 : 0.01852,
                        79 : 0.01859,
                        80 : 0.01866,
                        81 : 0.01873,
                        82 : 0.01880,
                        83 : 0.01887,
                        84 : 0.01895,
                        85 : 0.01902,
                        86 : 0.01909,
                        87 : 0.01917,
                        88 : 0.01925,
                        89 : 0.01932,
                        90 : 0.01940,
                        91 : 0.01948,
                        92 : 0.01956,
                        93 : 0.01964,
                        94 : 0.01973,
                        95 : 0.01981,
                        96 : 0.01989,
                        97 : 0.01998,
                        98 : 0.02006,
                        99 : 0.02015,
                        100 : 0.02023}

