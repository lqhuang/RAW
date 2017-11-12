'''
Created on Jul 11, 2010

@author: specuser

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

import os
import re
import json
import copy

import numpy as np


def createSASMFromImage(img_array, parameters = {}, x_c = None, y_c = None, mask = None,
                        readout_noise_mask = None, tbs_mask = None, dezingering = 0, dezing_sensitivity = 4):
    '''
        Load measurement. Loads an image file, does pre-processing:
        masking, radial average and returns a measurement object
    '''
    if mask != None:
        if mask.shape != img_array.shape:
            raise SASExceptions.MaskSizeError('Beamstop mask is the wrong size. Please' +
                            ' create a new mask or remove the old to make this plot.')

    if readout_noise_mask != None:
        if readout_noise_mask.shape != img_array.shape:
            raise SASExceptions.MaskSizeError('Readout-noise mask is the wrong size. Please' +
                            ' create a new mask or remove the old to make this plot.')

    if tbs_mask != None:
        if tbs_mask.shape != img_array.shape:
            raise SASExceptions.MaskSizeError('ROI Counter mask is the wrong size. Please' +
                            ' create a new mask or remove the old to make this plot.')

    try:
        [i_raw, q_raw, err_raw, qmatrix] = SASImage.radialAverage(img_array, x_c, y_c, mask, readout_noise_mask, dezingering, dezing_sensitivity)
    except IndexError, msg:
        print 'Center coordinates too large: ' + str(msg)

        x_c = img_array.shape[1]/2
        y_c = img_array.shape[0]/2

        [i_raw, q_raw, err_raw, qmatrix] = SASImage.radialAverage(img_array, x_c, y_c, mask, readout_noise_mask, dezingering, dezing_sensitivity)

        #wx.CallAfter(wx.MessageBox, "The center coordinates are too large for this image, used image center instead.",
        # "Center coordinates does not fit image", wx.OK | wx.ICON_ERROR)

    err_raw_non_nan = np.nan_to_num(err_raw)

    if tbs_mask != None:
        roi_counter = img_array[tbs_mask==1].sum()
        parameters['counters']['roi_counter'] = roi_counter

    sasm = SASM.SASM(i_raw, q_raw, err_raw_non_nan, parameters)

    return sasm


def loadMask(filename):
    ''' Loads a mask  '''

    if os.path.splitext(filename)[1] == 'msk':

        with open(filename, 'r') as FileObj:
            maskPlotParameters = cPickle.load(FileObj)

        i=0
        for each in maskPlotParameters['storedMasks']:
            each.maskID = i
            i = i + 1

        return SASImage.createMaskMatrix(maskPlotParameters), maskPlotParameters


def loadFile(filename, raw_settings, no_processing = False):
    ''' Loads a file an returns a SAS Measurement Object (SASM) and the full image if the
        selected file was an Image file

         NB: This is the function used to load any type of file in RAW
    '''
    try:
        file_type = checkFileType(filename)
        print(file_type)
    except IOError:
        raise
    except Exception:
        print(str(Exception.__str__), file=sys.stderr)
        file_type = None

    if file_type == 'image':
        try:
            sasm, img = loadImageFile(filename, raw_settings)
        except (ValueError, AttributeError):
            print('SASFileIO.loadFile : ' + str(Exception.__str__))
            raise SASExceptions.UnrecognizedDataFormat('No data could be retrieved from the file, unknown format.')

        if not RAWGlobals.usepyFAI_integration:
            try:
                sasm = SASImage.calibrateAndNormalize(sasm, img, raw_settings)
            except (ValueError, NameError):
                print(Exception.__str__)

        #Always do some post processing for image files
        if type(sasm) == list:
            for current_sasm in sasm:

                current_sasm.setParameter('config_file', raw_settings.get('CurrentCfg'))

                SASM.postProcessSasm(current_sasm, raw_settings)

                if not no_processing:
                    SASM.postProcessImageSasm(current_sasm, raw_settings)
        else:
            sasm.setParameter('config_file', raw_settings.get('CurrentCfg'))

            SASM.postProcessSasm(sasm, raw_settings)

            if not no_processing:
                SASM.postProcessImageSasm(sasm, raw_settings)

    else:
        sasm = loadAsciiFile(filename, file_type)
        img = None

        #If you don't want to post process asci files, return them as a list
        if type(sasm) != list:
            SASM.postProcessSasm(sasm, raw_settings)

    if type(sasm) != list and (sasm is None or len(sasm.i) == 0):
        raise SASExceptions.UnrecognizedDataFormat('No data could be retrieved from the file, unknown format.')

    return sasm, img


def loadImageFile(filename, raw_settings):

    img_fmt = raw_settings.get('ImageFormat')
    hdr_fmt = raw_settings.get('ImageHdrFormat')

    loaded_data, loaded_hdr = loadImage(filename, img_fmt)

    sasm_list = [None for i in range(len(loaded_data))]

    #Pre-load the flatfield file, so it's not loaded every time
    if raw_settings.get('NormFlatfieldEnabled'):
        flatfield_filename = raw_settings.get('NormFlatfieldFile')
        if flatfield_filename != None:
            flatfield_img, flatfield_img_hdr = loadImage(flatfield_filename, img_fmt)
            flatfield_hdr = loadHeader(flatfield_filename, flatfield_filename, hdr_fmt)
            flatfield_img = np.average(flatfield_img, axis=0)

    #Process all loaded images into sasms
    for i in range(len(loaded_data)):
        img = loaded_data[i]
        img_hdr = loaded_hdr[i]

        if len(loaded_data) > 1:
            temp_filename = os.path.split(filename)[1].split('.')
            if len(temp_filename) > 1:
                temp_filename[-2] = temp_filename[-2] + '_%05i' %(i)
            else:
                temp_filename[0] = temp_filename[0] + '_%05i' %(i)

            new_filename = '.'.join(temp_filename)
        else:
            new_filename = os.path.split(filename)[1]

        hdrfile_info = loadHeader(filename, new_filename, hdr_fmt)

        parameters = {'imageHeader' : img_hdr,
                      'counters'    : hdrfile_info,
                      'filename'    : new_filename,
                      'load_path'   : filename}

        for key in parameters['counters']:
            if key.lower().find('concentration') > -1 or key.lower().find('mg/ml') > -1:
                parameters['Conc'] = parameters['counters'][key]
                break

        x_c = raw_settings.get('Xcenter')
        y_c = raw_settings.get('Ycenter')

        ## Read center coordinates from header?
        if raw_settings.get('UseHeaderForCalib'):
            try:
                x_y = SASImage.getBindListDataFromHeader(raw_settings, img_hdr, hdrfile_info, keys = ['Beam X Center', 'Beam Y Center'])

                if x_y[0] != None: x_c = x_y[0]
                if x_y[1] != None: y_c = x_y[1]
            except ValueError:
                pass
            except TypeError:
                raise SASExceptions.HeaderLoadError('Error loading header, file corrupt?')

        # ********************
        # If the file is a SAXSLAB file, then get mask parameters from the header and modify the mask
        # then apply it...
        #
        # Mask should be not be changed, but should be created here. If no mask information is found, then
        # use the user created mask. There should be a force user mask setting.
        #
        # ********************

        masks = raw_settings.get('Masks')

        use_hdr_mask = raw_settings.get('UseHeaderForMask')

        if use_hdr_mask and img_fmt == 'SAXSLab300':
            try:
                mask_patches = SASImage.createMaskFromHdr(img, img_hdr, flipped = raw_settings.get('DetectorFlipped90'))
                bs_mask_patches = masks['BeamStopMask'][1]

                if bs_mask_patches != None:
                    all_mask_patches = mask_patches + bs_mask_patches
                else:
                    all_mask_patches = mask_patches

                bs_mask = SASImage.createMaskMatrix(img.shape, all_mask_patches)
            except KeyError:
                raise SASExceptions.HeaderMaskLoadError('bsmask_configuration not found in header.')

            dc_mask = masks['ReadOutNoiseMask'][0]
        else:
            bs_mask = masks['BeamStopMask'][0]
            dc_mask = masks['ReadOutNoiseMask'][0]


        tbs_mask = masks['TransparentBSMask'][0]

        # ********* WARNING WARNING WARNING ****************#
        # Hmm.. axes start from the lower left, but array coords starts
        # from upper left:
        #####################################################
        y_c = img.shape[0]-y_c

        if not RAWGlobals.usepyFAI_integration:
            # print('Using standard RAW integration')
            ## Flatfield correction.. this part gets moved to a image correction function later
            if raw_settings.get('NormFlatfieldEnabled'):
                if flatfield_filename != None:
                    img, img_hdr = SASImage.doFlatfieldCorrection(img, img_hdr, flatfield_img, flatfield_hdr)
                else:
                    pass #Raise some error

            dezingering = raw_settings.get('ZingerRemovalRadAvg')
            dezing_sensitivity = raw_settings.get('ZingerRemovalRadAvgStd')

            sasm = createSASMFromImage(img, parameters, x_c, y_c, bs_mask, dc_mask, tbs_mask, dezingering, dezing_sensitivity)

        else:
            sasm = SASImage.pyFAIIntegrateCalibrateNormalize(img, parameters, x_c, y_c, raw_settings, bs_mask, tbs_mask)

        sasm_list[i] = sasm


    return sasm_list, loaded_data


def loadOutFile(filename):

    five_col_fit = re.compile('\s*\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s+\d*[.]\d*[+eE-]*\d+\s+\d*[.]\d*[+eE-]*\d+\s+\d*[.]\d*[+eE-]*\d+\s*$')
    three_col_fit = re.compile('\s*\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s+\d*[.]\d*[+eE-]*\d+\s*$')
    two_col_fit = re.compile('\s*\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s*$')

    results_fit = re.compile('\s*Current\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s*\d*[.]?\d*[+eE-]*\d*\s*$')

    te_fit = re.compile('\s*Total\s+[Ee]stimate\s*:\s+\d*[.]\d+\s*\(?[A-Za-z\s]+\)?\s*$')
    te_num_fit = re.compile('\d*[.]\d+')
    te_quality_fit = re.compile('[Aa][A-Za-z\s]+\)?\s*$')

    p_rg_fit = re.compile('\s*Real\s+space\s*\:?\s*Rg\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*\+-\s*\d*[.]\d+[+eE-]*\d*')
    q_rg_fit = re.compile('\s*Reciprocal\s+space\s*\:?\s*Rg\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*')

    p_i0_fit = re.compile('\s*Real\s+space\s*\:?[A-Za-z0-9\s\.,+-=]*\(0\)\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*\+-\s*\d*[.]\d+[+eE-]*\d*')
    q_i0_fit = re.compile('\s*Reciprocal\s+space\s*\:?[A-Za-z0-9\s\.,+-=]*\(0\)\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*')

    qfull = []
    qshort = []
    Jexp = []
    Jerr  = []
    Jreg = []
    Ireg = []

    R = []
    P = []
    Perr = []

    outfile = []

    with open(filename) as f:
        for line in f:
            twocol_match = two_col_fit.match(line)
            threecol_match = three_col_fit.match(line)
            fivecol_match = five_col_fit.match(line)
            results_match = results_fit.match(line)
            te_match = te_fit.match(line)
            p_rg_match = p_rg_fit.match(line)
            q_rg_match = q_rg_fit.match(line)
            p_i0_match = p_i0_fit.match(line)
            q_i0_match = q_i0_fit.match(line)

            outfile.append(line)

            if twocol_match:
                # print line
                found = twocol_match.group().split()

                qfull.append(float(found[0]))
                Ireg.append(float(found[1]))

            elif threecol_match:
                #print line
                found = threecol_match.group().split()

                R.append(float(found[0]))
                P.append(float(found[1]))
                Perr.append(float(found[2]))

            elif fivecol_match:
                #print line
                found = fivecol_match.group().split()

                qfull.append(float(found[0]))
                qshort.append(float(found[0]))
                Jexp.append(float(found[1]))
                Jerr.append(float(found[2]))
                Jreg.append(float(found[3]))
                Ireg.append(float(found[4]))

            elif results_match:
                found = results_match.group().split()
                Actual_DISCRP = float(found[1])
                Actual_OSCILL = float(found[2])
                Actual_STABIL = float(found[3])
                Actual_SYSDEV = float(found[4])
                Actual_POSITV = float(found[5])
                Actual_VALCEN = float(found[6])

                if len(found) == 8:
                    Actual_SMOOTH = float(found[7])
                else:
                    Actual_SMOOTH = -1

            elif te_match:
                te_num_search = te_num_fit.search(line)
                te_quality_search = te_quality_fit.search(line)

                TE_out = float(te_num_search.group().strip())
                quality = te_quality_search.group().strip().rstrip(')').strip()


            if p_rg_match:
                found = p_rg_match.group().split()
                rg = float(found[-3])
                rger = float(found[-1])

            elif q_rg_match:
                found = q_rg_match.group().split()
                q_rg = float(found[-1])

            if p_i0_match:
                found = p_i0_match.group().split()
                i0 = float(found[-3])
                i0er = float(found[-1])

            elif q_i0_match:
                found = q_i0_match.group().split()
                q_i0 = float(found[-1])


    # Output variables not in the results file:
    #             'r'         : R,            #R, note R[-1] == Dmax
    #             'p'         : P,            #P(r)
    #             'perr'      : Perr,         #P(r) error
    #             'qlong'     : qfull,        #q down to q=0
    #             'qexp'      : qshort,       #experimental q range
    #             'jexp'      : Jexp,         #Experimental intensities
    #             'jerr'      : Jerr,         #Experimental errors
    #             'jreg'      : Jreg,         #Experimental intensities from P(r)
    #             'ireg'      : Ireg,         #Experimental intensities extrapolated to q=0

    name = os.path.basename(filename)

    results = { 'dmax'      : R[-1],        #Dmax
                'TE'        : TE_out,       #Total estimate
                'rg'        : rg,           #Real space Rg
                'rger'      : rger,         #Real space rg error
                'i0'        : i0,           #Real space I0
                'i0er'      : i0er,         #Real space I0 error
                'q_rg'      : q_rg,         #Reciprocal space Rg
                'q_i0'      : q_i0,         #Reciprocal space I0
                'out'       : outfile,      #Full contents of the outfile, for writing later
                'quality'   : quality,      #Quality of GNOM out file
                'chisq'     : Actual_DISCRP,#DISCRIP, chi squared
                'oscil'     : Actual_OSCILL,#Oscillation of solution
                'stabil'    : Actual_STABIL,#Stability of solution
                'sysdev'    : Actual_SYSDEV,#Systematic deviation of solution
                'positv'    : Actual_POSITV,#Relative norm of the positive part of P(r)
                'valcen'    : Actual_VALCEN,#Validity of the chosen interval in real space
                'smooth'    : Actual_SMOOTH,#Smoothness of the chosen interval? -1 indicates no real value, for versions of GNOM < 5.0 (ATSAS <2.8)
                'filename'  : name,         #GNOM filename
                'algorithm' : 'GNOM'        #Lets us know what algorithm was used to find the IFT
                    }

    iftm = SASM.IFTM(P, R, Perr, Jexp, qshort, Jerr, Jreg, results, Ireg, qfull)

    return [iftm]


def loadPrimusDatFile(filename):
    ''' Loads a Primus .dat format file '''

    iq_pattern = re.compile('\s*\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s+\d*[.]\d*[+eE-]*\d+\s*')

    i = []
    q = []
    err = []

    with open(filename) as f:
        lines = f.readlines()

    firstLine = lines[0]

    iq_match = iq_pattern.match(firstLine)

    if iq_match:
        firstLine = ''

    fileHeader = {'comment':firstLine}
    parameters = {'filename' : os.path.split(filename)[1],
                  'counters' : fileHeader}

    if len(lines) == 0:
        raise SASExceptions.UnrecognizedDataFormat('No data could be retrieved from the file.')

    if lines[1].find('model_intensity') > -1:
        #FoXS file with a fit! has four data columns
        is_foxs_fit=True
        comment = firstLine+'\n'+lines[0]+lines[1]
        parameters['comment']=comment
        lines = lines[2:]
        imodel = []

    else:
        is_foxs_fit = False


    for line in lines:

        if not is_foxs_fit:

            iq_match = iq_pattern.match(line)

            if iq_match:
                #print line
                found = iq_match.group().split()
                q.append(float(found[0]))
                i.append(float(found[1]))
                err.append(float(found[2]))
        else:
            found = line.split()
            q.append(float(found[0]))
            i.append(float(found[1]))
            imodel.append(float(found[2]))
            err.append(float(found[3]))


    #Check to see if there is any header from RAW, and if so get that.
    with open(filename) as f:   #Why does this need to open the file twice? Why did I do this?
        all_lines = f.readlines()

    header = []
    for j in range(len(all_lines)):
        if '### HEADER:' in all_lines[j]:
            header = all_lines[j+1:]

    hdict = None

    if len(header)>0:
        hdr_str = ''
        for each_line in header:
            hdr_str=hdr_str+each_line
        try:
            hdict = dict(json.loads(hdr_str))
            print('Loading RAW info/analysis...')
        except Exception:
            # print 'Unable to load header/analysis information. Maybe the file was not generated by RAW or was generated by an old version of RAW?'
            hdict = {}


    i = np.array(i)
    q = np.array(q)
    err = np.array(err)

    if hdict:
        for each in hdict.iterkeys():
            if each != 'filename':
                parameters[each] = hdict[each]

    sasm = SASM.SASM(i, q, err, parameters)

    if is_foxs_fit:
        parameters2 = copy.copy(parameters)
        parameters2['filename'] = os.path.splitext(os.path.split(filename)[1])[0]+'_FIT'

        sasm_model = SASM.SASM(imodel,q,err,parameters2)

        return [sasm, sasm_model]

    return sasm
