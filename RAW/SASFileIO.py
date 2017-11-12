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
