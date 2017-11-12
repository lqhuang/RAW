'''
Created on Jul 7, 2010

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

import numpy as np
from scipy import optimize
import SASExceptions, SASParser, SASCalib, SASM, RAWGlobals
import wx, sys, math

try:
    import pyFAI
    RAWGlobals.usepyFAI = True
except:
    RAWGlobals.usepyFAI = False

# If C extensions have not been built, build them:
if RAWGlobals.compiled_extensions:
    try:
        import ravg_ext

    except ImportError:
        import SASbuild_Clibs
        try:
            SASbuild_Clibs.buildAll()
            import ravg_ext

        except Exception, e:
            RAWGlobals.compiled_extensions = False
            print e

import polygonMasking as polymask

class Mask:
    ''' Mask super class. Masking is used for masking out unwanted regions
    of an image '''

    def __init__(self, id, img_dim, type, negative = False):

        self._is_negative_mask = negative
        self._img_dimension = img_dim            # need image Dimentions to get the correct fill points
        self._mask_id = id
        self._type = type
        self._points = None

    def setAsNegativeMask(self):
        self._is_negative_mask = True

    def setAsPositiveMask(self):
        self._is_negative_mask = False

    def isNegativeMask(self):
        return self._is_negative_mask

    def getPoints(self):
        return self._points

    def setPoints(self, points):
        self._points = points

    def setId(self, id):
        self._mask_id = id

    def getId(self):
        return self._mask_id

    def getType(self):
        return self._type

    def getFillPoints(self):
        pass    # overridden when inherited

class CircleMask(Mask):
    ''' Create a circular mask '''

    def __init__(self, center_point, radius_point, id, img_dim, negative = False):

        Mask.__init__(self, id, img_dim, 'circle', negative)

        self._points = [center_point, radius_point]
        self._radius = abs(self._points[1][0] - self._points[0][0])

    def getRadius(self):
        return self._radius

    def grow(self, pixels):
        ''' Grow the circle by extending the radius by a number
        of pixels '''

        xy_c, xy_r = self._points

        x_c, y_c = xy_c
        x_r, y_r = xy_r

        if x_r > x_c:
            x_r = x_r + pixels
        else:
            x_r = x_r - pixels

        self.setPoints([(x_c,y_c), (x_r,y_r)])

    def shrink(self, pixels):
        ''' Shrink the circle by shortening the radius by a number
        of pixels '''

        xy_c, xy_r = self._points

        x_c, y_c = xy_c
        x_r, y_r = xy_r

        if x_r > x_c:
            x_r = x_r - pixels
        else:
            x_r = x_r + pixels

        self.setPoints([(x_c,y_c), (x_r,y_r)])

    def setPoints(self, points):
        self._points = points
        self.radius = abs(points[1][0] - points[0][0])

    def getFillPoints(self):
        ''' Really Clumsy! Can be optimized alot! triplicates the points in the middle!'''

        radiusC = abs(self._points[1][0] - self._points[0][0])

        #P = bresenhamCirclePoints(radiusC, imgDim[0] - self._points[0][1], self._points[0][0])
        P = calcBresenhamCirclePoints(radiusC, self._points[0][1], self._points[0][0])

        fillPoints = []

        for i in range(0, int(len(P)/8) ):
            Pp = P[i*8 : i*8 + 8]

            q_ud1 = ( Pp[0][0], range( int(Pp[1][1]), int(Pp[0][1]+1)) )
            q_ud2 = ( Pp[2][0], range( int(Pp[3][1]), int(Pp[2][1]+1)) )

            q_lr1 = ( Pp[4][1], range( int(Pp[6][0]), int(Pp[4][0]+1)) )
            q_lr2 = ( Pp[5][1], range( int(Pp[7][0]), int(Pp[5][0]+1)) )

            for i in range(0, len(q_ud1[1])):
                fillPoints.append( (q_ud1[0], q_ud1[1][i]) )
                fillPoints.append( (q_ud2[0], q_ud2[1][i]) )
                fillPoints.append( (q_lr1[1][i], q_lr1[0]) )
                fillPoints.append( (q_lr2[1][i], q_lr2[0]) )

        return fillPoints

class RectangleMask(Mask):
    ''' create a retangular mask '''

    def __init__(self, first_point, second_point, id, img_dim, negative = False):

        Mask.__init__(self, id, img_dim, 'rectangle', negative)
        self._points = [first_point, second_point]

    def grow(self, pixels):

        xy1, xy2 = self._points

        x1, y1 = xy1
        x2, y2 = xy2

        if x1 > x2:
            x1 = x1 + pixels
            x2 = x2 - pixels
        else:
            x1 = x1 - pixels
            x2 = x2 + pixels

        if y1 > y2:
            y1 = y1 - pixels
            y2 = y2 + pixels
        else:
            y1 = y1 + pixels
            y2 = y2 - pixels

        self._points = [(x1,y1), (x2,y2)]

    def shrink(self):
        ''' NOT IMPLEMENTED YET '''
        pass

    def getFillPoints(self):

        self.startPoint, self.endPoint = self._points
        '''  startPoint and endPoint: [(x1,y1) , (x2,y2)]  '''

        startPointX = int(self.startPoint[1])
        startPointY = int(self.startPoint[0])

        endPointX = int(self.endPoint[1])
        endPointY = int(self.endPoint[0])

        fillPoints = []

        if startPointX > endPointX:

            if startPointY > endPointY:

                for c in range(endPointY, startPointY + 1):
                    for i in range(endPointX, startPointX + 1):
                        fillPoints.append( (i, c) )
            else:
                for c in range(startPointY, endPointY + 1):
                    for i in range(endPointX, startPointX + 1):
                        fillPoints.append( (i, c) )

        else:

            if startPointY > endPointY:

                for c in range(endPointY, startPointY + 1):
                    for i in range(startPointX, endPointX + 1):
                        fillPoints.append( (i, c) )
            else:
                for c in range(startPointY, endPointY + 1):
                    for i in range(startPointX, endPointX + 1):
                        fillPoints.append( (i, c) )

        return fillPoints

class PolygonMask(Mask):
    ''' create a polygon mask '''

    def __init__(self, points, id, img_dim, negative = False):

        Mask.__init__(self, id, img_dim, 'polygon', negative)

        self._points = points

    def getFillPoints(self):

        proper_formatted_points = []
        yDim, xDim = self._img_dimension

        for each in self._points:
            proper_formatted_points.append(list(each))

        proper_formatted_points = np.array(proper_formatted_points)

        pb = polymask.Polygeom(proper_formatted_points)

        grid = np.mgrid[0:xDim,0:yDim].reshape(2,-1).swapaxes(0,1)

        inside = pb.inside(grid)

        p = np.where(inside==True)

        coords = polymask.getCoords(p, (yDim, xDim))

        return coords



def calcCenterCoords(img, selected_points, tune = True):
        ''' Determine center from coordinates on circle peferie.

            Article:
              I.D.Coope,
              "Circle Fitting by Linear and Nonlinear Least Squares",
              Journal of Optimization Theory and Applications vol 76, 2, Feb 1993
        '''

        numOfPoints = len(selected_points)

        B = []
        d = []

        for each in selected_points:
            x = each[0]
            y = each[1]

            B.append(x)                   # Build B matrix as vector
            B.append(y)
            B.append(1)

            d.append(x**2 + y**2)

        B = np.matrix(B)                  # Convert to numpy matrix
        d = np.matrix(d)

        B = B.reshape((numOfPoints, 3))   # Convert 1D vector to matrix
        d = d.reshape((numOfPoints, 1))

        Y = np.linalg.inv(B.T*B) * B.T * d   # Solve linear system of equations

        x_c = Y[0] / 2                    # Get x and r from transformation variables
        y_c = Y[1] / 2
        r = np.sqrt(Y[2] + x_c**2 + y_c**2)

        x_c = x_c.item()
        y_c = y_c.item()
        r = r.item()
        finetune_success = True

        if tune:
            newPoints = []

            for each in selected_points:
                x = each[0]
                y = each[1]

                optimPoint = finetuneAgbePoints(img, int(x_c), int(y_c), int(x), int(y), r)

                if optimPoint == False:
                    optimPoint = (x,y)
                    finetune_success = False

                newPoints.append(optimPoint)

            selected_points = newPoints
            xy, r = calcCenterCoords(img, selected_points, tune = False)
            x_c = xy[0]
            y_c = xy[1]

        if finetune_success == False:
#            wx.MessageBox('Remember to set the points "outside" the AgBe ring, a circle will then be fitted to the first found ring behind them.', 'Center search failed', wx.OK | wx.ICON_ERROR)
            raise SASExceptions.CenterNotFound('Fine tune center search failed')

        return ( (x_c, y_c), r )


def createMaskMatrix(img_dim, masks):
    ''' creates a 2D binary matrix of the same size as the image,
    corresponding to the mask pattern '''

    negmasks = []
    posmasks = []
    neg = False

    for each in masks:
        if each.isNegativeMask() == True:
            neg = True
            negmasks.append(each)
        else:
            posmasks.append(each)

    if neg:
        for each in posmasks:
            negmasks.append(each)

            masks = negmasks
        mask = np.zeros(img_dim)
    else:
        mask = np.ones(img_dim)

    maxy = mask.shape[1]
    maxx = mask.shape[0]

    for each in masks:
        fillPoints = each.getFillPoints()

        if each.isNegativeMask() == True:
            for eachp in fillPoints:
                if eachp[0] < maxx and eachp[0] >= 0 and eachp[1] < maxy and eachp[1] >= 0:
                    mask[eachp] = 1
        else:
            for eachp in fillPoints:
                if eachp[0] < maxx and eachp[0] >= 0 and eachp[1] < maxy and eachp[1] >= 0:
                    mask[eachp] = 0

    #Mask is flipped (older RAW versions had flipped image)
    mask = np.flipud(mask)

    return mask


def radialAverage(in_image, x_cin, y_cin, mask = None, readoutNoise_mask = None, dezingering = 0, dezing_sensitivity = 4.0):
    ''' Radial averaging. and calculation of readout noise from a readout noise mask.
        It also returns the errorbars assuming possion distributed data

        in_image :     Input image
        dim:           Image dimentions
        x_c, y_c :     (x_c, y_c) Center coordinate in the image (Pixels)
        q_range :      q_range specifying [low_q high_q]

    '''

    in_image = np.float64(in_image)

    ylen, xlen = in_image.shape

    xlen = np.int(xlen)
    ylen = np.int(ylen)

    # If no mask is given, the mask is pure ones
    if mask == None:
        mask = np.ones(in_image.shape)

    if readoutNoise_mask == None:
        readoutNoiseFound = 0
        readoutNoise_mask = np.zeros(in_image.shape, dtype = np.float64)
    else:
        readoutNoiseFound = 1

    readoutN = np.zeros((1,4), dtype = np.float64)

    # Find the maximum distance to the edge in the image:
    maxlen1 = int(max(xlen - x_cin, ylen - y_cin, xlen - (xlen - x_cin), ylen - (ylen - y_cin)))

    diag1 = int(np.sqrt((xlen-x_cin)**2 + y_cin**2))
    diag2 = int(np.sqrt((x_cin**2 + y_cin**2)))
    diag3 = int(np.sqrt((x_cin**2 + (ylen-y_cin)**2)))
    diag4 = int(np.sqrt((xlen-x_cin)**2 + (ylen-y_cin)**2))

    maxlen = int(max(diag1, diag2, diag3, diag4, maxlen1))

    #print diag1, diag2, diag3, diag4, maxlen1

    # we set the "q_limits" (in pixels) so that it does radial avg on entire image (maximum qrange possible).
    q_range = (0, maxlen)

    ##############################################
    # Reserving memory for radial averaged output:
    ##############################################
    hist = np.zeros(q_range[1], dtype = np.float64)
    hist_count = np.zeros((3,q_range[1]), dtype = np.float64)  # -----" --------- for number of pixels in a circle at a certain q

    qmatrix = np.zeros((q_range[1], 4*xlen), dtype = np.float64)

    low_q = q_range[0]
    high_q = q_range[1]

    # This code is faulty.. x has been switched with y
    x_c = float(y_cin)
    y_c = float(x_cin)

    xlen_1 = ylen
    ylen_1 = xlen

    # print np.any(hist<0)
    # print np.any(hist<0)

    print 'Radial averaging in progress...',

    if RAWGlobals.compiled_extensions:
        ravg_ext.ravg(readoutNoiseFound,
                       readoutN,
                       readoutNoise_mask,
                       xlen_1, ylen_1,
                       x_c, y_c,
                       hist,
                       low_q, high_q,
                       in_image,
                       hist_count, mask, qmatrix, dezingering, dezing_sensitivity)
    else:
        ravg_python(readoutNoiseFound,
                       readoutN,
                       readoutNoise_mask,
                       xlen_1, ylen_1,
                       x_c, y_c,
                       hist,
                       low_q, high_q,
                       in_image,
                       hist_count, mask, qmatrix, dezingering, dezing_sensitivity)

    print 'done'

    # print np.any(hist<0)
    # print np.any(hist_count<0)

    hist_cnt = hist_count[2,:]    #contains x-mean

    hist_count = hist_count[0,:]  #contains N

    std_i = np.sqrt(hist_cnt/hist_count)

    std_i[np.where(np.isnan(std_i))] = 0

    iq = hist / hist_count

    # print iq

    if x_c > 0 and x_c < xlen and y_c > 0 and y_c < ylen:
        iq[0] = in_image[round(x_c), round(y_c)]  #the center is not included in the radial average, so it is set manually her


    #Estimated Standard deviation   - equal to the std of pixels in the area / sqrt(N)
    errorbars = std_i / np.sqrt(hist_count)


    if readoutNoiseFound:
        #Average readoutNoise
        readoutNoise = readoutN[0,1] /  readoutN[0,0]   ## sum(img(x,y)) / N
        print 'Readout Noise: ', readoutNoise

        #Estimated Standard deviation   - equal to the std of pixels in the area / sqrt(N)
        std_n = np.sqrt(readoutN[0,3] / readoutN[0,0])    # sqrt((X-MEAN)/N)
        errorbarNoise = std_n / np.sqrt(readoutN[0,0])

        print 'Readout Noise Err: ', errorbarNoise

        #Readoutnoise average subtraction
        iq = iq - readoutNoise
        errorbars = np.sqrt(np.power(errorbars, 2) + np.power(errorbarNoise, 2))


    iq[np.where(np.isnan(iq))] = 0
    errorbars[np.where(np.isnan(errorbars))] = 1e-10

    q = np.linspace(0, len(iq)-1, len(iq))

    if dezingering == 1:
        iq, errorbars = getIntensityFromQmatrix(qmatrix)
        iq[np.where(np.isnan(iq))] = 0
        errorbars[np.where(np.isnan(errorbars))] = 1e-10

    #Trim trailing zeros
    # iq = np.trim_zeros(iq, 'b')
    iq = iq[:-5]        #Last points are usually garbage they're very few pixels
                        #Cutting the last 5 points here.
    q = q[0:len(iq)]
    errorbars = errorbars[0:len(iq)]

    return [iq, q, errorbars, qmatrix]


def ravg_python(readoutNoiseFound, readoutN, readoutNoise_mask, xlen, ylen, x_c,
                y_c, hist, low_q, high_q, in_image, hist_count, mask, qmatrix,
                dezingering, dezing_sensitivity):

    WINDOW_LENGTH=30

    # double rel_x, rel_y, r, delta, deltaN, qmat_cnt, std;
    # int i, x, y, half_window_size, window_start_idx, q_idx, point_idx;

    window = np.empty(WINDOW_LENGTH)
    # double median;
    # double *window_ptr, *data;
    # int half_win_len;
    # int hist_length;

    hist_length = len(hist)

    data = qmatrix            #/* Pointer to the numpy array version of qmatrix */

    half_window_size = (WINDOW_LENGTH / 2.0)
    win_len = WINDOW_LENGTH

    for x in range(xlen):
            for y in range(ylen):
           # {
                rel_x = x-x_c
                rel_y = y_c-y

                r = int(((rel_y)**2. + (rel_x)**2.)**0.5)

                # //res(x,y) = r;

                if r < high_q and r > low_q and mask[x,y] == 1: #// && in_image(x,y) > 0)
                # {
                    q_idx = r

                    # /* res2(x,y) = r; */                                  /*  A test image, gives the included range image */

                    hist[r] = hist[r] + in_image[x,y]                    #/* Integration of pixel values */

                    qmat_cnt = hist_count[0, q_idx]                      #/* Number of pixels in a bin */
                    qmatrix[q_idx, int(qmat_cnt)] = in_image[x,y]        #/* Save pixel value for later analysis */

                    hist_count[0, q_idx] = hist_count[0, q_idx] + 1      #/* Number of pixels in a bin */

                    delta = in_image[x,y] - hist_count[1, q_idx]         #/* Calculation of variance start */

                    hist_count[1, q_idx] = hist_count[1, q_idx] + (delta / hist_count[0, q_idx])
                    hist_count[2, q_idx] = hist_count[2, q_idx] + (delta * (in_image[x,y]-hist_count[1, q_idx]))


                    # /* *******************   Dezingering   ******************** */

                    if hist_count[0, r] >= WINDOW_LENGTH and dezingering == 1:
                    # {
                        point_idx = int(hist_count[0, q_idx])
                        window_start_idx = point_idx - win_len

                        data = qmatrix[q_idx, window_start_idx:point_idx]

                        window = data

                        # window.sort()

                        # moveDataToWindow(window_ptr, data, win_len)
                        # quicksort(window, 0, WINDOW_LENGTH-1)

                        std = np.std(window)
                        # window_ptr = window                                                        #/* Reset pointers */
                        median = np.median(window)

                        # //printf("median: %f\\n", median);

                        half_win_len = point_idx - half_window_size


                        if qmatrix[q_idx, half_win_len] > (median + (dezing_sensitivity * std)): #{
                            qmatrix[q_idx, half_win_len] = median
                        # }

                    # } // end dezinger if case


                # }

                if readoutNoiseFound == 1 and r < high_q-1 and r > low_q and readoutNoise_mask[x,y] == 0:
                # {
                    readoutN[0,0] = readoutN[0,0] + 1
                    readoutN[0,1] = readoutN[0,1] + in_image[x,y]

                    deltaN = in_image[x,y] - readoutN[0,2]
                    readoutN[0,2] = readoutN[0,2] + (deltaN / readoutN[0,0]) #Running average
                    readoutN[0,3] = readoutN[0,3] + (deltaN * (in_image[x,y]-readoutN[0,2]))

                # }
            # }

    # /* *********************************************  */
    # /* Remove zingers at the first (window/2) points  */
    # /* *********************************************  */

    if dezingering == 1:

        half_window_size = int(WINDOW_LENGTH / 2.0)
        win_len = WINDOW_LENGTH

        for q_idx in range(hist_length):
            # {

            if hist_count[0, q_idx] > (win_len + half_window_size):
                # {

                for i in range(win_len+half_window_size, win_len, -1):
                    # {

                    point_idx = i
                    window_start_idx = point_idx - win_len

                    data = qmatrix[q_idx, window_start_idx:point_idx]

                    window = data

                    # moveDataToWindow(window_ptr, data, win_len)
                    # quicksort(window, 0, WINDOW_LENGTH-1)

                    std = np.std(window)
                    # window_ptr = window                      #/* Reset pointers */
                    median = np.median(window)

                    half_win_len = point_idx - win_len

                    if qmatrix[q_idx, half_win_len] > (median + (dezing_sensitivity * std)):
                    # {
                        qmatrix[q_idx, half_win_len] = median

                    # }
                # }
            # }
        # }
    # }

    print "\n\n********* Radial Averaging and dezingering ********\n"
    print "Done!"

    # """

    # code2 = """