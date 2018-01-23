
from reproject import reproject_interp, reproject_exact
from astropy.io import fits

import math
from ibcommon.strings import id_generator
# from qb_common.strings import id_generator

import operator
import os
import traceback
# from shutil import rmtree
# from shutil import copyfile

# import MySQLdb

import astrolib.coords as coords
import numpy as np
# import qb_common.mysqls as qbmysqls
from coordtrans import pix2coord, coord2pix, coordConv, natCoords
from astropy.io.fits import getheader, getdata, writeto, setval
from pyraf import iraf
from ibexternal import sextractor
from scipy import spatial

# from scipy import stats
# import montage_wrapper as montage
# from montages import rebinImage

__author__ = "ibojicic"
__date__ = "$Oct 26, 2015 2:52:09 PM$"


def reprojectOnTemplate(inImage, outImage, templateImage):
    """
    Reproject one image based on image template
    :param templateImage:
    :param inImage: full path to the input image (fits)
    :param outImage: full path to the output image (fits)
    """
    hdu1 = fits.open(templateImage)[0]
    hdu2 = fits.open(inImage)[0]
    array, footprint = reproject_interp(hdu2, hdu1.header)
    fits.writeto(outImage, array, hdu1.header, clobber=True)
    return True


def doSubtraction(leftimage, rightmage, subimage, ratio=-1, lefthdrNo=0, righthdrNo=0, mindist=4, alignment=False):
    rnd = id_generator()
    ## TODO create temp files to be destroyed in case script dies
    temp_left = "/tmp/temp_left" + rnd + ".fits"
    temp_right = "/tmp/temp_right" + rnd + ".fits"

    if alignment:
        align_images(rightmage, leftimage, temp_right, temp_left)
    else:
        shiftto0header(leftimage, temp_left, lefthdrNo)
        shiftto0header(rightmage, temp_right, righthdrNo)

    tempimage = "/tmp/tmprc" + rnd + ".fits"
    temp_sub = "/tmp/temp_sub" + rnd + ".fits"

    if os.path.exists(tempimage):
        os.remove(tempimage)
    if os.path.exists(temp_sub):
        os.remove(temp_sub)

    if ratio == -1:
        ratios = findRatioForSubtracion(temp_left, temp_right, mindist)
        ratio = ratios['AVG']

    correctImageByFactor(temp_right, tempimage, ratio)
    contsubstrImage(temp_left, tempimage, subimage, 'y')

    if os.path.exists(tempimage):
        os.remove(tempimage)
    if os.path.exists(temp_left):
        os.remove(temp_left)
    if os.path.exists(temp_right):
        os.remove(temp_right)

    return True


def checkValidVal(inImage, X, Y, minval, maxval, boxsize=3):
    box = range(0 - boxsize, 1 + boxsize)
    totpix = 0
    try:
        for xtmp in box:
            for ytmp in box:
                listline = "[%s:%s,%s:%s]" % (int(X + xtmp), int(X + xtmp), int(Y + ytmp), int(Y + ytmp))
                checkpixtmp = iraf.listpix(inImage + listline, Stdout=1)
                checkpix = checkpixtmp[0].split()
                totpix += float(checkpix[1])
        result = (minval < (totpix / pow(2 * boxsize + 1, 2)) < maxval) and totpix != 0.0
    except Exception:
        # traceback.print_exc()  # print full exception
        result = False
    return result


def align_images(inImage_1, inImage_2, outImage_1, outImage_2):
    rnd_1 = id_generator()
    rnd_2 = id_generator()

    corr_image_1 = "/tmp/corr_image" + rnd_1 + ".fits"
    corr_image_2 = "/tmp/corr_image" + rnd_2 + ".fits"
    try:
        shiftto0header(inImage_1, corr_image_1, 1)
        shiftto0header(inImage_2, corr_image_2, 1)
    except:
        corr_image_1 = inImage_1
        corr_image_2 = inImage_2

    blank_image_1 = "/tmp/temp_image" + rnd_1 + ".fits"
    blank_image_2 = "/tmp/temp_image" + rnd_2 + ".fits"

    iraf.imarith(operand1=corr_image_1, operand2=0, op='*', result=blank_image_1)
    iraf.imarith(operand1=corr_image_2, operand2=0, op='*', result=blank_image_2)

    set_1 = "{},{}".format(corr_image_1, blank_image_2)
    set_2 = "{},{}".format(corr_image_2, blank_image_1)

    iraf.imcombine(set_1, outImage_1, combine='sum', offset='wcs')
    iraf.imcombine(set_2, outImage_2, combine='sum', offset='wcs')

    if os.path.exists(corr_image_1):
        os.remove(corr_image_1)
    if os.path.exists(corr_image_2):
        os.remove(corr_image_2)

    if os.path.exists(blank_image_1):
        os.remove(blank_image_1)
    if os.path.exists(blank_image_2):
        os.remove(blank_image_2)


def shiftto0header(inImage, outImage, hdrNo):
    hdr_new = getheader(inImage, hdrNo)
    data_new = getdata(inImage, hdrNo)
    writeto(outImage, data_new, hdr_new)
    return


def findRatioForSubtracion(haimage, rimage, mindist=4.):
    ha_sex = runSex(haimage)
    r_sex = runSex(rimage)

    average_array = np.array([])
    positions_array = []

    ha_positions = [[ha.get('X_IMAGE'), ha.get('Y_IMAGE')] for ha in ha_sex]
    r_positions = [[r.get('X_IMAGE'), r.get('Y_IMAGE')] for r in r_sex]
    results = spatial.KDTree(ha_positions).query_ball_tree(spatial.KDTree(r_positions), mindist)

    for key, ha_res in enumerate(ha_sex):
        if results[key]:
            for r_res_key in results[key]:
                r_res = r_sex[r_res_key]
                try:
                    ratio = ha_res['FLUX_BEST'] / r_res['FLUX_BEST']
                    average_array = np.insert(average_array, 0, ratio)
                    pair = {'ha': ha_res, 'r': r_res}
                    # makeKernel(haimage,rimage,pair)
                    positions_array.append(pair)
                except Exception:
                    continue

    avg = np.average(average_array)
    med = np.median(average_array)
    std = np.std(average_array)

    return {'AVG': avg, 'STD': std, 'MED': med, 'positions': positions_array}


def correctImageByFactor(inImage, outImage, factor, hdrNo=0):
    rnd = id_generator()
    temp_image = "/tmp/temp_image" + rnd + ".fits"
    shiftto0header(inImage, temp_image, hdrNo)

    if os.path.exists(outImage):
        os.remove(outImage)
    iraf.imarith(operand1=temp_image, operand2=factor, op='*', result=outImage)
    if os.path.exists(temp_image):
        os.remove(temp_image)

    return


def contsubstrImage(left, right, result, allowshift='y', lefthdrNo=0, righthdrNo=0):
    if os.path.exists(result):
        os.remove(result)

    hdr_dvd = getheader(left, lefthdrNo)
    maxX_dvd = int(hdr_dvd['NAXIS1'])
    maxY_dvd = int(hdr_dvd['NAXIS2'])

    hdr_dvs = getheader(right, righthdrNo)
    maxX_dvs = int(hdr_dvs['NAXIS1'])
    maxY_dvs = int(hdr_dvs['NAXIS2'])

    if maxX_dvd == maxX_dvs and maxY_dvd == maxY_dvs:
        iraf.imarith(operand1=left, operand2=right, op='-', result=result)
    elif allowshift == 'y':
        maxX = min(maxX_dvd, maxX_dvs)
        maxY = min(maxY_dvd, maxY_dvs)
        cutout = '[1:%s,1:%s]' % (maxX, maxY)
        iraf.imarith(operand1=left + cutout, operand2=right + cutout, op='-', result=result)
    return


# input Xcoord and Ycoord RA (degrees) and DEC (degrees)
# boxCoords in arcsec
def FITScutout(inImage, outImage, Xcoord, Ycoord, box_x_coord, box_y_coord = None, incoords='', hdrNo=0, extracomm="",
               epoch="J2000"):
    if os.path.exists(outImage):
        os.remove(outImage)
    if box_y_coord is None:
        box_y_coord = box_x_coord

    if incoords == '':

        hdr = getheader(inImage, hdrNo)
        max_x = int(hdr['NAXIS1'])
        max_y = int(hdr['NAXIS2'])

        pos_x, pos_y = natCoords(inImage, Xcoord, Ycoord)

        if epoch == "1950":
            temp_coord = coords.Position((pos_x, pos_y))  # ,equinox='b1950')
            calc_coord = temp_coord.b1950()
            pos_x = calc_coord[0]
            pos_y = calc_coord[1]

        temp_offset_x_1, temp_offset_y_1 = coord2pix(inImage, pos_x + float(box_x_coord) / 3600., pos_y)
        temp_offset_x_2, temp_offset_y_2 = coord2pix(inImage, pos_x, float(box_y_coord) / 3600. + pos_y)

        coord_x_center, coord_y_center = coord2pix(inImage, pos_x, pos_y)

        if extracomm == 'msxcheckoffset':
            if coord_x_center > 216000:
                coord_x_center = coord_x_center - 216000

        coord_offset = max(abs(coord_x_center - temp_offset_x_1), abs(coord_y_center - temp_offset_y_1),
                           abs(coord_x_center - temp_offset_x_2), abs(coord_y_center - temp_offset_y_2))

        coord_x_min = coord_x_center - coord_offset
        coord_y_min = coord_y_center - coord_offset

        coord_x_max = coord_x_center + coord_offset
        coord_y_max = coord_y_center + coord_offset

        if coord_x_center < 0 or coord_x_center > max_x or coord_y_center < 0 or coord_y_center > max_y:
            return 0

        if coord_x_max >= max_x:
            coord_x_max = max_x - 1
        if coord_y_max >= max_y:
            coord_y_max = max_y - 1
        if coord_x_min <= 0:
            coord_x_min = 1
        if coord_y_min <= 0:
            coord_y_min = 1

        imcopy_coords = '[{}:{},{}:{}]'.format(int(coord_x_min), int(coord_x_max), int(coord_y_min), int(coord_y_max))
    else:
        imcopy_coords = incoords

    imcopy_line = inImage + imcopy_coords
    try:
        iraf.imcopy(imcopy_line, outImage, verbose="no")
        new_hdr = getheader(outImage, hdrNo)
        if 'CRVAL1' not in new_hdr:
            iraf.hedit(outImage, fields='CRVAL1', value=0., add='yes', verify='no')
        if 'CRVAL2' not in new_hdr:
            iraf.hedit(outImage, fields='CRVAL2', value=0., add='yes', verify='no')
    except:
        traceback.print_exc()  # print full exception
        print 'Prooblem'
        exit()

    return imcopy_coords

def valueAtCoord(inImage, X, Y):
    hdu = fits.open(inImage)
    data = hdu[0].data
    return data[:, Y, X]

def updateHeader(inImage, newvals):
    hdu = fits.open(inImage)
    header = hdu[0].header
    for key,val in newvals.iteritems():
        header[key] = val
    hdu.writeto(inImage, clobber=True)

def getHeaderItems(inImage, items, hdrNo=0):
    hdr = getheader(inImage, hdrNo)
    if not isinstance(items,list):
        items = [items]
    keys = set(items).intersection(hdr)
    result = {key:hdr[key] for key in keys}
    return result


def runSex(inImage):
    # Create a SExtractor instance
    sex = sextractor.SExtractor()

    # Modify the SExtractor configuration

    #    sex.config['CATALOG_NAME'] =     'test.ASC'
    #    sex.config['ASSOC_NAME'] =       'testcatalog.csv'
    sex.config['ASSOC_PARAMS'] = '1,2'
    sex.config['ASSOC_RADIUS'] = 3
    sex.config['PARAMETERS_NAME'] = 'default.param'
    sex.config['CATALOG_TYPE'] = 'ASCII_HEAD'
    sex.config['DETECT_MINAREA'] = 4
    sex.config['THRESH_TYPE'] = "RELATIVE"
    sex.config['DETECT_THRESH'] = 3.
    sex.config['ANALYSIS_THRESH'] = 3.
    sex.config['FILTER'] = 'Y'
    # sex.config['FILTER_NAME'] = '/usr/local/star-namaka/bin/extractor/config/default.conv'
    sex.config['DEBLEND_NTHRESH'] = 32
    sex.config['DEBLEND_MINCONT'] = 0.005
    sex.config['CLEAN'] = 'Y'
    sex.config['CLEAN_PARAM'] = 1.0
    sex.config['MAG_ZEROPOINT'] = 0.0
    sex.config['PHOT_APERTURES'] = 5
    sex.config['PHOT_AUTOPARAMS'] = 2.5, 3.5
    # sex.config['PHOT_PETROPARAMS'] = 2.0,3.5
    # sex.config['PHOT_FLUXFRAC'] =    0.5
    sex.config['MASK_TYPE'] = 'CORRECT'
    sex.config['DETECT_TYPE'] = 'CCD'
    sex.config['PIXEL_SCALE'] = 0.3
    sex.config['SATUR_LEVEL'] = 50000.0
    sex.config['GAIN'] = 1.0
    sex.config['MAG_GAMMA'] = 4.
    sex.config['SEEING_FWHM'] = 1.2
    # sex.config['STARNNW_NAME'] = '/usr/local/star-namaka/bin/extractor/config/default.nnw'
    sex.config['BACK_SIZE'] = 32
    sex.config['BACK_FILTERSIZE'] = 3
    sex.config['BACK_TYPE'] = 'AUTO'
    sex.config['BACK_VALUE'] = 0.0
    sex.config['BACKPHOTO_TYPE'] = 'LOCAL'
    sex.config['BACKPHOTO_THICK'] = 12
    sex.config['CHECKIMAGE_TYPE'] = 'APERTURES'
    sex.config['CHECKIMAGE_NAME'] = 'check.fits'
    sex.config['MEMORY_OBJSTACK'] = 2000
    sex.config['MEMORY_PIXSTACK'] = 300000
    sex.config['MEMORY_BUFSIZE'] = 1024
    sex.config['VERBOSE_TYPE'] = 'NORMAL'

    # Add a parameter to the parameter list

    sex.config['PARAMETERS_LIST'].append('BACKGROUND')
    sex.config['PARAMETERS_LIST'].append('THRESHOLD')
    sex.config['PARAMETERS_LIST'].append('CLASS_STAR')
    # sex.config['PARAMETERS_LIST'].append('VECTOR_ASSOC')
    # sex.config['PARAMETERS_LIST'].append('NUMBER_ASSOC')
    sex.config['PARAMETERS_LIST'].append('FLUX_RADIUS')
    sex.config['PARAMETERS_LIST'].append('ELLIPTICITY')

    # Lauch SExtractor on a FITS file
    sex.run(inImage)

    # Read the resulting catalog [first method, whole catalog at once]
    catalog = sex.catalog()

    # total_stars = len(catalog)

    catalog.sort(key=operator.itemgetter('FWHM_IMAGE'))

    # for star in catalog:

    # star = catalog[0]

    # xcrd = int(star['X_IMAGE'])

    return catalog

#
#
#
#
#
# # def regridMIR(inimage_reg, inimage_mod, outimage):
# #
# #     from mirpy import miriad
# #     tmp_reg = "/tmp/" + id_generator() + 'tmpreg.mir'
# #     tmp_mod = "/tmp/" + id_generator() + 'tmpmod.mir'
# #     tmp_out = "/tmp/" + id_generator() + 'tmpout.mir'
# #
# #     miriad.fits(IN=inimage_reg, out=tmp_reg, op='xyin')
# #     miriad.fits(IN=inimage_mod, out=tmp_mod, op='xyin')
# #     miriad.regrid(IN=tmp_reg, axes='1,2', out=tmp_out, tin=tmp_mod)
# #     miriad.fits(IN=tmp_out, out=outimage, op='xyout')
# #
# #     if (os.path.exists(tmp_reg)): rmtree(tmp_reg)
# #     if (os.path.exists(tmp_mod)): rmtree(tmp_mod)
# #     if (os.path.exists(tmp_out)): rmtree(tmp_out)
#
#
# def reprojectMIR(inimage_reg, outimage):
#     from mirpy import miriad
#
#     tmp_reg = "/tmp/" + id_generator() + 'tmpreg.mir'
#     tmp_out = "/tmp/" + id_generator() + 'tmpout.mir'
#     miriad.fits(IN=inimage_reg, out=tmp_reg, op='xyin')
#     miriad.regrid(IN=tmp_reg, axes='1,2', out=tmp_out, project='SIN')
#     miriad.fits(IN=tmp_out, out=outimage, op='xyout')
#     if (os.path.exists(tmp_reg)): rmtree(tmp_reg)
#     if (os.path.exists(tmp_out)): rmtree(tmp_out)
#
#
# # COPIED CHECK AND DELETE
#
#
# def cutout_size(inImage, hdrNo=0):
#     hdr = getheader(inImage, hdrNo)
#
#     coord_type = hdr['CTYPE1'].split('-')[0]
#
#     pix1 = int(round(float(hdr['NAXIS1'])))
#     pix2 = int(round(float(hdr['NAXIS2'])))
#
#     maxX, maxY = pix2coord(inImage, pix1, pix2)
#
#     minX, minY = pix2coord(inImage, 0, 0)
#
#     if coord_type == 'GLON':
#         newMax = coordConv(maxX, maxY, 'galactic j2000', 'fk5 j2000', 'glonglat', 'radec')
#         maxX = float(newMax[0])
#         maxY = float(newMax[1])
#
#         newMin = coordConv(minX, minY, 'galactic j2000', 'fk5 j2000', 'glonglat', 'radec')
#         minX = float(newMin[0])
#         minY = float(newMin[1])
#
#     pos2 = coords.Position((minX, minY))
#     pos3 = coords.Position((maxX, minY))
#     pos4 = coords.Position((minX, maxY))
#
#     Xsep = str(pos2.angsep(pos3)).split(' ')
#     Ysep = str(pos2.angsep(pos4)).split(' ')
#
#     if Xsep[1] == 'degrees':
#         X = float(Xsep[0]) * 3600.
#         Y = float(Ysep[0]) * 3600.
#     else:
#         X = -1.
#         Y = -1.
#
#     return X, Y
#
#
# def doGalexG(ndimage, fdimage, outimage):
#     if (os.path.exists(outimage)): os.remove(outimage)
#     tmpim1 = "/tmp/tmpsmooth" + id_generator() + ".fits"
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#     iraf.gauss(input=fdimage, output=tmpim1, sigma=5)
#     contAvgImage(ndimage, tmpim1, outimage)
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#
#
# def contAvgImage(left, right, result, allowshift='y', lefthdrNo=0, righthdrNo=0):
#     if (os.path.exists(result)): os.remove(result)
#     tmpim1 = "/tmp/tmpsqrt1" + id_generator() + ".fits"
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#
#     hdr_dvd = getheader(left, lefthdrNo)
#     maxX_dvd = int(hdr_dvd['NAXIS1'])
#     maxY_dvd = int(hdr_dvd['NAXIS2'])
#
#     hdr_dvs = getheader(right, righthdrNo)
#     maxX_dvs = int(hdr_dvs['NAXIS1'])
#     maxY_dvs = int(hdr_dvs['NAXIS2'])
#     if maxX_dvd == maxX_dvs and maxY_dvd == maxY_dvs:
#         iraf.imarith(operand1=left, operand2=right, op='+', result=tmpim1)
#         iraf.imarith(operand1=tmpim1, operand2="2", op='/', result=result)
#     elif allowshift == 'y':
#         maxX = min(maxX_dvd, maxX_dvs)
#         maxY = min(maxY_dvd, maxY_dvs)
#         cutout = '[1:%s,1:%s]' % (maxX, maxY)
#         iraf.imarith(operand1=left + cutout, operand2=right + cutout, op='-', result=result)
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#
#     return
#
#
# def quotientImage(dividend, divisor, result, allowshift='y', dvdhdrNo=0, dvsrhdrNo=0, do_align=False):
#     if (os.path.exists(result)): os.remove(result)
#
#     rnd = id_generator()
#
#     temp_dividend = "/tmp/temp_dividend" + rnd + ".fits"
#     temp_divisor = "/tmp/temp_divisor" + rnd + ".fits"
#
#     if do_align:
#         align_images(dividend, divisor, temp_dividend, temp_divisor)
#     else:
#         shiftto0header(dividend, temp_dividend, dvdhdrNo)
#         shiftto0header(divisor, temp_divisor, dvsrhdrNo)
#
#     hdr_dvd = getheader(temp_dividend, 0)
#     maxX_dvd = int(hdr_dvd['NAXIS1'])
#     maxY_dvd = int(hdr_dvd['NAXIS2'])
#
#     hdr_dvs = getheader(temp_divisor, 0)
#     maxX_dvs = int(hdr_dvs['NAXIS1'])
#     maxY_dvs = int(hdr_dvs['NAXIS2'])
#
#     if maxX_dvd == maxX_dvs and maxY_dvd == maxY_dvs:
#         iraf.imarith(operand1=temp_dividend, operand2=temp_divisor, op='/', result=result)
#     elif allowshift == 'y':
#         maxX = min(maxX_dvd, maxX_dvs)
#         maxY = min(maxY_dvd, maxY_dvs)
#         cutout = '[1:%s,1:%s]' % (maxX, maxY)
#         iraf.imarith(operand1=temp_dividend + cutout, operand2=temp_divisor + cutout, op='/', result=result)
#
#     return
#
#
#
# def correctMissAllign(inImage, outImage, positionsArray):
#     window_size = 10
#
#     hdulist = fits.open(inImage)
#     scidata = hdulist[0].data
#
#     for found_stars in positionsArray:
#
#         ha_star = found_stars['ha']
#
#         window = scidata[
#                  int(ha_star['Y_IMAGE']) - window_size: int(ha_star['Y_IMAGE']) + window_size,
#                  int(ha_star['X_IMAGE']) - window_size: int(ha_star['X_IMAGE']) + window_size
#                  ]
#
#         avg = np.average(window)
#         med = np.median(window)
#         std = np.std(window)
#
#         print avg, med, std
#
#         window_shape = window.shape
#
#         for i in range(window_size * window_size / 2):
#             try:
#                 max_key = np.argmax(window)
#                 max_index = np.unravel_index(max_key, window_shape)
#                 max_val = np.amax(window)
#                 min_key = np.argmin(window)
#                 min_index = np.unravel_index(min_key, window_shape)
#                 min_val = np.amin(window)
#
#                 if min_val > 0. or max_val < 0.:
#                     break
#                 new_val = max_val + min_val
#                 window[min_index[0]][min_index[1]] = new_val
#                 window[max_index[0]][max_index[1]] = new_val
#             except:
#                 continue
#
#         scidata[
#         int(ha_star['Y_IMAGE']) - 10: int(ha_star['Y_IMAGE']) + 10,
#         int(ha_star['X_IMAGE']) - 10: int(ha_star['X_IMAGE']) + 10
#         ] = window
#
#     hdulist.writeto(outImage)
#
#
#
#
#
#
# def addImageByFactor(inImage, outImage, factor, hdrNo):
#     # rnd = id_generator()
#     # temp_image = "/tmp/temp_image" + rnd + ".fits"
#     # shiftto0header(inImage, temp_image, hdrNo)
#
#     if (os.path.exists(outImage)): os.remove(outImage)
#     iraf.imarith(operand1=inImage, operand2=factor, op='+', result=outImage)
#     return
#
#
#
#
# def contsqrtImage(left, right, result, allowshift='y', lefthdrNo=0, righthdrNo=0):
#     if (os.path.exists(result)): os.remove(result)
#     tmpim1 = "/tmp/tmpsqrt1" + id_generator() + ".fits"
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#
#     hdr_dvd = getheader(left, lefthdrNo)
#     maxX_dvd = int(hdr_dvd['NAXIS1'])
#     maxY_dvd = int(hdr_dvd['NAXIS2'])
#
#     hdr_dvs = getheader(right, righthdrNo)
#     maxX_dvs = int(hdr_dvs['NAXIS1'])
#     maxY_dvs = int(hdr_dvs['NAXIS2'])
#
#     if maxX_dvd == maxX_dvs and maxY_dvd == maxY_dvs:
#         iraf.imarith(operand1=left, operand2=right, op='*', result=tmpim1)
#         iraf.imfunction(input=tmpim1, output=result, function="sqrt")
#     elif allowshift == 'y':
#         maxX = min(maxX_dvd, maxX_dvs)
#         maxY = min(maxY_dvd, maxY_dvs)
#         cutout = '[1:%s,1:%s]' % (maxX, maxY)
#         iraf.imarith(operand1=left + cutout, operand2=right + cutout, op='-', result=result)
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#
#     return
#
#
# def addIntImage(left, right, result, allowshift='y', lefthdrNo=0, righthdrNo=0):
#     if (os.path.exists(result)): os.remove(result)
#     tmpim1 = "/tmp/tmpsqrt1" + id_generator() + ".fits"
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#
#     hdr_dvd = getheader(left, lefthdrNo)
#     maxX_dvd = int(hdr_dvd['NAXIS1'])
#     maxY_dvd = int(hdr_dvd['NAXIS2'])
#
#     hdr_dvs = getheader(right, righthdrNo)
#     maxX_dvs = int(hdr_dvs['NAXIS1'])
#     maxY_dvs = int(hdr_dvs['NAXIS2'])
#
#     if maxX_dvd == maxX_dvs and maxY_dvd == maxY_dvs:
#         iraf.imarith(operand1=left, operand2=right, op='*', result=tmpim1)
#         iraf.imfunction(input=tmpim1, output=result, function="sqrt")
#     elif allowshift == 'y':
#         maxX = min(maxX_dvd, maxX_dvs)
#         maxY = min(maxY_dvd, maxY_dvs)
#         cutout = '[1:%s,1:%s]' % (maxX, maxY)
#         iraf.imarith(operand1=left + cutout, operand2=right + cutout, op='*', result=tmpim1)
#         iraf.imfunction(input=tmpim1, output=result, function="sqrt")
#     if (os.path.exists(tmpim1)): os.remove(tmpim1)
#
#     return
#
#
# def decompressFitsImage(inFile, outPath, compressMethod, outPref='decom'):
#     newfile = "%s_%s.fits" % (outPref, id_generator())
#     command = "%s -O %s/%s %s" % (compressMethod, outPath.rstrip("/"), newfile, inFile)
#     os.system(command)
#     return newfile
#
#
# def compressFitsImage(inFile, outFile, outPath, compressMethod, outPref='decom'):
#     command = "%s -O %s/%s %s" % (compressMethod, outPath.rstrip("/"), newfile, inFile)
#     os.system(command)
#     return newfile
#
#
# def parse_header(in_image, return_fields, header_no=0):
#     return_list = {}
#     header_no = int(header_no)
#     if header_no == -1:
#         header_no = 0
#     for i in range(header_no + 1):
#         try:
#             hdr = getheader(in_image, i)
#         except:
#             return return_list
#         for ind in return_fields:
#             if ind in hdr:
#                 return_list[ind] = hdr[ind]
#             elif ind not in return_list:
#                 return_list[ind] = "null"
#
#     return return_list
#
#
# def find_band(filename, bands, headers, band_header=False):
#     if len(bands) == 1: return bands[0]
#     if band_header != "null":
#         if band_header in headers: return headers[band_header].replace("#", "")
#     for band in bands:
#         if band in filename: return band
#     return "null"
#
#
# def corect_VISTA_headers(in_image, out_image, hdr_no=0):
#     # try:
#     hdr_new = getheader(in_image, int(hdr_no))
#     data_new = getdata(in_image, int(hdr_no))
#     writeto(out_image, data_new, hdr_new)
#     pix1 = int(round(float(hdr_new['NAXIS1']) / 2.))
#     pix2 = int(round(float(hdr_new['NAXIS2']) / 2.))
#     RA, DEC = pix2coord(out_image, pix1, pix2)
#     iraf.hedit(out_image, fields='RA_CENT', value=RA, add='yes', verify='no')
#     iraf.hedit(out_image, fields='DEC_CENT', value=DEC, add='yes', verify='no')
#     if os.path.exists(out_image): return 1
#     # except:
#     #    return 0
#
#
# def check_UKIDSS(in_image, hdr_no=0):
#     # database = MySQLdb.connect(host="niksicko", user="gpneadmin", passwd="(g.pne.admin)")
#     database = qbmysqls.hash_database()
#     cursor = database.cursor()
#     hdr = getheader(in_image, int(hdr_no))
#     sdsuid = str(hdr['SDSUID'])
#     runid = str(hdr['RUNID'])
#     ra = str(hdr['RA_CENT'])
#     dec = str(hdr['DEC_CENT'])
#     sql = "SELECT * FROM `ImagesSources`.`UKIDSS` \
#             WHERE `MainGPN`.`GPNspherDist_ib`(`RA_CENT`,`DEC_CENT`," + ra + "," + dec + ") * 3600 < 60 \
#             AND `RUNID` = '" + runid + "' AND `SDSUID` = '" + sdsuid + "';"
#     cursor.execute(sql)
#     ids = cursor.fetchall()
#     if ids:
#         return True
#     else:
#         return False
#
#
# def check_VVVE(in_image, in_hdr=0):
#     # database = MySQLdb.connect(host="niksicko", user="gpneadmin", passwd="(g.pne.admin)")
#     database = qbmysqls.hash_database()
#     cursor = database.cursor()
#     hdr = getheader(in_image, int(in_hdr))
#     this_object = str(hdr['OBJECT'])
#     vsa_mfid = str(hdr['VSA_MFID'])
#     ra = str(hdr['RA_CENT'])
#     dec = str(hdr['DEC_CENT'])
#     sql = "SELECT * FROM `ImagesSources`.`VVVE` \
#             WHERE `MainGPN`.`GPNspherDist_ib`(`RA_CENT`,`DEC_CENT`," + ra + "," + dec + ") * 3600 < 60 \
#             AND `OBJECT` = '" + this_object + "' AND `VSA_MFID` = '" + vsa_mfid + "';"
#     cursor.execute(sql)
#     ids = cursor.fetchall()
#     if ids:
#         return True
#     else:
#         return False
#
#
#
#
#
#
#
#
#
# def findPercentile(inImage, levels, dataset='all'):
#     result = {}
#     data_point = 2
#     if dataset == 'nvss': data_point = 4
#     try:
#         A = iraf.listpixels(inImage, Stdout=1)
#         B = [float(i.split()[data_point]) for i in A if
#              (float(levels['v']['min']) < float(i.split()[data_point]) < float(
#                  levels['v']['max']))]
#         B.sort()
#         result['min'] = stats.scoreatpercentile(B, float(levels['r']['min']))
#         result['max'] = stats.scoreatpercentile(B, float(levels['r']['max']))
#         result['max_amp'] = max(B)
#         result['min_amp'] = min(B)
#         result['n_points'] = len(B)
#         result['description'] = stats.describe(B)
#     except Exception as e:
#         return {str(e)}
#     return result
#
#
# def checkCoverage(inImage, Xcoord, Ycoord, box_x_coord, box_y_coord, hdrNo=0, extracomm="",
#                   epoch="J2000"):
#     hdr = getheader(inImage, hdrNo)
#     min_x = 0
#     max_x = int(hdr['NAXIS1'])
#     min_y = 0
#     max_y = int(hdr['NAXIS2'])
#
#     box_x_coord = float(box_x_coord) / 3600.
#     box_y_coord = float(box_y_coord) / 3600.
#
#     nat_cords = natCoords(inImage, Xcoord, Ycoord)
#     pos_x = nat_cords[0]
#     pos_y = nat_cords[1]
#
#     if epoch == "1950":
#         temp_coord = coords.Position((pos_x, pos_y))  # ,equinox='b1950')
#         calc_coord = temp_coord.b1950()
#         pos_x = calc_coord[0]
#         pos_y = calc_coord[1]
#
#     temp_offset_x_1, temp_offset_y_1 = coord2pix(inImage, pos_x + box_x_coord, pos_y)
#     temp_offset_x_2, temp_offset_y_2 = coord2pix(inImage, pos_x, box_y_coord + pos_y)
#
#     coord_x_center, coord_y_center = coord2pix(inImage, pos_x, pos_y)
#
#     if extracomm == 'msxcheckoffset':
#         if coord_x_center > 216000: coord_x_center = coord_x_center - 216000
#
#     coordoffset = max(abs(coord_x_center - temp_offset_x_1), abs(coord_y_center - temp_offset_y_1),
#                       abs(coord_x_center - temp_offset_x_2),
#                       abs(coord_y_center - temp_offset_y_2))
#
#     coordXmin = coord_x_center - coordoffset
#     coordYmin = coord_y_center - coordoffset
#
#     coordXmax = coord_x_center + coordoffset
#     coordYmax = coord_y_center + coordoffset
#
#     if coord_x_center < 0 or coord_x_center > max_x or coord_y_center < 0 or coord_y_center > max_y:
#         return 'n'
#
#     if coordXmax >= max_x or coordYmax >= max_y or coordXmin <= min_x or coordYmin <= min_y:
#         return 'n'
#
#     return 'y'
#
#
# # def findRatioForSubtracion(haimage, rimage, mindist = 4, maxobjects = 100):
# #     ha_sex = runSex(haimage)
# #     r_sex = runSex(rimage)
# #
# #     avarray = np.array([])
# #
# #     no_step = ceil(len(ha_sex) / maxobjects)
# #     avg = 2.
# #     med = 1
# #     std = 0.
# #     k = 6
# #     flagpass = -1
# #
# #     while (abs(avg - med) > std or flagpass == -1) and k > 5:
# #         # while (std > avg / 10. or flagpass == -1) and k > 5:
# #         k = 0
# # #        for ha_star in ha_sex:
# #         for ha_step in range(0,len(ha_sex),int(no_step)):
# #             ha_star = ha_sex[ha_step]
# #             if ha_star['FWHM_IMAGE'] > 3 and ha_star['FWHM_IMAGE'] < 7 and ha_star['ELLIPTICITY'] < 0.2 and ha_star[
# #                 'FLAGS'] == 0:
# #                 ha_x = ha_star['X_IMAGE']
# #                 ha_y = ha_star['Y_IMAGE']
# #                 for r_star in r_sex:
# #                     if r_star['FWHM_IMAGE'] > 3 and r_star['FWHM_IMAGE'] < 7 and r_star['ELLIPTICITY'] < 0.2 and r_star[
# #                         'FLAGS'] == 0:
# #                         r_x = r_star['X_IMAGE']
# #                         r_y = r_star['Y_IMAGE']
# #
# #                         dist = sqrt(pow(ha_x - r_x, 2) + pow(ha_y - r_y, 2))
# #
# #                         if dist < mindist:
# #
# #                             ratio = ha_star['FLUX_BEST'] / r_star['FLUX_BEST']
# #                             if flagpass == -1 or abs(ratio - avg) < std:
# #                                 k = k + 1
# #                                 avarray = np.insert(avarray, 0, ratio)
# #         flagpass = 1
# #         avg = np.average(avarray)
# #         med = np.median(avarray)
# #         std = np.std(avarray)
# #
# #     avg = {'AVG': avg, 'STD': std, 'MED': med}
# #     return avg
#
#
#
#
# def makeKernel(haimage, rimage, results, mindist=4):
#     print haimage
#     print rimage
#     exit()
#
#     window_size = 10
#
#     hdulist_ha = fits.open(haimage)
#     scidata_ha = hdulist_ha[0].data
#
#     hdulist_r = fits.open(rimage)
#     scidata_r = hdulist_r[0].data
#
#     window_ha = scidata_ha[
#                 int(results['ha']['Y_IMAGE']) - window_size: int(results['ha']['Y_IMAGE']) + window_size,
#                 int(results['ha']['X_IMAGE']) - window_size: int(results['ha']['X_IMAGE']) + window_size
#                 ]
#
#     window_r = scidata_r[
#                int(results['ha']['Y_IMAGE']) - window_size: int(results['ha']['Y_IMAGE']) + window_size,
#                int(results['ha']['X_IMAGE']) - window_size: int(results['ha']['X_IMAGE']) + window_size
#                ]
#
#     subtract = np.subtract(window_ha, window_r)
#
#     avg = np.average(subtract)
#     med = np.median(subtract)
#     std = np.std(subtract)
#
#     print 'kernel', avg, med, std
#     shift_window_ha = window_ha
#     shift_window_r = window_r
#
#     for axis in range(2):
#         for dist in range(mindist):
#             shift_window_ha = np.delete(shift_window_ha, 0, axis)
#
#             # ha_star = found_stars['ha']
#             #
#             # window = scidata[
#             #       int(ha_star['Y_IMAGE']) - window_size: int(ha_star['Y_IMAGE']) + window_size,
#             #       int(ha_star['X_IMAGE']) - window_size: int(ha_star['X_IMAGE']) + window_size
#             #       ]
#             #
#             # avg = np.average(window)
#             # med = np.median(window)
#             # std = np.std(window)
#             #
#             # print avg, med, std
#             #
#             # window_shape = window.shape
#             #
#             # for i in range(window_size * window_size / 2):
#             #     try:
#             #         max_key = np.argmax(window)
#             #         max_index = np.unravel_index(max_key,window_shape)
#             #         max_val = np.amax(window)
#             #         min_key = np.argmin(window)
#             #         min_index = np.unravel_index(min_key,window_shape)
#             #         min_val = np.amin(window)
#             #
#             #         if min_val > 0. or max_val < 0.:
#             #             break
#             #         new_val = max_val + min_val
#             #         window[min_index[0]][min_index[1]] = new_val
#             #         window[max_index[0]][max_index[1]] = new_val
#             #     except:
#             #         continue
#             #
#             # scidata[
#             # int(ha_star['Y_IMAGE']) - 10: int(ha_star['Y_IMAGE']) + 10,
#             # int(ha_star['X_IMAGE']) - 10: int(ha_star['X_IMAGE']) + 10
#             # ] = window
#
#             # hdulist.writeto(outImage)
#
#

# def correct_iphas_corner(inImage):
#     rnd = id_generator()
#     corner = "/tmp/corner" + rnd + ".fits"
#     iraf.imcopy("{}:[1:300,1:300] {}".format(inImage, corner))
#
# def correct_header(inImage, tmpImage, hdrKey, hdrVal, hdrNo = 0, hdrDel = False):
#     imHeader = getheader(inImage, hdrNo)
#
#     if hdrKey in imHeader:
#         copyfile(inImage,tmpImage)
#         iraf.hedit(tmpImage, fields= hdrKey, value=hdrVal, add='yes', verify='no')
#         return tmpImage
#     return inImage
#
#
# def correct_shs_WATs(inFile):
#     setval(inFile, 'WAT0_001', value='system=image')
#     setval(inFile, 'WAT1_001', value='wtype=tan axtype=ra')
#     setval(inFile, 'WAT2_001', value='wtype=tan axtype=dec')
#
#
# def rebinOnTemplate(inImage,outImage,templateImage,templateFolder = '/tmp/'):
#
#     tempF = "{}{}/".format(templateFolder,id_generator())
#     if not os.path.exists(tempF):
#         os.makedirs(tempF)
#     else:
#         return False
#     copyfile(templateImage,tempF + "template.fits")
#     tmpltbl = 'templatetable' + id_generator() + ".tbl"
#     montage.mImgtbl(tempF,tmpltbl)
#     tmplhdr = 'templatehdt' + id_generator() + ".hdr"
#     montage.mMakeHdr(tmpltbl,tmplhdr)
#     montage.mProject(inImage,outImage,tmplhdr)
#     if (os.path.exists(tempF)): rmtree(tempF)
#     if (os.path.exists(tmpltbl)): os.remove(tmpltbl)
#     if (os.path.exists(tmplhdr)): os.remove(tmplhdr)
#

# # input Xcoord and Ycoord RA (degrees) and DEC (degrees)
# # boxCoords in arcsec
# def FITScutoutMir(inImage, outImage, Xcoord, Ycoord, boxXcoord, boxYcoord, incoords='', hdrNo=0):
#     from mirpy import miriad
#
#     tmpnvssmir = "/tmp/" + id_generator() + 'tmpnvss.mir'
#     tmpnvsssub = "/tmp/" + id_generator() + 'tmpnvss_sub.mir'
#
#     if (os.path.exists(outImage)): os.remove(outImage)
#
#     hdr = getheader(inImage, hdrNo)
#
#     minX = 0
#     maxX = int(hdr['NAXIS1'])
#     minY = 0
#     maxY = int(hdr['NAXIS2'])
#
#     boxXcoord = float(boxXcoord) / 3600.
#     boxYcoord = float(boxYcoord) / 3600.
#
#     ntcords = natCoords(inImage, Xcoord, Ycoord)
#     decpos1 = ntcords[0]
#     decpos2 = ntcords[1]
#
#     offsetPixX, offsetPixY = coord2pix(inImage, decpos1, decpos2 + boxYcoord)
#
#     coordXcent, coordYcent = coord2pix(inImage, decpos1, decpos2)
#
#     coordoffset = abs(coordYcent - abs(offsetPixY))
#
#     coordXmin = coordXcent - coordoffset
#     coordYmin = coordYcent - coordoffset
#
#     coordXmax = coordXcent + coordoffset
#     coordYmax = coordYcent + coordoffset
#
#     if coordXcent < 0 or coordXcent > maxX or coordYcent < 0 or coordYcent > maxY:
#         print "outside of the image?"
#         return 0
#
#     if coordXmax >= maxX: coordXmax = maxX - 1
#     if coordYmax >= maxY: coordYmax = maxY - 1
#     if coordXmin <= minX: coordXmin = minX + 1
#     if coordYmin <= minY: coordYmin = minY + 1
#
#     imsubregion = 'boxes(%s,%s,%s,%s)' % (coordXmin, coordYmin, coordXmax, coordYmax)
#
#     miriad.fits(IN=inImage, out=tmpnvssmir, op='xyin')
#     miriad.imsub(IN=tmpnvssmir, out=tmpnvsssub, region=imsubregion)
#     miriad.fits(IN=tmpnvsssub, out=outImage, op='xyout')
#
#     if (os.path.exists(tmpnvssmir)): rmtree(tmpnvssmir)
#     if (os.path.exists(tmpnvsssub)): rmtree(tmpnvsssub)
#
#     return
