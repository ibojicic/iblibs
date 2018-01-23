
# import os
# from shutil import rmtree

from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
# from qb_common.strings import id_generator
import montages
from astropy.io.fits import getheader
import math


def radec2hmsdms(crd1, crd2, frame_in='fk5', frame_out='fk5'):
    """
    Convert decimal RA,DEC to hmsdms RA,DEC
    :param crd1: float RA 
    :param crd2: float DEC
    :param frame_in: string
    :param frame_out: string
    :return: strings
    """
    crd1 = float(crd1)
    crd2 = float(crd2)
    c = SkyCoord(ra=crd1 * u.degree, dec=crd2 * u.degree, frame=frame_in)
    c = change_frame(c, frame_in, frame_out)
    RA = c.ra.to_string('h', sep=':', alwayssign=False, precision=1, pad=True)
    DEC = c.dec.to_string(sep=':', alwayssign=False, precision=2, pad=True)
    return RA, DEC


def coord2pix(inImage, Xcoord, Ycoord):
    try:
        hdulist = fits.open(inImage)
        w = wcs.WCS(hdulist[0].header)
        Xpix, Ypix = w.wcs_world2pix([[Xcoord, Ycoord]], 1)[0].astype(int)
        Xpix = Xpix + 1
        Ypix = Ypix + 1
    except:
        Xpix, Ypix = oldcoord2pix(inImage, Xcoord, Ycoord)
    return Xpix, Ypix


def pix2coord(inImage, Xcoord, Ycoord, hduid=0):
    try:
        ra, dec = montages.pix2coord(inImage, Xcoord, Ycoord)
        if math.isnan(float(ra)) or math.isnan(float(dec)) or not isinstance(ra, float) or not isinstance(dec, float):
            ra, dec = oldpix2coord(inImage, Xcoord, Ycoord, hduid)
    except:
        ra, dec = oldpix2coord(inImage, Xcoord, Ycoord, hduid)
    return ra, dec


def coordConv(Xcrd, Ycrd, inSys, outSys, inFormat, outFormat):
    transformation = [0, 0]
    from pyraf.iraf import skyctran

    skyctran.setParam('input', 'STDIN')
    skyctran.setParam('output', 'STDOUT')
    skyctran.setParam('insystem', inSys)

    if inFormat == 'radec' or inFormat == 'glonglat':
        skyctran.setParam('ilngunits', 'degrees')
        skyctran.setParam('ilatunits', 'degrees')
        skyctran.setParam('ilngformat', '%11.7f')
        skyctran.setParam('ilatformat', '%11.7f')
    elif inFormat == 'hmsdms':
        skyctran.setParam('ilngunits', 'hours')
        skyctran.setParam('ilatunits', 'degrees')
        skyctran.setParam('ilngformat', '%12.2h')
        skyctran.setParam('ilatformat', '%12.2h')

    if outFormat == 'radec' or outFormat == 'glonglat':
        skyctran.setParam('olngunits', 'degrees')
        skyctran.setParam('olatunits', 'degrees')
        skyctran.setParam('olngformat', '%11.7f')
        skyctran.setParam('olatformat', '%11.7f')
    elif outFormat == 'hmsdms':
        skyctran.setParam('olngunits', 'hours')
        skyctran.setParam('olatunits', 'degrees')
        skyctran.setParam('olngformat', '%12.2h')
        skyctran.setParam('olatformat', '%12.2h')

    skyctran.setParam('outsystem', outSys)

    coordsList = [str(Xcrd) + " " + str(Ycrd)]
    tmp_transformation = skyctran(Stdin=coordsList, Stdout=1, verbose='no', mode='h')
    new_transformation = tmp_transformation[8].split()

    transformation[0] = new_transformation[2]
    transformation[1] = new_transformation[3]
    return transformation


def natCoords(inImage, RA, DEC, hdrNo=0):
    position = False
    hdr = getheader(inImage, hdrNo)
    tmp_coordinates_1 = str(hdr['CTYPE1']).split("-")
    tmp_coordinates_2 = str(hdr['CTYPE2']).split("-")
    coordinates_1 = tmp_coordinates_1[0]
    coordinates_2 = tmp_coordinates_2[0]
    if coordinates_1 == 'RA' and coordinates_2 == 'DEC':
        position = [float(RA), float(DEC)]
    elif coordinates_1 == 'GLON' and coordinates_2 == 'GLAT':
        position = radec2gal(float(RA), float(DEC))
    return position


def change_frame(c, frame_in, frame_out):
    if frame_in != frame_out:
        c = getattr(c, frame_out)
    return c


def oldcoord2pix(inImage, Xcoord, Ycoord, Zcoord=0):
    from pyraf.iraf import wcsctran
    wcsctran.unlearn()
    wcsctran.setParam('input', 'STDIN')
    wcsctran.setParam('output', 'STDOUT')
    wcsctran.setParam('image', inImage)
    wcsctran.setParam('inwcs', 'world')
    wcsctran.setParam('outwcs', 'logical')
    coordsList = [str(Xcoord) + " " + str(Ycoord) + " " + str(Zcoord) + " " + "1"]
    tmp_transformation = wcsctran(Stdin=coordsList, Stdout=1, verbose='no', mode='h')
    new_transformation = tmp_transformation[0].split()
    return int(round(float(new_transformation[0]))), int(round(float(new_transformation[1])))


def radec2gal(crd1, crd2, frame_in='fk5'):
    crd1 = float(crd1)
    crd2 = float(crd2)
    c = SkyCoord(ra=crd1 * u.degree, dec=crd2 * u.degree, frame=frame_in)
    galactic = c.galactic
    return galactic.l.deg, galactic.b.deg


def oldpix2coord(inImage, Xcoord, Ycoord, hduid=0):
    hdu = fits.open(inImage)
    w = wcs.WCS(hdu[hduid].header)
    ra, dec = w.all_pix2world(int(Xcoord) - 1, int(Ycoord) - 1, ra_dec_order=True)
    return float(ra), float(dec)

# def hmsdms2radec(crd1, crd2, frame_in='fk5', frame_out='fk5'):
#     c = SkyCoord("%s %s" % (crd1, crd2), unit=(u.hourangle, u.deg), frame=frame_in)
#     c = change_frame(c, frame_in, frame_out)
#     ra = c.ra.deg
#     dec = c.dec.deg
#     return ra, dec
#
#
#
#
# def hmsdms2gal(crd1, crd2, frame_in='fk5', frame_out='fk5'):
#     ra, dec = hmsdms2radec(crd1, crd2, frame_in, frame_out)
#     glon, glat = radec2gal(ra, dec, frame_out)
#     return glon, glat
#
#
# def gal2radec(crd1, crd2, frame_in='galactic', frame_out='fk5'):
#     crd1 = float(crd1)
#     crd2 = float(crd2)
#     c = SkyCoord(l=crd1 * u.degree, b=crd2 * u.degree, frame='galactic')
#     radec = getattr(c, frame_out)
#     return radec.ra.deg, radec.dec.deg
#
#
# def gal2hmsdms(crd1, crd2, frame_in='galactic', frame_out='fk5'):
#     crd1 = float(crd1)
#     crd2 = float(crd2)
#     dra, ddec = gal2radec(crd1, crd2, 'galactic', frame_out)
#     ra, dec = radec2hmsdms(dra, ddec, frame_out)
#     return ra, dec
#
#
#
#
#
#
#
#
#
#
#
#
# def angular_distance(Xcoord1, Ycoord1, frame1, Xcoord2, Ycoord2, frame2):
#     c1 = SkyCoord(Xcoord1, Ycoord1, unit=(u.deg, u.deg), frame=frame1)
#     c2 = SkyCoord(Xcoord2, Ycoord2, unit=(u.deg, u.deg), frame=frame2)
#     return c1.separation(c2).arcsecond
#
# # ========== UNUSED ?? ====================
#
# def coord2pixMIR(inImage, Xcoord, Ycoord):
#     from mirpy import miriad
#
#     tmp_reg = "/tmp/" + id_generator() + 'tmpreg.mir'
#     miriad.fits(IN=inImage, out=tmp_reg, op='xyin')
#     coordstr = '%s,%s' % (Xcoord, Ycoord)
#     result = miriad.impos(IN=tmp_reg, coord=coordstr, type='absdeg,absdeg')
#     if os.path.exists(tmp_reg): rmtree(tmp_reg)
#     return result
#
