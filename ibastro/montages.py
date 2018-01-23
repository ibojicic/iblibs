
import os
import shutil
import traceback

import montage_wrapper as montage
from qb_common.strings import id_generator


def download_postage_stamps(survey_name, band_name, DRAJ2000, DDECJ2000, Xwidth, outfile, workspace=False,
                            delete_workspace=True):
    coords = '%s %s' % (DRAJ2000, DDECJ2000)
    tmp_region_header = '/tmp/' + id_generator() + ".hdr"
    if not workspace:
        workspace = '/tmp/MOSAIC_' + id_generator()
    montage.mHdr(coords, width=Xwidth, out_file=tmp_region_header, system='eq')
    montage.mExec(survey=survey_name, band=band_name, region_header=tmp_region_header, output_image=outfile,
                  workspace_dir=workspace)
    if os.path.exists(tmp_region_header): os.remove(tmp_region_header)
    if delete_workspace and os.path.exists(workspace): shutil.rmtree(workspace)


def FITScutout(inImage, outImage, Xcoord, Ycoord, box_x_coord, box_y_coord  = False):
    Xcoord = float(Xcoord)
    Ycoord = float(Ycoord)

    box_x_coord = float(box_x_coord)

    if not box_y_coord:
        box_y_coord = box_x_coord
    box_y_coord = float(box_y_coord)
    print inImage, outImage, Xcoord, Ycoord, box_x_coord / 3600., box_y_coord / 3600.
    try:
        res = montage.mSubimage(inImage, outImage, Xcoord, Ycoord, xsize=box_x_coord / 3600., ysize=box_y_coord / 3600.)
        results = parse_montage_out(str(res))
        if results['stat'] != 'OK' or results['content'] != 'normal':
            if os.path.exists(outImage):
                os.remove(outImage)
            return False
    except Exception:
        traceback.print_exc()  # print full exception
        return False
    return True


def rebinImage(inImage, outImage, factor):
    montage.mShrink(inImage, outImage, factor)


def pix2coord(inImage, Xcoord, Ycoord):
    result = montage.mPix2Coord(inImage, Xcoord, Ycoord)
    return result.lon, result.lat


def parse_montage_out(mon_out):
    result = {}
    tmp1 = mon_out.split("\n")
    for line in tmp1:
        tmp2 = line.split(":")
        result[tmp2[0].strip()] = tmp2[1].strip()
    return result
