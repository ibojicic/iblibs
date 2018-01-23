import os

# import sys

os.environ['MPLCONFIGDIR'] = '/tmp'
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.ioff()

import aplpy
# import numpy as np
from matplotlib.path import Path

# import math

__author__ = "ivan"
__date__ = "$May 18, 2012 10:45:26 AM$"


def photFigure(inputFile, outputFile, inRA, inDEC,
               inR_ann=None,
               inR_dann=None,
               inR_app=None,
               PNDiam=None,
               cLevels=None,
               excludeSkyPixels=-1,
               exboxSize=1,
               anColour='yellow',
               apColour='cyan',
               diamColour='blue',
               anLineSt='dashed',
               apLineSt='solid',
               diamLineSt='dashed',
               rLevels=None,
               percmin=2,
               percmax=98,
               valmin=None,
               valmax=None,
               imsize=False,
               crosses=None):
    if cLevels is None:
        cLevels = [0.]
    if rLevels is None:
        rLevels = [0.]

    fig = plt.figure(figsize=(8.67, 7.79))  # ,tight_layout=True)
    fig.subplots_adjust(bottom=0.1, left=0.12, top=0.95, right=0.92)

    markrs = aplpy.FITSFigure(inputFile, figure=fig, convention='calabretta')
    if not imsize or inR_ann is not None:
        if not imsize:
            imsize = inR_dann * 1.3
        markrs.recenter(inRA, inDEC, radius=imsize)
    elif imsize > 0:
        markrs.recenter(inRA, inDEC, radius=imsize)

    markrs.show_grayscale(invert='True', pmin=percmin, pmax=percmax, vmin=valmin, vmax=valmax)
    markrs.add_colorbar()
    if inR_ann is not None:
        markrs.show_circles(inRA, inDEC, inR_ann, edgecolor=anColour, linestyle=anLineSt)
    if inR_dann is not None:
        markrs.show_circles(inRA, inDEC, inR_dann, edgecolor=anColour, linestyle=anLineSt)
    if inR_app is not None:
        markrs.show_circles(inRA, inDEC, inR_app, edgecolor=apColour, linestyle=apLineSt)
    if PNDiam is not None:
        markrs.show_circles(inRA, inDEC, PNDiam / (2. * 3600.), edgecolor=diamColour, linestyle=diamLineSt)
    if cLevels is not None and cLevels.all > 0:
        markrs.show_contour(inputFile, levels=cLevels, colors='#00FF00', smooth=1, convention='calabretta')
    if crosses is not None:
        symSize = 5000.
        crosspath = markCross(250.)
        for X, Y, cross_id in crosses:
            markrs.show_markers(X, Y, marker=crosspath, s=symSize, lw=1., edgecolor='r')
            markrs.add_label(X + 250. / 3600., Y + 250. / 3600., text=cross_id, size=16)

    if rLevels is not None and rLevels.all > 0:
        markrs.show_contour(inputFile, levels=rLevels, colors='red', smooth=1, convention='calabretta',
                            linestyles='dotted')

    if excludeSkyPixels is not None and excludeSkyPixels != -1:
        for excl in excludeSkyPixels:
            worldCoords = markrs.pixel2world(float(excl[0]), float(excl[1]))
            markrs.show_markers(worldCoords[0], worldCoords[1], marker='s', s=exboxSize)

    markrs.save(outputFile)
    plt.close(fig)
    return True


def markCross(size):
    verts = [(-size, 0.), (size, 0.), (0., size), (0., -size)]
    codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
    return Path(verts, codes)

# def makeImage(ImageType,
#               InImages={},
#               outImages={},
#               imSize=False,
#               rgbCube_name=False,
#               tempRGB=False,
#               Xcoord=False,
#               Ycoord=False,
#               majDiam=False,
#               perc_levels={},
#               addGrid=False,
#               imlevels='v',
#               drawBox=[0, 0],
#               boxColour='white',
#               **kwargs):
#
#     Xcoord = float(Xcoord)
#     Ycoord = float(Ycoord)
#
#     fig = plt.figure(figsize=(8.67, 7.79))  # ,tight_layout=True)
#     fig.subplots_adjust(bottom=0.1, left=0.12, top=0.95, right=0.92)
#
#     if ImageType == 'rgb' and (InImages['R'] and InImages['G'] and InImages['B']):
#         # TODO: IF JUST COLOUR LEVELS HAS BEEN CHANGED SKIP MAKING RGB CUBE
#         aplpy.make_rgb_cube([InImages['R'], InImages['G'], InImages['B']], rgbCube_name + ".fits", north=True)
#         if imlevels == 'v':
#             aplpy.make_rgb_image(rgbCube_name + ".fits", tempRGB + ".png",
#                                  vmax_r=float(perc_levels['R']['max']),
#                                  vmax_g=float(perc_levels['G']['max']),
#                                  vmax_b=float(perc_levels['B']['max']),
#                                  vmin_r=float(perc_levels['R']['min']),
#                                  vmin_g=float(perc_levels['G']['min']),
#                                  vmin_b=float(perc_levels['B']['min']),
#                                  embed_avm_tags=False)
#         else:
#             aplpy.make_rgb_image(rgbCube_name + ".fits", tempRGB + ".png",
#                                  pmax_r=float(perc_levels['R']['max']),
#                                  pmax_g=float(perc_levels['G']['max']),
#                                  pmax_b=float(perc_levels['B']['max']),
#                                  pmin_r=float(perc_levels['R']['min']),
#                                  pmin_g=float(perc_levels['G']['min']),
#                                  pmin_b=float(perc_levels['B']['min']),
#                                  embed_avm_tags=False)
#
#         mainCompos = aplpy.FITSFigure(rgbCube_name + "_2d.fits", figure=fig)
#         if imSize != False:
#             mainCompos.recenter(Xcoord, Ycoord, radius=imSize / 2.)
#         mainCompos.show_rgb(tempRGB + ".png")
#
#     elif ImageType == 'intensity' and InImages['in']:
#         mainCompos = aplpy.FITSFigure(InImages['in'], figure=fig, north=True)
#         if imSize != False:
#             mainCompos.recenter(Xcoord, Ycoord, radius=imSize / 2.)
#         mainCompos.show_colorscale(cmap='cubehelix')
#
#         if imlevels == 'v':
#             mainCompos.show_grayscale(vmin=perc_levels['in']['min'], vmax=perc_levels['in']['max'])
#         else:
#             mainCompos.show_grayscale(pmin=perc_levels['in']['min'], pmax=perc_levels['in']['max'])
#     else:
#         return False;
#
#     if 'setaxlabsize' in kwargs:
#         mainCompos.axis_labels.set_font(size=kwargs['setaxlabsize'])
#         mainCompos.tick_labels.set_font(size=kwargs['setaxlabsize'])
#
#     # mainCompos.tick_labels.set_style('colons')
#     mainCompos.tick_labels.set_xformat('hh:mm:ss')
#     mainCompos.tick_labels.set_yformat('dd:mm:ss')
#
#     ax = fig.gca()
#
#     if 'main' in outImages:
#
#         mainCompos.axis_labels.show()
#         mainCompos.tick_labels.show()
#         if addGrid: mainCompos.add_grid()
#         if drawBox != [0, 0]: mainCompos.show_rectangles(Xcoord, Ycoord, drawBox[0] / (3600.), drawBox[1] / (3600.),
#                                                          layer='boxes', lw=1., linestyle='dotted', edgecolor=boxColour)
#         fig.savefig(outImages['main']['filename'], dpi=100)
#
#     if 'thumb' in outImages:
#         mainCompos.axis_labels.hide()
#         mainCompos.tick_labels.hide()
#         mainCompos.ticks.hide()
#         mainCompos.frame.set_linewidth(0)
#         try:
#             mainCompos.remove_layer('boxes')
#         except:
#             sys.exc_clear()
#         try:
#             mainCompos.hide_grid()
#         except:
#             sys.exc_clear()
#         fig.savefig(outImages['thumb']['filename'], dpi=16, bbox_inches='tight', frameon=False)
#
#     if any(overlay in outImages for overlay in ['diameter', 'centroid', 'CS_pos']):
#         mainCompos.axis_labels.show()
#         mainCompos.tick_labels.show()
#         mainCompos.ticks.hide()
#         mainCompos.frame.set_linewidth(0)
#         try:
#             mainCompos.remove_layer('boxes')
#         except:
#             sys.exc_clear()
#         try:
#             mainCompos.remove_grid()
#         except:
#             sys.exc_clear()
#
#
#         if ImageType == 'rgb':
#             mainCompos.hide_colorscale()
#         else:
#             mainCompos.hide_grayscale()
#
#         ax.axes.get_xaxis().set_label_text("RA (J2000)", alpha=0)
#         ax.axes.get_yaxis().set_label_text("DEC (J2000)", alpha=0)
#         labels = [item.get_text() for item in ax.get_xticklabels()]
#         ax.set_xticklabels(labels, alpha=0)
#         labels = [item.get_text() for item in ax.get_yticklabels()]
#         ax.set_yticklabels(labels, alpha=0)
#
#         '''
#         markers :
#             centroid    -> centroid of the nebula: halfcross (corner), pa = extPA (default pa = 0 deg)
#             diameter    -> diameter of the nebula: ellipse, pa = PA (default circle, pa = 0 deg)
#             CS_pos      -> position of the CS: cross, pa = 0deg
#             maxext      -> maximum extension of the nebula: ellipse, pa = extPA (default circle, pa = 0 deg)
#             fchart      -> centroid, diameter and orientation for finding charts
#         '''
#
#         if 'CS_pos' in outImages:
#             if os.path.exists(outImages['CS_pos']['filename']): os.remove(outImages['CS_pos']['filename'])
#             symSize = 5000.
#             crosspath = markCross(250.)
#             mainCompos.show_markers(Xcoord, Ycoord, marker=crosspath, s=symSize, lw=1.,
#                                     edgecolor=outImages['CS_pos']['colour'], layer='CS_pos')
#             fig.savefig(outImages['CS_pos']['filename'], transparent=True, dpi=100)
#             mainCompos.remove_layer('CS_pos')
#
#         if 'diameter' in outImages and float(majDiam) > 0. and float(majDiam) > 0.:
#             if os.path.exists(outImages['diameter']['filename']): os.remove(outImages['diameter']['filename'])
#             mainCompos.show_ellipses(Xcoord, Ycoord, float(majDiam) / 3600.,
#                                      float(majDiam) / 3600., float(0.), lw=1.,
#                                      edgecolor=outImages['diameter']['colour'], layer='diameter',
#                                      linestyle=outImages['diameter']['linetype'])
#             fig.savefig(outImages['diameter']['filename'], transparent=True, dpi=100)
#             mainCompos.remove_layer('diameter')
#
#         if 'centroid' in outImages:
#             if os.path.exists(outImages['centroid']['filename']): os.remove(outImages['centroid']['filename'])
#             symSize = 5000.
#             crosspath = markCorn(250.)
#             mainCompos.show_markers(Xcoord, Ycoord, marker=crosspath, s=symSize, lw=1.,
#                                     edgecolor=outImages['centroid']['colour'], layer='centroid')
#             fig.savefig(outImages['centroid']['filename'], transparent=True, dpi=100)
#             mainCompos.remove_layer('centroid')
#
#     if 'fchart' in outImages:
#         if os.path.exists(outImages['fchart']['filename']): os.remove(outImages['fchart']['filename'])
#
#         mainCompos.axis_labels.show()
#         mainCompos.tick_labels.show()
#
#         symSize = 5000.
#         crosspath = markCross(250.)
#         mainCompos.show_markers(Xcoord, Ycoord, marker=crosspath, s=symSize, lw=1.,
#                                 edgecolor='r', layer='fchart')
#
#         if majDiam > 20:
#             mainCompos.show_ellipses(Xcoord, Ycoord, float(majDiam) / 3600.,
#                                      float(majDiam) / 3600., float(0.), lw=1.,
#                                      edgecolor='r')
#
#         xoffset = imSize / 2. - imSize / 30.
#         yoffset = imSize / 2. - imSize / 10.
#         xarrow = Xcoord - xoffset / math.cos(math.radians(Ycoord))
#         if Ycoord < 0: yoffset = -1 * yoffset
#         yarrow = Ycoord - yoffset
#         # centarrow = mainCompos.pixel2world(xarrow,yarrow)
#         mainCompos.show_arrows(xarrow, yarrow, 0, imSize / 16., color='g')
#         mainCompos.show_arrows(xarrow, yarrow, imSize / (math.cos(math.radians(Ycoord)) * 16.), 0, color='g')
#         mainCompos.add_label(xarrow, yarrow + imSize / 13., "N", color='g')
#         mainCompos.add_label(xarrow + imSize / (math.cos(math.radians(Ycoord)) * 13.), yarrow, "E", color='g')
#
#         fig.savefig(outImages['fchart']['filename'], transparent=True, dpi=100)
#         mainCompos.remove_layer('fchart')
#
#     plt.close()
#
#     return True
#
#
# def removeTempFiles(tempRGB, out_format_main, out_format_thumb):
#     if os.path.exists("%s.%s" % (tempRGB, out_format_main)): os.remove("%s.%s" % (tempRGB, out_format_main))
#     if os.path.exists("%s_thumb.%s" % (tempRGB, out_format_thumb)): os.remove(
#         "%s_thumb.%s" % (tempRGB, out_format_thumb))
#
#
# def makeFindChart(IntImage=False, rImage=False, gImage=False, bImage=False, \
#                   imSize=False, \
#                   out_name=False, rgbCube_name=False, tempRGB=False, \
#                   Xcoord=False, Ycoord=False, majDiam=False, \
#                   r_imSc_max=99., g_imSc_max=99., b_imSc_max=99., \
#                   r_imSc_min=5., g_imSc_min=5., b_imSc_min=5., \
#                   imSc_min=5., imSc_max=99., \
#                   outformatmain='png', \
#                   addGrid=True, \
#                   redoRGB='n', \
#                   imlevels='v', \
#                   addlabel=False, \
#                   fchartmpos='y', fchartmdiam='y', fchartmorien='y', \
#                   **kwargs):
#     if (os.path.exists(tempRGB + "." + outformatmain)): os.remove(tempRGB + "." + outformatmain)
#
#     Xcoord = float(Xcoord)
#     Ycoord = float(Ycoord)
#
#     if not IntImage and (rImage and gImage and bImage):
#         making = 'rgb'
#         compos = [rImage, gImage, bImage]
#         if (not (os.path.exists(rgbCube_name + ".fits")) or not (
#                 os.path.exists(rgbCube_name + "_2d.fits")) or redoRGB == 'y'):
#             aplpy.make_rgb_cube(compos, rgbCube_name + ".fits", north=True)
#             # HERE
#         if imlevels == 'v':
#             aplpy.make_rgb_image(rgbCube_name + ".fits", tempRGB + "." + outformatmain, vmax_r=float(r_imSc_max),
#                                  vmax_g=float(g_imSc_max), vmax_b=float(b_imSc_max), vmin_r=float(r_imSc_min),
#                                  vmin_g=float(g_imSc_min), vmin_b=float(b_imSc_min), embed_avm_tags=False)
#         else:
#             aplpy.make_rgb_image(rgbCube_name + ".fits", tempRGB + "." + outformatmain, pmax_r=float(r_imSc_max),
#                                  pmax_g=float(g_imSc_max), pmax_b=float(b_imSc_max), pmin_r=float(r_imSc_min),
#                                  pmin_g=float(g_imSc_min), pmin_b=float(b_imSc_min), embed_avm_tags=False)
#     else:
#         making = 'int'
#
#     fig = plt.figure(figsize=(8.67, 7.79), tight_layout=True)
#
#     # fig = plt.figure(figsize=(3.3,3),tight_layout=True)
#
#     if making == 'rgb':
#         mainCompos = aplpy.FITSFigure(rgbCube_name + "_2d.fits", figure=fig)
#         mainCompos.recenter(Xcoord, Ycoord, radius=imSize / 2.)
#         mainCompos.show_rgb(tempRGB + "." + outformatmain)
#     else:
#         mainCompos = aplpy.FITSFigure(IntImage, north=True, figure=fig)
#         mainCompos.recenter(Xcoord, Ycoord, radius=imSize / 2.)
#         if imlevels == 'v':
#             mainCompos.show_grayscale(vmin=imSc_min, vmax=imSc_max, invert=True)
#         else:
#             mainCompos.show_grayscale(pmin=imSc_min, pmax=imSc_max, invert=True)
#
#     if 'setaxlabsize' in kwargs:
#         mainCompos.axis_labels.set_font(size=kwargs['setaxlabsize'])
#         mainCompos.tick_labels.set_font(size=kwargs['setaxlabsize'])
#
#     # mainCompos.tick_labels.set_style('colons')
#     mainCompos.tick_labels.set_xformat('hh:mm:ss')
#     mainCompos.tick_labels.set_yformat('dd:mm:ss')
#     mainCompos.axis_labels.show()
#     mainCompos.tick_labels.show()
#
#     if addlabel: mainCompos.add_label(0.0, -0.06, addlabel, relative=True, color='red')
#     if addGrid: mainCompos.add_grid()
#
#     '''
#         markers :
#             centroid    -> centroid of the nebula: cross, pa = 0deg
#             diameter    -> diameter of the nebula: ellipse, pa = PA (default circle, pa = 0 deg)
#     '''
#     if fchartmpos == 'y':
#         symSize = 5000.
#         crosspath = markCross(250.)
#         mainCompos.show_markers(Xcoord, Ycoord, marker=crosspath, s=symSize, lw=1.)
#
#     if fchartmdiam == 'y' and majDiam > 20:
#         mainCompos.show_circles(Xcoord, Ycoord, majDiam / (2. * 3600.), lw=1.)
#
#     if fchartmorien == 'y':
#         # xarrow = 15.4 * fig.dpi * (imSize * 3600.) / 300.
#         # yarrow = 15.0 * fig.dpi * (imSize * 3600.) / 300.
#         xoffset = imSize / 2. - imSize / 30.
#         yoffset = imSize / 2. - imSize / 10.
#
#         xarrow = Xcoord - xoffset / math.cos(math.radians(Ycoord))
#         if Ycoord < 0: yoffset = -1 * yoffset
#         yarrow = Ycoord - yoffset
#
#         # centarrow = mainCompos.pixel2world(xarrow,yarrow)
#         mainCompos.show_arrows(xarrow, yarrow, 0, imSize / 16., color='g')
#         mainCompos.show_arrows(xarrow, yarrow, imSize / (math.cos(math.radians(Ycoord)) * 16.), 0, color='g')
#         mainCompos.add_label(xarrow, yarrow + imSize / 13., "N", color='g')
#         mainCompos.add_label(xarrow + imSize / (math.cos(math.radians(Ycoord)) * 13.), yarrow, "E", color='g')
#
#     fig.savefig(out_name + '_fchart.png', dpi=100)
#     plt.close()
#
#     if (os.path.exists(tempRGB + "." + outformatmain)): os.remove(tempRGB + "." + outformatmain)
#
#     return
#
#
# def pointsAndFit(x_obs, y_obs, x_err, y_err, x_fit, y_fit, pltname, x_label, y_label, addtext=None,
#                  axlimits=[None, None, None, None]):
#     fig = plt.figure(figsize=(5, 5))
#
#     ax = fig.add_subplot(1, 1, 1)
#
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#
#     ax.set_xlim(axlimits[0], axlimits[1])
#     # ax.set_ylim(axlimits[2],axlimits[3])
#
#     # fig.gca().set_xlabel(r'$\nu$')
#     ax.set_xlabel(x_label)
#     ax.set_ylabel(y_label)
#
#     ax.errorbar(x_obs, y_obs, fmt='b.', yerr=y_err, color='black')
#
#     ax.plot(x_fit, y_fit, color='black', ls='dashed')
#
#     axx = ax.axis()
#
#     if addtext is not None:
#         ax.text(0.7, axx[3] - 6 * (axx[3] - axx[2]) / 10, addtext)
#
#     plt.savefig(pltname)
#
#
# def get_histo_no_bins(N):
#     # sigmaplot
#     return int(3. + math.log10(float(N)) * math.log10(float(N)) / math.log10(2.))
#
#
# def plotHisto(input_list, h_min, h_max, x_label='X', y_label='Y'):
#     # h_min,h_max = None,None
#     # for h_list in input_list:
#     #     if h_min is None or h_min > min(h_list['data']):
#     #         h_min = min(h_list['data'])
#     #     if h_max is None or h_max < max(h_list['data']):
#     #         h_max = max(h_list['data'])
#
#     max_bins = None
#     for h_list in input_list:
#         hist, edges = np.histogram(h_list['data'], bins='auto', range=(h_min, h_max))
#         if max_bins is None or len(max_bins) > len(edges):
#             max_bins = edges
#
#     h_params = dict(
#         bins=max_bins,
#         range=(h_min, h_max),
#         alpha=1.  # 0.3,
#         # normed = True
#     )
#
#     # fig = plt.figure()
#     # ax = fig.add_subplot(111)
#     # NoBins = get_histo_no_bins(len(inputList['data']))
#     for h_list in input_list:
#         n, bins, patches = ax.hist(h_list['data'], histtype=h_list['htype'], **h_params)
#
#     # ax.set_xlabel(x_label)
#     # ax.set_ylabel(y_label)
#     #
#     # xminorLocator = MultipleLocator(1)
#     # yminorLocator = MultipleLocator(1)
#     #
#     # ax.xaxis.set_minor_locator(xminorLocator)
#     # ax.yaxis.set_minor_locator(yminorLocator)
#
#
#     return plt
#
#
#
#
#
#
# def markCorn(size):
#     verts = [(-size, 0.), (0., 0.), (0., -size)]
#     codes = [Path.MOVETO, Path.LINETO, Path.LINETO]
#     return Path(verts, codes)
#
#
# def overlayImage(inRA, inDEC, intImage, contImage, contours1, resImage, inR_app):
#     image = aplpy.FITSFigure(intImage)
#     image.show_grayscale()
#     image.show_contour(data=contImage, levels=contours1, colors='yellow', smooth=1)
#     image.show_circles(inRA, inDEC, inR_app, edgecolor='cyan')
#     image.save(resImage + ".png")
#     return
#
#
#
#
# # Set up axes and plot some awesome science
# class XYData:
#     def __init__(self, sqlres, colourin, colourout, char, size, leglabel, plottype='scatter', annotate='no'):
#
#         self.x, self.y, self.fora = self.xydata(sqlres)
#         self.colourin = colourin
#         self.colourout = colourout
#         self.char = char
#         self.size = size
#         self.leglabel = leglabel
#         self.plottype = plottype
#         self.annotate = annotate
#
#     def xydata(self, sqlres):
#         x = []
#         y = []
#         fora = []
#         for tmpdata in sqlres:
#
#             try:
#                 x1 = tmpdata[0]
#             except:
#                 print 'errrrrrr'
#                 sys.exit()
#             x.append(x1)
#
#             try:
#                 y1 = tmpdata[1]
#             except:
#                 print 'errrrrrr'
#                 sys.exit()
#             y.append(y1)
#
#             try:
#                 id = tmpdata[2]
#             except:
#                 id = 'dummy'
#
#             try:
#                 xoff = tmpdata[3]
#             except:
#                 xoff = 0.02
#
#             try:
#                 yoff = tmpdata[4]
#             except:
#                 yoff = 0.02
#
#             fora.append((x1, y1, id, xoff, yoff))
#
#         return x, y, fora
