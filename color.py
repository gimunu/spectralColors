#!/usr/bin/env python

# Copyright (C) 2013 Umberto De Giovannini 
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


import numpy as np
from scipy.interpolate import interp1d #interpolation 
from scipy import integrate
try:
    import pylab as pl
except ImportError:
    pl = None
try:    
    from colormath.color_objects import XYZColor
except ImportError:
    XYZColor = None

#standard modules for python > 2.7
import math
import argparse
import os.path
import sys
import logging

###########################################################
#
# Convert a tristimulus color value XYZ to target 
# color-space
#
###########################################################

def XYZ_to(X,Y,Z,outc, which = None, si= None):
    "Converts an XYZ color to a target color space"

    if XYZColor: color = XYZColor(X, Y, Z)
        
    if outc == 'xyY':    
        out = [X/(X+Y+Z),Y/(X+Y+Z),Y]
    elif 'RGB' in outc   :   
        if XYZColor and which not in [None, 'srgb-Nosi']:
            #THIS part gives fishy results since there is a whitepoint transformation 
            # due to an assumed default illuminant (d65) no matter what
            color = XYZColor(X, Y, Z,)
            color.illuminant = si
            print "illuminant",color.illuminant
            print "color",color
            rgb = color.convert_to('rgb', target_rgb=which)
            out = [rgb.rgb_r/255., rgb.rgb_g/255., rgb.rgb_b/255.]
        else:
            #See http://en.wikipedia.org/wiki/SRGB
            M = np.matrix('3.2404542, -1.5371385, -0.4985314;\
                          -0.9692660,  1.8760108,  0.0415560;\
                           0.0556434, -0.2040259,  1.0572252')
            # M = np.matrix('3.2406, -1.5372, -0.4986;\
            #               -0.9689,  1.8758,  0.0415;\
            #                0.0557, -0.2040,  1.0570')
        
                           
            XYZ = np.matrix('%f %f %f' % (X,Y,Z))               
            RGB = M * (XYZ.T) 
            # print M
            # print RGB
            if RGB[0] > 0.0031308:
                if RGB[0] >= 1.0 :
                    RGB[0] = 1.0
                else:
                    RGB[0] = (1.055 * math.pow(RGB[0], 1.0 / 2.4)) - 0.055
            elif RGB[0] >= 0.0:
                RGB[0] = RGB[0] * 12.92
            else:
                RGB[0] = 0.0

            if RGB[1] > 0.0031308:
                if RGB[1] >= 1.0 :
                    RGB[1] = 1.0
                else:
                    RGB[1] = (1.055 * math.pow(RGB[1], 1.0 / 2.4)) - 0.055
            elif RGB[1] >= 0.0:
                RGB[1] = RGB[1] * 12.92
            else:
                RGB[1] = 0.0
   
            if RGB[2] > 0.0031308:
                if RGB[2] >= 1.0 :
                    RGB[2] = 1.0
                else:
                    RGB[2] = (1.055 * math.pow(RGB[2], 1.0 / 2.4)) - 0.055
            elif RGB[2] >= 0.0:
                RGB[2] = RGB[2] * 12.92
            else:
                RGB[2] = 0.0

            out = RGB
            

    elif outc == 'XYZ'  :   
        out = [X,Y,Z]
    else:
        out = NUL
        
            
    return out    




######################## MAIN ################################


VERSION = '0.1 alpha'

RESOURCES = './Resources'


desc="""
This utility calculates the color associated with a spectral distribution under 
several different illumination conditions.
Colors from reflection, transmission and direct light spectral distributions can 
be calculated.
"""

parser = argparse.ArgumentParser(version='%s version %s' %(sys.argv[0],VERSION),
                                 description=desc)

parser.add_argument('-q', '--quiet', action='store_true',
                    help='Print out the color in the chosen color-space only. This option is for scripting purposes.')

parser.add_argument('-g', '--graphic', action='store_true',
                    help='Plot some fancy graphic output.')    

parser.add_argument('-a', '--absorb', action='store_true',
                    help='The input file is to be considered as an absorbance spectrum (i.e reflectance or transmittance). \
                    In order to obtain the color an illuminant source has to be specified [-s]. \
                    By default --absorb=%(default)s.')                      

cspaces= ['XYZ', 'xyY', 'RGB']
csp_specification= ['RGB'] #Which color space accepts additional specifications
class ValidateCspace(argparse.Action):
    CHOICES = cspaces
    def __call__(self, parser, namespace, values, option_string=None):
        specify = False
        spaces = []
        idx = -1
        if values:
            for value in values:
                # print 'idx',idx, value, specify
                if specify and value not in self.CHOICES:
                    spaces[idx][1] = value
                    specify = False
                    continue
                else:
                    spaces.append([value, None]) # Default rgb is the one calculated without colormath
                    idx+=1

                if value in csp_specification:
                    specify = True
                
                    
                if value not in self.CHOICES and not specify:
                    message = ("invalid choice: {0!r} (choose from {1})"
                               .format(value,
                                       ', '.join([repr(action)
                                                  for action in self.CHOICES])))

                    raise argparse.ArgumentError(self, message)

                    
            setattr(namespace, self.dest, spaces)
                    
parser.add_argument('-c', '--cspace', action=ValidateCspace,  nargs ='*', default=[['XYZ', None]] ,
                    help='Output color space. Options (default value [%(default)s] ): XYZ, xyY, RGB.\
                    The kind RGB color space can be specified as argument:\'-c RGB apple_rgb\'.\
                    If specified, the RGB values will be calculated with the Python-Colormath external \
                    library (http://code.google.com/p/python-colormath/). If no specification is provided the RGB values\
                    are in the standard sRGB colorspace.')

silist =['d50','d55','d65','d75'] 

parser.add_argument('-s', '--source', action='store',  default='d65', choices = silist,
                    help='illuminant light source type. Options (default value [%(default)s]): d{50,55,65,75} family. ')

parser.add_argument('-o', '--observer', 
                    help='Standard observer color matching functions. Options (default value [%(default)s]):\
                          ciexyz31, ciexyzj, ciexyzjv, ciexyz64, lin2012xyz2e_5_7sf',
     action='store', default='ciexyz64')

parser.add_argument('-f','--file', nargs='+', type=argparse.FileType('r'), default=[sys.stdin],
                    help='The input spectral power distribution I(lambda) must be in nm and defined at least in the\
                          (visible) range lambda=[380, 780] nm. Multiple input files or direct stdin stream are accepted.')


args = parser.parse_args()

# print  "colorspace",args.cspace

#Define logged print 'lprint' in order to control output
logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                    format="%(message)s")
lprint = logging.info


#####################
# Options error check
#####################
# Illuminants 

# if args.source not in silist:
#     # raise ValueError('invalid color matching function'.format(s=args.cmf))
#     print args.source," not valid illuminant source. Possible options are: \n",silist
#     exit(1)

# Color matching functions 
cmflist=['ciexyz31', 'ciexyzj', 'ciexyzjv', 'ciexyz64','lin2012xyz2e_5_7sf']
if args.observer not in cmflist:
    # raise ValueError('invalid color matching function'.format(s=args.cmf))
    print args.observer," not valid standard observer cmf. Possible options are: \n",cmflist
    exit(1)
#####################

# CMFs
lprint ('Color matching functions: %s' % args.observer)

CMF = np.genfromtxt('%s/%s.csv' %(RESOURCES, args.observer), delimiter=',')

xbar = interp1d(CMF[:,0], CMF[:,1], kind='linear')
ybar = interp1d(CMF[:,0], CMF[:,2], kind='linear')
zbar = interp1d(CMF[:,0], CMF[:,3], kind='linear')

# The grid of wavelenghts we are going to use 
# is bounded by cmfs range
# xrange = np.linspace(np.amin(CMF[:,0]), np.amax(CMF[:,0]), 1000)
xrange = np.linspace(380, 780, 800)


#Illuminant
if args.absorb:
    Sdata = np.genfromtxt('%s/si%s.csv' %(RESOURCES, args.source), delimiter=',')
    S = interp1d(Sdata[:,0], Sdata[:,1], kind='linear')

    Xs = integrate.trapz(S(xrange) * xbar(xrange), x= xrange, dx = xrange[1]-xrange[0])
    Ys = integrate.trapz(S(xrange) * ybar(xrange), x= xrange, dx = xrange[1]-xrange[0])
    Zs = integrate.trapz(S(xrange) * zbar(xrange), x= xrange, dx = xrange[1]-xrange[0])
    lprint ('Source illuminant : %s (XYZ=[%1.4f, %1.4f, %1.4f])' % (args.source, Xs/Ys,Ys/Ys,Zs/Ys))
    


if args.graphic:
    pl.rc('font', family='serif')
    pl.rc('font', size=12)
    # pl.rc('axes', labelsize=14)    
    pl.rc('legend', fontsize=10)    
    
    fig = pl.figure()
    
    p1 = fig.add_subplot(2, 2, 1, title = 'Color matching functions \'%s\'' % args.observer)
    
    #CMF panel
    p1.plot(CMF[:,0], CMF[:,1], lw=2, color='r', label='xbar')
    p1.plot(CMF[:,0], CMF[:,2], lw=2, color='g', label='ybar')
    p1.plot(CMF[:,0], CMF[:,3], lw=2, color='b', label='zbar')
    p1.set_xlabel('$\lambda$ [nm]') 
    # p1.ylabel('$Y$') 
    p1.legend()
    
    #Spectra
    p2 = fig.add_subplot(2, 2, 2, title='Spectra')
    p2.set_xlabel('$\lambda$  [nm]') 
    if args.absorb: p2.plot(xrange, S(xrange)/np.amax(S(xrange)), c='r', lw=2, label='Source %s/%1.1e'%(args.source,np.amax(S(xrange))))

    #Plot the point on the chromaticity diagram 
    p3 = fig.add_subplot(2, 2, 3,  title='Chromaticity diagram',aspect = 1.0)
    p3.set_xlabel('$x$') 
    p3.set_ylabel('$y$') 
    #diagram edges
    Cx = xbar(xrange) / (xbar(xrange) + ybar(xrange) + zbar(xrange))
    Cy = ybar(xrange) / (xbar(xrange) + ybar(xrange) + zbar(xrange))
    p3.plot(Cx, Cy, lw=1, color ='k')
    p3.plot([Cx[0],Cx[-1]], [Cy[0],Cy[-1]], lw=1, color='k')
    if args.absorb: p3.plot(Xs/(Xs+Ys+Zs),Ys/(Xs+Ys+Zs), 'r+',label=file)    
    
    #Color palette
    display_color = None # By default do not use colormath to calculate rgb color to display 
    for space in args.cspace:
        if space[0] == 'RGB':
            display_color=space[1]
    p4 = fig.add_subplot(2, 2, 4, title='Color(s) - RGB (%s)' %(display_color), aspect = 1.0)
    p4.set_ylim(0,1)
    p4.set_xlim(0,1)
    p4.set_xticks([ ])
    p4.set_yticks([ ])
    #find the nicer arrangement
    delta=0.00 #empty space between elements
    ncolors=len(args.file)
    nrows = round(np.sqrt(ncolors))
    ncols = np.ceil(ncolors/nrows)
    w = 1./ncols
    h = 1./nrows
    ic = 0
    ir = 0



#Loop though the file list from command
for file in args.file:
    
        
    lprint ("Color for: %s" % (file.name))

    # The Input spectral density 
    Idata = np.genfromtxt(file)
    I = interp1d(Idata[:,0], Idata[:,1], kind='linear')
    Itab = I(xrange)
    Itab /= np.amax(Itab) # rescale to have values in [0,1] - numerically stabler   
    
    if args.absorb:

        #XYZ for the absorbed light S*I 
        X = integrate.trapz(Itab * S(xrange) * xbar(xrange), x= xrange, dx = xrange[1]-xrange[0])
        Y = integrate.trapz(Itab * S(xrange) * ybar(xrange), x= xrange, dx = xrange[1]-xrange[0])
        Z = integrate.trapz(Itab * S(xrange) * zbar(xrange), x= xrange, dx = xrange[1]-xrange[0])
        #We have to normalize all the values by the photopic response, 
        #i.e. the total energy of the illuminant seen by the cones
        X /=Ys
        Z /=Ys
        Y /=Ys
        
    else:
        X = integrate.trapz(Itab * xbar(xrange), x= xrange, dx = xrange[1]-xrange[0])
        Y = integrate.trapz(Itab * ybar(xrange), x= xrange, dx = xrange[1]-xrange[0])
        Z = integrate.trapz(Itab * zbar(xrange), x= xrange, dx = xrange[1]-xrange[0])
        #If we are calculating color of direct light we need to normalize 
        #by the total energy on the cones Y 
        X /=Y
        Z /=Y
        Y /=Y

    #######################################
    ############ TEXT OUTPUT ##############
    #######################################
    for space in args.cspace:
        v = XYZ_to(X,Y,Z,space[0],which=space[1], si=args.source)
        if args.quiet:
            print '%1.4f  %1.4f  %1.4f  ' % (v[0],v[1],v[2]),
        else:    
            if space[1]:
                spname= '%s(%s)' %(space[0],space[1])
            else: 
                spname= '%s' %(space[0])
            print '%s = [%1.4f, %1.4f, %1.4f]' % (spname, v[0],v[1],v[2])

    if args.quiet: print'' #newline        

    #######################################
    ########## GRAPHIC OUTPUT #############
    #######################################
    if args.graphic:
        #Plot the point on the chromaticity diagram 
        xyY=XYZ_to(X,Y,Z,'xyY')
        rgb=XYZ_to(X,Y,Z,'RGB', which=display_color, si = args.source)
        p3.plot(xyY[0], xyY[1], c=(rgb[0],rgb[1],rgb[2]), marker='o', label=file.name)
        
        #Plot the (normalized) input spectra
        p2.plot(xrange, Itab, c=(rgb[0],rgb[1],rgb[2]), lw=2, label=file.name)
        p2.legend()
        
        #The colors
        ix = ic * w + delta
        iy = (1 - (ir + 1) * h) + delta
        p4.broken_barh([(ix , w - 2 * delta)], 
                        (iy , h - 2 * delta), 
                        facecolors=(rgb[0],rgb[1],rgb[2]))
        ic +=1
        if ic == ncols: #newline
            ic = 0
            ir += 1
   
    
if args.graphic: 
    pl.tight_layout() #No overlapping subplots  
    pl.show()   
    

lprint ("Done")


