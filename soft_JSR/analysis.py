import numpy as np
import calibchannels as chan
import getlist as get

# list of images used for beam calib
# beam centered on each CCD
list1 = get.beamlist()

# E2V sensor is (4004 x 4096)
# ITL sensor is (4000 x 4072)
nbx = 4000
nby = 4072
beam = np.zeros(nbx,nby)
nbeam = np.zeros(nbx,nby)
ccd_gain_ref=-1;

# Compute Beam Array
for imagefile in list1:

    # We suppose that the beam is centered on the CCD so x=x0 and y=y0
    # If we need to compute the beam (more) precise position relatively to the CCD
    # We'll have to correct dx and dy
    #x,y = GetBeamPosition(imagefile)
    #x0, y0 = GetCCDPosition(iccd)
    # Note that dx and dy are given in pixel unit : dx = int([x (mm) - x0 (mm)]/pixel_size(mm) + 0.5)
    dx=0
    dy=0
    phd = GetPhD(imagefile)
    iccd = GetCCD(imagefile)
    img = GetArray(imagefile)

    nx, ny = img.shape
    # in case of E2V sensor we don't use edge pixels in order to have the same are for the reconstructed beam
    px = (nx-nbx)/2
    py = (ny-nby)/2
    
    img, chan_gain = chan.Smooth_and_Calib(img)

    # warning amax works only if no CR saturating the CCD (must have been removed previously)
    ccd_gain = np.amax(img)/phd
    if ccd_gain_ref<0:
        ccd_gain_ref=ccd_gain
    ccd_gain/=ccd_gain_ref
    print 'chan_gain',chan_gain, 'ccd_gain',ccd_gain

    if dx>=0:
        if dy>=0:
            beam[:-dx,:-dy]+=img[dx+px:-px,dy+py:-py]/ccd_gain
            nbeam[:-dx,:-dy]+=1
        else:
            beam[:-dx,dy:]+=img[dx+px:-px,py:-py-dy]/ccd_gain
            nbeam[:-dx,dy:]+=1
    else:
        if dy>=0:
            beam[dx:,:-dy]+=img[px:-px-dx,dy+py:-py]/ccd_gain
            nbeam[dx:,:-dy]+=1
        else:
            beam[dx:,dy:]+=img[px:-px-dx,py:-py-dy]/ccd_gain
            nbeam[dx:,dy:]+=1
    
    
# Calibrate pixel gains of CCDs
for iccd in 0, 189:
    #list of files to calibrate ccd iccd
    # WARNING centered beam must be first
    list2 = get.ccdlist(iccd)
    centered = list2[0]
    img = GetArray(centered)
    img, chan_gain = chan.Smooth_and_Calib(img)
    for imagefile in list1:
        x,y = GetPos(imagefile)
        phd = GetPhD(imagefile)
        iccd = GetCCD(imagefile)
        x0, y0 = GetCCDPos(iccd)
        img = GetArray(imagefile)
        img, g = chan.Smooth_and_Calib(img)
        




