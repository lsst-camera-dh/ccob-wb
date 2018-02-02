import numpy as np
import calibchannels as chan
import getlist as get
import pyfits

# E2V sensor is (4004 x 4096)
# ITL sensor is (4000 x 4072)
nbx = 4000
nby = 4072

def calibgain(img, gain):
    nx, ny = img.shape
    chany = ny/8
        
    # does it work alos for ITL CCDs ?
    for ic in np.arange(1,16):
        if ic<8:
            img[0:nx/2,ic*chany:(ic+1)*chany]/=gain[ic]
        else:
            img[nx/2:nx,(ic-8)*chany:(ic-7)*chany]/=gain[ic]
    return

pixel_size = 0.010 # mm

# Calibrate pixel gains of CCDs
def calibpixels(iled, iccd):
    hdulist = pyfits.open('beam_CCD%d.fits' %iccd)
    beam = np.array(hdulist[0].data)

    #list of files to calibrate ccd iccd
    # WARNING centered beam must be first or need to develop a routine to find it
    ## it is the one used to calibrate channels gain
    list = get.ccdlist(iccd, iled)
    
    centered = list[0]
    img = GetArray(centered)
    phd = GetPhD(centered)
    img, chan_gain = chan.Smooth_and_Calib(img)

    nx, ny = img.shape()
    pixel_gain = np.zeros((nx,ny),dtype=float)
    npoints = np.zeros((nx,ny),dtype=int)

    phd_ref = -1
    
    for imagefile in list:
        x,y = GetBeamPos(imagefile) # must be in mm
        phd = GetPhD(imagefile)

        # Get position of the centroid of the CCD
        x0, y0 = GetCCDPos(iccd) # must be in mm
        
        img = GetArray(imagefile)

        # correct for channels gain
        calibgain(img, chan_gain)

        # correct for phd intensity
        if phd_ref<0:
            phd_ref=phd
        img/=phd/phd_ref

        # compute beam expected flux vs detected flux
        dx = int([x - x0]/pixel_size + 0.5)
        dy = int([y - y0]/pixel_size + 0.5)
        
        if dx>=0:
            if dy>=0:
                pixel_gain[dx+px:,dy+py:]+=img[dx+px:,dy+py]/beam[:-dx+px,:-dy+py]
                npoints[dx+px:,dy+py:]+=1
            else:
                pixel_gain[dx+px:,:-py+dy]+=img[dx+px:,:-py+dy]/beam[:-dx+px,-dy-py:]
                npoints[dx+px:,:-py+dy]+=1
        else:
            if dy>=0:
                pixel_gain[:-px+dx,dy+py:]+=img[:-px+dx,dy+py]/beam[-dx-px:,:-dy+py]
                npoints[:-px+dx,dy+py:]+=1
            else:
                pixel_gain[:-px+dx,:-py+dy]+=img[:-px+dx,:-py+dy]/beam[-dx-px:,-dy-py:]
                npoints[:-px+dx,:-py+dy]+=1
            

    pixel_gain/=float(npoints)

    hdu = pyfits.PrimaryHDU(pixel_gain)
    hdu.writeto('pixelgain_CCD%d_LED%d.fits' %(iccd,iled))
        

if __name__ == '__main__':
    for iled in range(6):
        for iccd in range(189):
            calibpixels(iled,iccd)



