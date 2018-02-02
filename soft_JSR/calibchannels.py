import numpy as np
import scipy.ndimage as ndimage
from scipy import interpolate
from scipy import optimize

ccobfile = "/sps/lsst/DataBE/lpnws5203/data/frames//20170616/10_CCD1_20170616140354.fz"

bx = 40
by = 5
edge = 20

# several functions used to get ccd data without border areas #

def safe(nx,ny):
    mx = np.append(np.arange(edge,nx/2-bx),np.arange(nx/2+bx,nx-edge))
    my = np.array([])
    vecy = np.arange(ny)
    chany = ny/8
    for i in np.arange(8):
        for val in vecy[chany*i+by:chany*(i+1)-by]:
          my = np.append(my,int(val))
    my=my[edge-by:-(edge-by)].astype(int)
    return mx,my

lfit = 150

def getvxvy(ic, nx, ny):
     return getvx(ic,nx),getvy(ic,ny)

def getvx(ic, nx):
     if (ic>7):
          return np.arange(edge,nx/2-bx)
     else:
          return np.arange(nx/2+bx,nx-edge)

def getvy(ic, ny):
     chany = ny/8
     if ic>7:
          ic-=8
     zr = by
     zl = by
     if (ic==0):
          zl=edge
     if (ic==7):
          zr=edge
     vy = np.arange(chany*ic+zl,chany*(ic+1)-zr)
     return vy

def getsafevy(ic, ny):
    chany = ny/8
    if ic>7:
        ic-=8
    vy = np.append(np.arange(chany*ic-lfit,chany*ic-by), np.arange(chany*ic+by,chany*ic+lfit))
    return vy

def getbordervy(ic, ny):
    chany = ny/8
    if ic>7:
        ic-=8
    vy = np.arange(chany*ic-by, chany*ic+by)
    return vy

def poly(p, x):
    nx = x.size
    g = np.ones(nx)
    g[nx/2:]=p[3]
    return g*(p[0] * x*x + p[1] * x + p[2])

def errfuncy(p, x, y, err):
    return (y - poly(p,x)) / err

# calib ccd_data with channels gain
def calibgain(img, gain):
    nx, ny = img.shape
    chany = ny/8
    array = np.copy(img)
    for ic in np.arange(1,16):
        if ic<8:
            array[0:nx/2,ic*chany:(ic+1)*chany]=img[0:nx/2,ic*chany:(ic+1)*chany]/gain[ic]
        else:
            array[nx/2:nx,(ic-8)*chany:(ic-7)*chany]=img[nx/2:nx,(ic-8)*chany:(ic-7)*chany]/gain[ic]
    return array

# compute channels gain
def CalibGainPerChannel(img):
    array = ndimage.gaussian_filter(img, sigma=(9, 9), order=0)

    nx, ny = array.shape
    safe_x, safe_y = safe(nx,ny)
    
    chany = ny/8

    gain = np.ones(16)
    
    step = 100

    # compute gain of channel ic with respect to channel ic-1
    # exept for channel 0 and 8
    # gain[1:7] relative to channel 0
    # gain[9:15] relative to channel 8
    # NOTE that in this code channels are ordered
    # 0  1  2  3  4  5  6  7
    # 8  9 10 11 12 13 14 15
    # so my numbers do not correspond to the real channel numbers
    
    for ic in np.arange(1,16):
        if (ic==8):
            continue
        xs = 0
        xe = nx/2

        if ic>7:
            xs=nx/2
            xe=nx

        vx = safe_x[np.logical_and(safe_x>=xs,safe_x<xe)][::step]
                    
        vy = getsafevy(ic,ny)
        
        hg = np.array([])

        pinit = [1, 1, 1, 1.0]
        for ix in vx:
            vay = array[ix,vy]
                        
            vae = np.sqrt(vay)
            out = optimize.leastsq(errfuncy, pinit,
                                   args=(vy, vay, vae), full_output=1)
            pout = out[0]
            hg = np.append(hg,pout[3])
            cout = out[1]
        gain[ic] = np.median(hg)
        if ic<8:
            array[xs:xe,ic*chany:(ic+1)*chany]/=gain[ic]
        else:
            array[xs:xe,(ic-8)*chany:(ic-7)*chany]/=gain[ic]

    # then compare 0-8, 1-9, ..., 7-15 and calib all channels with respect to channel 0
    vy=np.array([])
    for ic in np.arange(8):
        vy=np.append(vy,getsafevy(ic,ny))
    vy = vy[::step]
    vx = np.append(np.arange(nx/2-bx-lfit,nx/2-bx),np.arange(nx/2+bx,nx/2+bx+lfit))
    hg = np.array([])
    pinit = [1, 1, 1, 1.0]
    for iy in vy:
        vax = array[vx,iy]
        vae = np.sqrt(vax)
        out = optimize.leastsq(errfuncy, pinit,
                               args=(vx, vax, vae), full_output=1)
        pout = out[0]
        hg = np.append(hg,pout[3])
        cout = out[1]
    g = np.median(hg)
    #print 'up,down g = %5.3f +/- %5.3f' %(np.median(hg),np.std(hg))
    gain[8:]*=g
    return gain
    

def Smooth_and_Calib(img):
     chan_gain = CalibGainPerChannel(img)
     calib_img = calibgain(img,chan_gain)
     new_img = ndimage.gaussian_filter(calib_img, sigma=(9, 9), order=0)
     chan_gain = CalibGainPerChannel(new_img)
     calib_img = calibgain(new_img,chan_gain)
     new_img = ndimage.gaussian_filter(calib_img, sigma=(9, 9), order=0)
     return new_img, chan_gain


