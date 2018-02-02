
fiducialbeam_halfx = 19.5
fiducialbeam_halfy = 20.22

def is_ccd_in_beam(xb, yb, iccd):
    # Get position of the centroid of the CCD
    x0, y0 = GetCCDPos(iccd)

    # size (in mm) of ccd pixels = (40.00, 40.72) ITL, (40.04, 40.96) E2V
    ccd_xsize, ccd_ysize = GetCCDSize(iccd)
    
    lengthx = ccd_xsize/2+fiducialbeam_halfx
    lengthy = ccd_ysize/2+fiducialbeam_halfy
    if abs(xb-x0)<=lengthx and abs(yb-y0)<=lengthy:
        return True
    else:
        return False
