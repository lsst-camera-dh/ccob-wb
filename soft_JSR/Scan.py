import numpy as np

# scan step
step = 8 # mm

# To be define
raft_size =  127 # (mm)

# size of the beam that can be safely use to calib pixel (<1/1000 accuracy)
# need to be properly estimated with nbeam in beamreco.py

fiducialbeam_halfx = 19.5 # (mm)
fiducialbeam_halfy = 20.22 # (mm)

# origin can be 0,0 or FP center
x0 = 0 # mm
y0 = 0 # mm

# scan limits
lminx = -fiducialbeam_halfx
lmaxx = 5*raft_size+fiducialbeam_halfx
lminy = -fiducialbeam_halfy
lmaxy = 5*raft_size+fiducialbeam_halfy


# if we plan to save a bit of scan time we can avoid scaning empty corners
# and define ymin and ymax as below
def getylimits(x):
    ymin = lminy
    ymax = lmaxy
    if x<raft_size-fiducialbeam_halfx or x>4*raft_size+fiducialbeam_halfx:
        ymin = raft_size-fiducialbeam_halfy
        ymax = 4*raft_size+fiducialbeam_halfy
    return ymin,ymax

def CreatePositions(filename):
    vx = []
    vy = []

    # Fist positions are centered CCD
    for iccd in range(189):
        x,y = GetCCDPos(iccd)
        vx.append(x)
        vy.append(y)
    
    # Then Scaning of the whole FP
    for x in np.arange(lminx,lmaxx+0.1,step):
        ymin, ymax = getylimits(x)
        for y in np.arange(ymin,ymax):
            vx.append(x-x0)
            vy.append(y-y0)

    pos = np.array((vx,vy),dtype='f')
    np.savetxt(filename,pos.T,fmt='%.3f')



if __name__ == '__main__':
    filename = 'scan_positions.txt'
    CreatePositions(filename)
