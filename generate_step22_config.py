import numpy as np
from lsst.obs.lsst.cameraTransforms import LsstCameraTransforms
from lsst.obs.lsst import LsstCamMapper as camMapper

####################################################################################
# This script generates the configuration file for step22 of the CCOB_WB acquisition
# The output file name is on line 63.
# The offsets can be adjusted below
####################################################################################

# Offsets
xoffset = 5.0 # mm
yoffset = -7.0 # mm

# Scan definition: Delta_x x Delta_y cm area, scanned using N_x x N_y points
Delta_x = 80 # mm; should not be changed
Delta_y = 80 # mm; should not be changed
N_x = 60 # change to 12 for the quick scan
N_y = 60 # change to 12 for the quick scam

# LED setup: current and exposure times
led = ['uv', 'blue', 'red', 'nm750', 'nm850', 'nm960']
current = [0.1, 0.006, 0.009, 0.01, 0.01, 0.01]
exp_time_orig = [0.02, 0.05, 0.05, 0.09, 0.06, 0.12]
exp_time_corr_factor = [1.33, 1.69, 1.46, 1.4, 1.5, 2.7]
exp_time = np.array(exp_time_orig) * np.array(exp_time_corr_factor)
led_config = list(zip(led, current, np.around(exp_time,2)))
print(led_config)

# Find the scan center, in camera coordinates
# Here, choosing the middle of Segment14 of R22_S11
raft = 'R22'
sensor = 'S11'
reb = 'Reb1'
ccdid = raft + '_' + sensor
camera = camMapper._makeCamera()
lct = LsstCameraTransforms(camera)
det = lct.getDetector(ccdid)
ref_sensor = det.getSerial()
print(ref_sensor)

if 'E2V' in ref_sensor:
    ycenter, xcenter = lct.ccdPixelToFocalMm(2304,3003,ccdid) #e2V
else:
    ycenter, xcenter = lct.ccdPixelToFocalMm(2286,3000,ccdid) #ITL
print(xcenter, ycenter)

# Scan boundaries
xmin = xcenter - Delta_x / 2.
xmax = xcenter + Delta_x / 2.
ymin = ycenter - Delta_y / 2.
ymax = ycenter + Delta_y / 2.
print(xmin, xmax, ymin, ymax)

# Scan positions, within the boundaries
xarr = np.linspace(xmin, xmax, N_x)
yarr = np.linspace(ymin, ymax, N_y)
xall = np.repeat(xarr, N_y)
yall = np.broadcast_to(yarr, (N_y, N_x)).flatten()
np.savetxt('tmp.txt', list(zip(xall, yall)),fmt='%.3f', delimiter=('  '))

# Print into configuration file
config_filename = 'ccob_'+str(N_x)+'x'+str(N_y)+'_'+raft+'_'+sensor+'_all.cfg'

with open('tmp.txt', 'r') as istr:
    with open(config_filename, 'w') as ostr:
        print('[ACQUIRE]', file=ostr)
        print('bias', file=ostr)
        print('ccob', file=ostr)
        print('bias\n', file=ostr)

        print('[BIAS]', file=ostr)
        print('COUNT=10', file=ostr)
        print('LOCATIONS = '+raft+'/'+reb+'\n', file=ostr)

        print('[CCOB]', file=ostr)
        print('BCOUNT = 0', file=ostr)
        print('IMCOUNT = 1', file=ostr)
        print('XOFFSET = '+str(xoffset), file=ostr)
        print('YOFFSET = '+str(yoffset), file=ostr)
        print('LOCATIONS = '+raft+'/'+reb+'\n', file=ostr)
        led_number = 0
        print('expose = '+ str(led_config[led_number][0])+' '
              + str(led_config[led_number][1])+' '+str(led_config[led_number][2])+',', file=ostr)
        led_number = 1
        print('         '+ str(led_config[led_number][0])+' '
              + str(led_config[led_number][1])+' '+str(led_config[led_number][2])+',', file=ostr)
        led_number = 2
        print('         '+ str(led_config[led_number][0])+' '
              + str(led_config[led_number][1])+' '+str(led_config[led_number][2])+',', file=ostr)
        led_number = 3
        print('         '+ str(led_config[led_number][0])+' '
              + str(led_config[led_number][1])+' '+str(led_config[led_number][2])+',', file=ostr)
        led_number = 4
        print('         '+ str(led_config[led_number][0])+' '
              + str(led_config[led_number][1])+' '+str(led_config[led_number][2])+',', file=ostr)
        led_number = 5
        print('         '+ str(led_config[led_number][0])+' '
              + str(led_config[led_number][1])+' '+str(led_config[led_number][2])+'\n', file=ostr)

        for i,line in enumerate(istr):
            if i==0: 
                line = 'point =  ' + line.rstrip('\n') + ','
            elif i!=len(xall)-1:
                line = '         ' + line.rstrip('\n') + ','
            else:
                line = '         ' + line.rstrip('\n')
            print(line, file=ostr)
