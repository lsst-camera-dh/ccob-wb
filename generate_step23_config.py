import numpy as np
from lsst.obs.lsst.cameraTransforms import LsstCameraTransforms
from lsst.obs.lsst import LsstCamMapper as camMapper

####################################################################################
# This script generates the configuration file for step22 of the CCOB_WB acquisition
# The output file name is on line 53.
# The offsets can be adjusted below
####################################################################################

# Offsets
xoffset = 5.0 # mm
yoffset = -7.0 # mm

# LED setup: current and exposure times
led = ['uv', 'blue', 'red', 'nm750', 'nm850', 'nm960']
current = [0.1, 0.006, 0.009, 0.01, 0.01, 0.01]
exp_time_orig = [0.02, 0.05, 0.05, 0.09, 0.06, 0.12]
exp_time_corr_factor = [1.33, 1.69, 1.46, 1.4, 1.5, 2.7]
exp_time = np.array(exp_time_orig) * np.array(exp_time_corr_factor)
led_config = list(zip(led, current, np.around(exp_time,2)))
print(led_config)

# Full focal plane pointings.
# 5 pointings/CCD. One in the center and one in the center of each quadrant
camera = camMapper._makeCamera()
lct = LsstCameraTransforms(camera)
raft_list = ['R01','R02','R03',
             'R10','R11','R12','R13','R14',
             'R20','R21','R22','R23','R24',
             'R30','R31','R32','R33','R34',
             'R41','R42','R43']
sensor_list = ['S00','S01','S02','S10','S11','S12','S20','S21','S22']
fp_pos_e2v = [(1024,3003),(3072,3003),(2048,2002),(1024,1001),(3072,1001)] 
fp_pos_itl = [(1018,3000),(3054,3000),(2036,2000),(1018,1000),(3054,1000)] 

xall=np.zeros((len(raft_list),len(sensor_list)*len(fp_pos_e2v)))
yall=np.zeros((len(raft_list),len(sensor_list)*len(fp_pos_e2v)))

for i,raft in enumerate(raft_list):
    for j,sensor in enumerate(sensor_list):
        ccdid = raft + '_' + sensor
        det = lct.getDetector(ccdid)
        fp_pos = fp_pos_e2v
        if 'ITL' in det.getSerial(): fp_pos = fp_pos_itl
        for k,pos in enumerate(fp_pos):
            ycenter, xcenter = lct.ccdPixelToFocalMm(pos[0],pos[1],ccdid)
            l = j*len(fp_pos)+k
            xall[i][l] = xcenter
            yall[i][l] = ycenter

# Print into configuration file
config_filename = 'ccob_qe_all_byREB.cfg'
print(config_filename)

with open(config_filename, 'w') as ostr:
    print('[ACQUIRE]', file=ostr)
    print('BIAS', file=ostr)
    print('CCOB', file=ostr)
    print('BIAS\n', file=ostr)
 

    print('[BIAS]', file=ostr)
    print('COUNT=10\n', file=ostr)

    print('[CCOB]', file=ostr)
    print('BCOUNT = 0', file=ostr)
    print('IMCOUNT = 5', file=ostr)
    print('XOFFSET = '+str(xoffset), file=ostr)
    print('YOFFSET = '+str(yoffset)+'\n', file=ostr)

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

    for ii,raft in enumerate(raft_list):
        np.savetxt('tmp.txt', list(zip(xall[ii], yall[ii])),fmt='%8.3f', delimiter=('     '))
        with open('tmp.txt', 'r') as istr:
            for i,line in enumerate(list(istr)):
                if i==0 and ii==0: 
                    line = 'point =  ' + line.rstrip('\n') + '     ' +raft+'/Reb0,'
                elif i<15:
                    line = '         ' + line.rstrip('\n') + '     ' +raft+'/Reb0,'
                elif i>= 15 and i<30:
                    line = '         ' + line.rstrip('\n') + '     ' +raft+'/Reb1,'
                elif i>=30 and i<44:
                    line = '         ' + line.rstrip('\n') + '     ' +raft+'/Reb2,'
                elif raft!='R43':
                    line = '         ' + line.rstrip('\n') + '     ' +raft+'/Reb2,'
                else:
                    line = '         ' + line.rstrip('\n') + '     ' +raft+'/Reb2'
                print(line, file=ostr)

