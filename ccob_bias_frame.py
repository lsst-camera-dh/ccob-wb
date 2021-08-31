"""
All that is needed to reconstruct the CCOB beam and generate the beam model
used to produce the synthetic flat fields.
"""

import os
import glob
import numpy as np
import ccob_utils as u

def main():
    """
    First, finds the list of files relevant 
    """
    
basedir1 = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_stage-test/LCA-10134_Cryostat/LCA-10134_Cryostat-0001/'
basedir2 = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_stage/LCA-10134_Cryostat/LCA-10134_Cryostat-0001/'

path_to_data = {'R01': basedir1 + '6848D/BOT_acq/v0/48087/',
                'R02': basedir1 + '6849D/BOT_acq/v0/48093/',
                'R10': None,
                'R11': basedir1 + '6851D/BOT_acq/v0/48108/',
                'R12': basedir1 + '6852D/BOT_acq/v0/48113/',
                'R20': basedir1 + '6853D/BOT_acq/v0/48118/',
                'R21': None,
                'R22': basedir2 + '11974/BOT_acq/v0/93868/',
                'R30': basedir1 + '6843D/BOT_acq/v0/48047/'
               }

outdir = '/gpfs/slac/lsst/fs1/u/combet/DATA/CCOB_Bias/Run3'
sensor_list = ['S00','S01','S02','S10','S11','S12','S20','S21','S22']
raft_list = ['R01', 'R02', 'R11', 'R12', 'R20', 'R21', 'R22', 'R30']

for raft in ['R01', 'R10, ''R11', 'R12', 'R20', 'R30']:
    for sensor in sensor_list:
        ccdid = raft + '_' + sensor
        print(ccdid)
        bias_files = glob.glob(f"{path_to_data[raft]}/bias_bias*/*{ccdid}*")
#        print(sorted(bias_files))
        u.make_superbias_frame(raft, sensor, bias_files, outdir, file_pattern='_sbias.fits')

      
if __name__ == '__main__':
    main()
 