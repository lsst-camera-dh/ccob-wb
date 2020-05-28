import sys
import os
import glob
import numpy as np
import pickle as pkl
import ccob_beam as beam
import ccob_qe_analysis as qe
import ccob_utils as u

 
def main():
    # print command line arguments
    config = sys.argv[0]
    print(config)
    config = u.load_ccob_config("beam_config.yaml")
    dirlist = []
    for d in config['rundir']:
        dirlist += glob.glob(config['rootdir']+d+'ccob_'+config['led_name']+'*')
    
    print(len(dirlist))
    assert (len(dirlist)/config['scan_size']).is_integer, 'Wrong number of scan locations'
        
    print(f'Scan over {len(dirlist)} locations')
  

    b = beam.CcobBeam(config)

    for i in np.arange(config['scan_size']): 
        # that's just to allow saving the data at every line of the scan in case
        # the job is interrupted.
        start = i*config['scan_size']
        end = (i+1)*config['scan_size']
        b.read_multibunch(config, dirlist=dirlist[start:end])
        b.save(os.path.join(config['tmpdir'],'beam_object_'+config['ref_raft']+'_'
                            +config['ref_slot']+'_'+config['led_name']+'.pkl'))
        
if __name__ == "__main__":
    main()
