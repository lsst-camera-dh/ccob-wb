import os
import glob
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import lsst.eotest.raft as raft
import matplotlib.pyplot as plt
import numpy as np
import yaml
import sys
import ccob_utils as u
import ccob_beam as b
import pickle



def make_all_beams(config_file_beam):
    config = u.load_ccob_config(config_file_beam)
    led_names = ['nm960', 'nm850', 'nm750', 'red', 'blue', 'uv']
    ref_amp = 13

    for led in led_names:
        config['led_name'] = led
        beam = b.CcobBeam(config)
        beam.recons(ref_slot='11', ref_amp=ref_amp, ref_pix_x=1000,ref_pix_y=256)
        beam.make_image()
        beam.find_max()
        filename = beam.config['led_name']+'_beam_slot'+beam.properties['ref_slot']+'_amp'+str(ref_amp)+\
                   '_refx'+str(beam.properties['ref_pix_x'])+'_refy'+str(beam.properties['ref_pix_y'])+'.pkl'
        print('Printing to ',filename)
        beam.save(os.path.join(beam.config['tmp_dir'],filename))

def load_flat(config_file_data, config_file_beam):
    
    #config_file_data = '../ccob_config_RTM-006.yaml'
    config_data = u.load_ccob_config(config_file_data)
    config_beam = u.load_ccob_config(config_file_beam)
    
    slot = config_beam['slot']
    file_list=sorted(u.find_files(config, slot=slot))
 
    mean_slot_file = slot+'_mean_ccob_image.fits'
    imutils.fits_mean_file(file_list, os.path.join(config['tmp_dir'],mean_slot_file))

    fits_file = os.path.join(config_data['tmp_dir'],mean_slot_file)
    gains_dict={}
    ccd_dict={}

    #bias_frames = glob.glob(os.path.join(config['path'], slot+'_bias*'))
    mean_bias_file = slot+'_mean_bias_image_RTM-006_new.fits'
    #imutils.fits_mean_file(bias_frames, os.path.join(config['tmp_dir'],mean_bias_file))
    ccd_dict = sensorTest.MaskedCCD(fits_file, bias_frame=os.path.join(config['tmp_dir'],mean_bias_file))
    eotest_results_file = os.path.join(config_data['eo_data_path'], '{}_eotest_results.fits'.format(ccd_dict.md('LSST_NUM')))
    gains_dict = u.gains(eotest_results_file)
    
    mosaic, amp_coord = u.make_ccd_2d_array(fits_file, gains=gains_dict)
    return mosaic, amp_coord