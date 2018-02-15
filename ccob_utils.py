import os
import glob
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import lsst.eotest.raft as raft
import numpy as np
import yaml

def gains(eotest_results_file):
    """
    Extract Fe55 gains from the results file of some eo testing.
    """
    results = sensorTest.EOTestResults(eotest_results_file)
    return {amp: gain for amp, gain in zip(results['AMP'], results['GAIN'])}

def load_ccob_config(config_file):
    """
    Loads ccob configuration (led, current, exp_time and position)
    from a config yaml file
    """
    config = yaml.load(open(config_file))
    if 'path' not in config.keys(): config['path'] = './'
    if 'led_name' not in config.keys(): config['led_name'] = '*'
    if 'current' not in config.keys(): config['current'] = '*'
    if 'exp_time' not in config.keys(): config['exp_time'] = '*'
    if 'xpos' not in config.keys(): config['xpos'] = '*'
    if 'ypos' not in config.keys(): config['ypos']= '*'
    if 'slot' not in config.keys(): config['slot'] ='*'
    
    return config

def find_files(config):
    """
    Find all the files matching a given ccob configuration
    """
    f_pattern = os.path.join(config['path'], config['slot'] + '_CCOB_' + config['led_name'] + '_'
                             + config['current'] + '_' + config['exp_time'] + '_X'
                             + config['xpos'] + '_Y' + config['ypos'] + '*')
    print f_pattern
    return glob.glob(f_pattern)

def build_mean_bias_frame(config, slot, mean_frame_pattern='_mean_bias_image.fits'):
    """
    Builds and save mean bias frames. Only need to do it once for each slot.
    make_image will look into config['tmp_dir'] to find them.
    """
    bias_frames = glob.glob(os.path.join(config['path'], slot+'_Bias*'))
    mean_bias_file = slot + mean_frame_pattern
    imutils.fits_mean_file(bias_frames, os.path.join(config['tmp_dir'],mean_bias_file))
    
def define_ccd_pos(ccd_pos_dict, raft_name, slot_names, xpos, ypos):
    """
    Updates a given dictionary of raft/slots center positions with that of a new raft.
    Inputs: raft_name (str), slot_names (list of str), xpos, ypos (lists of floats)           
    """
    ccd_pos_dict[raft_name] = {slot:[xpos[i],ypos[i]] for i,slot in enumerate(slot_names)}    
    
def make_image(config, slot_names, mean_frame_pattern='_mean_bias_image.fits'):
    """
    Make the mosaic image of the entire raft, when illuminated by the CCOB 
    according to config. Returns:
    - the raw image of the raft
    - the corrected image where mean bias frame has been removed and gains applied
    """
    file_list = sorted(find_files(config))
    fits_files_dict = {slot_names[i] : file_list[i] for i in range(len(file_list))}
    ccd_dict = {}
    ccd_dict_wbias = {}
    gains_dict = {}
    for slot in slot_names:

        mean_bias_file = slot + mean_frame_pattern

        ccd_dict[slot] = sensorTest.MaskedCCD(fits_files_dict[slot])
        outfile = os.path.join(config['tmp_dir'],'ccd' + slot + '.fits')
        image={}
     
        ccd_dict_wbias[slot]=sensorTest.MaskedCCD(fits_files_dict[slot],\
                                                 bias_frame=os.path.join(config['tmp_dir'],mean_bias_file))
        outfile_wbias = os.path.join(config['tmp_dir'],'ccd' + slot + '_wbias.fits')
        image_wbias={}
        
        eotest_results_file = os.path.join(config['eo_data_path'],'{}_eotest_results.fits'.format(ccd_dict[slot].md('LSST_NUM')))
        gains_dict[slot] = gains(eotest_results_file)
        
        for amp in ccd_dict[slot]:
            image[amp] = ccd_dict[slot].bias_subtracted_image(amp)
            image[amp] *= gains_dict[slot][amp]
            image_wbias[amp] = ccd_dict_wbias[slot].bias_subtracted_image(amp)
            image_wbias[amp] *= gains_dict[slot][amp]
      
        imutils.writeFits({amp: image_wbias[amp].getImage() for amp in ccd_dict_wbias[slot]}, 
                          outfile_wbias, fits_files_dict[slot])
        imutils.writeFits({amp: image[amp].getImage() for amp in ccd_dict[slot]}, 
                          outfile, fits_files_dict[slot])

    fits_files_dict_corr={slot : os.path.join(config['tmp_dir'],'ccd'+slot+'.fits') for slot in slot_names}
    fits_files_dict_corr_wbias={slot : os.path.join(config['tmp_dir'],'ccd'+slot+'_wbias.fits') for slot in slot_names}
    
    im_corr = raft.RaftMosaic(fits_files_dict_corr, bias_subtract=False)
    im_corr_wbias = raft.RaftMosaic(fits_files_dict_corr_wbias, bias_subtract=False)
    im_raw = raft.RaftMosaic(fits_files_dict, bias_subtract=False)
    
    return im_raw, im_corr, im_corr_wbias
