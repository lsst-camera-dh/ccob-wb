import os
import glob
import warnings
import lsst.eotest.image_utils as imutils
import lsst.eotest.fitsTools as fitsTools
import lsst.eotest.sensor as sensorTest
import lsst.eotest.raft as raft
from lsst.eotest.sensor.AmplifierGeometry import parse_geom_kwd
import numpy as np
import yaml
import astropy.io.fits as fits
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning
import pdb

def gains(eotest_results_file, is_PTC=False):
    """
    Extract Fe55 gains from the results file of some eo testing.
    """
    results = sensorTest.EOTestResults(eotest_results_file)
    if is_PTC:
        return {amp: gain for amp, gain in zip(results['AMP'], results['PTC_GAIN'])}
    else:
        return {amp: gain for amp, gain in zip(results['AMP'], results['GAIN'])}


def load_ccob_config(config_file):
    """
    Loads ccob configuration (led, current, exp_time and position)
    from a config yaml file
    """
    config = yaml.load(open(config_file), Loader=yaml.FullLoader)
    if 'path' not in config.keys(): config['path'] = './'
    if 'led_name' not in config.keys(): config['led_name'] = '*'
    if 'current' not in config.keys(): config['current'] = '*'
    if 'exp_time' not in config.keys(): config['exp_time'] = '*'
    if 'xpos' not in config.keys(): config['xpos'] = '*'
    if 'ypos' not in config.keys(): config['ypos']= '*'
    if 'slot' not in config.keys(): config['slot'] ='*'
    
    return config

def find_files(config, slot='*'):
    """
    Find all the files matching a given ccob configuration
    """
    f_pattern = os.path.join(os.path.join(config['path'],config['led_name']), slot+'*' + config['led_name'] + '*'
                             + config['current'] + '*' + config['exp_time'] + '*'
                             + config['xpos'] + '*' + config['ypos'] + '*')
    print(f_pattern)
    return glob.glob(f_pattern)

def build_mean_bias_frame(config, slot, file_patterm='_mean_bias_image.fits'):
    """
    Builds and save mean bias frames. Only need to do it once for each slot.
    make_image will look into config['tmp_dir'] to find them.
    """
    bias_frames = glob.glob(os.path.join(config['path'], slot+'_Bias*'))
    outfile = slot + file_pattern
    imutils.fits_mean_file(bias_frames, os.path.join(config['tmp_dir'],outfile))

def make_superbias_frame(config, bias_files, raft, slot, file_pattern='_sbias_image.fits'):
    """
    Make and save a super biasframes. Only need to do it once for each slot.
    make_image will look into config['tmp_dir'] to find them.
    """
    amp_geom = sensorTest.makeAmplifierGeometry(bias_files[0])
    outfile = raft+'_'+slot + file_pattern
    imutils.superbias_file(bias_files, amp_geom.serial_overscan, os.path.join(config['tmp_dir'],outfile)) 

def define_ccd_pos(ccd_pos_dict, raft_name, slot_names, xpos, ypos):
    """
    Updates a given dictionary of raft/slots center positions with that of a new raft.
    Inputs: raft_name (str), slot_names (list of str), xpos, ypos (lists of floats)           
    """
    ccd_pos_dict[raft_name] = {slot:[xpos[i],ypos[i]] for i,slot in enumerate(slot_names)}    

def glob_files(root_dir, raft, run, *args):
    return sorted(glob.glob(os.path.join(root_dir, raft, run, *args)))
        
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

def fill_seg_dict(ccd):
    seg_dict_slot = {}
    for amp in ccd.keys():
        datasec = [ccd.amp_geom[amp]['DATASEC'].strip('[]').split(',')[0].split(':')[0],
                   ccd.amp_geom[amp]['DATASEC'].strip('[]').split(',')[0].split(':')[1],
                   ccd.amp_geom[amp]['DATASEC'].strip('[]').split(',')[1].split(':')[0],
                   ccd.amp_geom[amp]['DATASEC'].strip('[]').split(',')[1].split(':')[1]]
        detsec = [ccd.amp_geom[amp]['DETSEC'].strip('[]').split(',')[0].split(':')[0],
                  ccd.amp_geom[amp]['DETSEC'].strip('[]').split(',')[0].split(':')[1],
                  ccd.amp_geom[amp]['DETSEC'].strip('[]').split(',')[1].split(':')[0],
                  ccd.amp_geom[amp]['DETSEC'].strip('[]').split(',')[1].split(':')[1]]
        detsize = [ccd.amp_geom[amp]['DETSIZE'].strip('[]').split(',')[0].split(':')[0],                      
                   ccd.amp_geom[amp]['DETSIZE'].strip('[]').split(',')[0].split(':')[1],
                   ccd.amp_geom[amp]['DETSIZE'].strip('[]').split(',')[1].split(':')[0],
                   ccd.amp_geom[amp]['DETSIZE'].strip('[]').split(',')[1].split(':')[1]]
     
        seg_dict_slot[amp] = {'datasec':list(map(int,datasec)),
                         'detsec':list(map(int,detsec)),
                         'detsize':list(map(int,detsize))
                        }
    return seg_dict_slot

def make_ccd_image(seg_dict, ccd_dict, slot):
    img_tot = np.empty(shape=(seg_dict[slot][1]['detsize'][3],seg_dict[slot][1]['detsize'][1]))
    for seg in seg_dict[slot].keys():
        xmin = seg_dict[slot][seg]['datasec'][0]
        xmax = seg_dict[slot][seg]['datasec'][1]
        ymin = seg_dict[slot][seg]['datasec'][2]
        ymax = seg_dict[slot][seg]['datasec'][3]
        xmin_detset = seg_dict[slot][seg]['detsec'][0]
        xmax_detset = seg_dict[slot][seg]['detsec'][1]
        ymin_detset = seg_dict[slot][seg]['detsec'][2]
        ymax_detset = seg_dict[slot][seg]['detsec'][3]
    
        if (xmin_detset < xmax_detset) & (ymin_detset < ymax_detset):
                img_tot[ymin_detset-1:ymax_detset-1,xmin_detset-1:xmax_detset-1] = \
                  ccd_dict[slot][seg].getArrays()[0][ymin-1:ymax-1,xmin-1:xmax-1]

        if (xmin_detset < xmax_detset) & (ymin_detset > ymax_detset):
                img_tot[ymax_detset-1:ymin_detset-1,xmin_detset-1:xmax_detset-1] =\
                  ccd_dict[slot][seg].getArrays()[0][ymin-1:ymax-1,xmin-1:xmax-1][::-1,:]

        if (xmin_detset > xmax_detset) & (ymin_detset < ymax_detset):
                img_tot[ymin_detset-1:ymax_detset-1,xmax_detset-1:xmin_detset-1] =\
                  ccd_dict[slot][seg].getArrays()[0][ymin-1:ymax-1,xmin-1:xmax-1][:,::-1]

        if (xmin_detset > xmax_detset) & (ymin_detset > ymax_detset):
                img_tot[ymax_detset-1:ymin_detset-1,xmax_detset-1:xmin_detset-1] =\
                  ccd_dict[slot][seg].getArrays()[0][ymin-1:ymax-1,xmin-1:xmax-1][::-1,::-1]

    return img_tot

def make_ccd_2d_array(infile, biasfile=None, gains=None):
    '''
    Generate a 2d array of a sensor image (trimmed, bias subtracted) from 
    an input fits file and a gain dictionary. 
    Function adapted from sensorTest.plot_flat()
    '''
    ccd = sensorTest.MaskedCCD(infile, bias_frame=biasfile)
    foo = fits.open(infile)
    datasec = parse_geom_kwd(foo[1].header['DATASEC'])
    # Specialize to science sensor or wavefront sensor geometries.
    nx_segments = 8
    ny_segments = len(ccd)//nx_segments
    nx = nx_segments*(datasec['xmax'] - datasec['xmin'] + 1)
    ny = ny_segments*(datasec['ymax'] - datasec['ymin'] + 1)
    mosaic = np.zeros((ny, nx), dtype=np.float)
    amp_coords = {}
    for ypos in range(ny_segments):
        for xpos in range(nx_segments):
            amp = ypos*nx_segments + xpos + 1
            #
            # Determine subarray boundaries in the mosaicked image array
            # from DETSEC keywords for each segment.
            detsec = parse_geom_kwd(foo[amp].header['DETSEC'])
            datasec = parse_geom_kwd(foo[amp].header['DATASEC'])
            xmin = nx - max(detsec['xmin'], detsec['xmax'])
            xmax = nx - min(detsec['xmin'], detsec['xmax']) + 1
            ymin = ny - max(detsec['ymin'], detsec['ymax'])
            ymax = ny - min(detsec['ymin'], detsec['ymax']) + 1
            #
            #
            # Extract the bias-subtracted masked image for this segment.
            segment_image = ccd.unbiased_and_trimmed_image(amp)
            subarr = segment_image.getImage().getArray()
            #
            # Determine flips in x- and y-directions in order to
            # get the (1, 1) pixel in the lower right corner.
            flipx = False
            flipy = False
            if detsec['xmax'] > detsec['xmin']:  # flip in x-direction
                subarr = subarr[:, ::-1]
                flipx = True
            if detsec['ymax'] > detsec['ymin']:  # flip in y-direction
                subarr = subarr[::-1, :]
                flipy = True
            #
            # Convert from ADU to e-
            if gains is not None:
                subarr *= gains[amp]
            #
            # Save coordinates of segment for later use
            amp_coords[(xpos, ypos)] = {'amp':amp,
                                        'segment':foo[amp].header['EXTNAME'], 
                                        'xmin':xmin, 
                                        'xmax':xmax, 
                                        'ymin':ymin, 
                                        'ymax':ymax, 
                                        'flipx':flipx, 
                                        'flipy':flipy,
                                        'detsec':detsec,
                                        'datasec':datasec}
            # Set the subarray in the mosaic.
            mosaic[ymin:ymax, xmin:xmax] = subarr
    return mosaic, amp_coords

def pix_coord_in_mosaic(amp_coord, amp, posx_in_amp, posy_in_amp):
    '''
    Given a pixel position in a given segment of a given sensor 
    and a amp_coords dictionary (obtained from make_ccd_2d_array)
    returns the pixel location in the full mosaicked image of the sensor.
    '''
    subdict = {k:v for k, v in amp_coord.items() if v['amp'] == amp}    
    v = subdict[list(subdict.keys())[0]]
    idx_x = v['xmin'] + posx_in_amp
    idx_y = v['ymin'] + posy_in_amp
    if v['flipx'] == True: # flip-x
        idx_x = v['xmax'] - posx_in_amp
    if v['flipy'] == True: # flip-y
        idx_y = v['ymax'] - posy_in_amp
    return idx_x, idx_y


def writeFits_from_dict(amp_dict, outfile, template_file, bitpix=32):
    '''
    Same as eotest imutils writeFits but takes a dictionary of amplifier as input
    rather than a list of afwImage images
    '''
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    all_amps = imutils.allAmps()
    for amp in all_amps:
        if bitpix < 0:
            output.append(fits.ImageHDU(data=amp_dict[amp]))
        else:
            output.append(fits.CompImageHDU(data=amp_dict[amp],
                                            compression_type='RICE_1'))
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyUserWarning,
                                append=True)

    with fits.open(template_file) as template:
        output[0].header.update(template[0].header)
        output[0].header['FILENAME'] = outfile
        for amp in all_amps:
            output[amp].header.update(template[amp].header)
            imutils.set_bitpix(output[amp], bitpix)
#            print(np.median(output[amp].data.ravel()))
        for i in (-3, -2, -1):
            output.append(template[i])
        imutils.fitsWriteto(output, outfile, overwrite=True, checksum=True)


