import os
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import numpy as np
import ccob_utils as u
import ccob_beam as b
import pickle
import matplotlib.pyplot as plt
import pdb 
from lsst.obs.lsst.cameraTransforms import LsstCameraTransforms
from lsst.obs.lsst import LsstCamMapper as camMapper


def load_ccob_beam(path_to_beam='/home/combet/tmp_9rafts/60x60' ,led_name='red', ref_amp=5, ccdid='R22_S11', ref_pix_x=1000,ref_pix_y=256):
    """
    Load a CCOB beam object from a scan, find the maximum

    Parameters
    ----------
        led_name: string
            Choice of CCOB LED
        ref_amp: int
            Amplifier where the bunch of pixels is located
        ref_slot: string
            Sensor where the bunch of pixels is located
        ref_pix_x: int
            x position in pixel coordinates where the bunch of pixels is located
        ref_pix_y:
            y position in pixel coordinates where the bunch of pixels is located

    """

    fl = glob.glob(os.path.join(path_to_beam),'*'+ccdid+'_'+led_name+'*')
    beam_file = fl[0]
    b = pkl.load(open(beam_file,'rb'))
    b.interp_beam_BOT(amp=ref_amp, pd_corr=True)
    b.make_image_BOT()
    b.find_max_from_avg()
    return b


def amp_to_seg_dict():

    d = {
        1:'10',
        2:'11',
        3:'12',
        4:'13',
        5:'14',
        6:'15',
        7:'16',
        8:'17', 
        9:'07',
        10:'06',
        11:'05',
        12:'04',
        13:'03',
        14:'02',
        15:'01',
        16:'00'
    }
    return d

def compute_offsets(beam, lct, ccdid='R22_S11'):
    
    det = lct.getDetector(ccdid)
    tmp_ccdid, xmax, ymax = lct.focalMmToCcdPixel(bb.properties['max_yccob'], bb.properties['max_xccob'])
    if 'E2V' in det.getSerial():
        width_seg = 512
        length_seg = 2002
    else:
        width_seg = 508
        length_seg = 2000

    ref_pix_x =  + beam.properties['ref_pix_x']
    ref_pix_y =  + beam.properties['ref_pix_y']
        
    delta_x = ref_pix_x - xmax
    delta_y = ref_pix_y -ymax
    return delta_x, delta_y # in FP coordinates

def load_ccd(self, led_name='red'):
    """
    Load and generate mosaic image from a given sensor illuminated 
    by the CCOB. The path to the data is provided in the self.config_file_data file.

    Parameters
    ----------
        led_name: choice of the CCOB LED. Either one of ['nm960','nm850','nm750,'red','blue,'uv'] 
    """

    #config_file_data = '../ccob_config_RTM-006.yaml'
#        config_data = u.load_ccob_config(self.config_file_data)
#        config_beam = u.load_ccob_config(self.config_file_beam)

    self.config_data['led_name'] = led_name
    slot = self.beam.properties['ref_slot']
    file_list=sorted(u.find_files(self.config_data, slot=slot))

    mean_slot_file = slot+'_mean_ccob_image.fits'
    imutils.fits_mean_file(file_list, os.path.join(self.config_data['tmp_dir'],mean_slot_file))

    fits_file = os.path.join(self.config_data['tmp_dir'],mean_slot_file)
    gains_dict={}
    ccd_dict={}

    #bias_frames = glob.glob(os.path.join(config['path'], slot+'_bias*'))
    #mean_bias_file = slot+'_mean_bias_image_RTM-006_new.fits'
    #imutils.fits_mean_file(bias_frames, os.path.join(config['tmp_dir'],mean_bias_file))
    #ccd_dict = sensorTest.MaskedCCD(fits_file, bias_frame=os.path.join(self.config_data['tmp_dir'],mean_bias_file))

    superbias_frame = make_superbias_frame(self.config_data, slot=slot)
    ccd_dict = sensorTest.MaskedCCD(fits_file, bias_frame=os.path.join(self.config_data['tmp_dir'],superbias_frame))

    eotest_results_file = os.path.join(self.config_data['eo_data_path'],\
                                       '{}_eotest_results.fits'.format(ccd_dict.md('LSST_NUM')))

    gains_dict = u.gains(eotest_results_file)
    self.ccd = {}
    self.ccd['mosaic'], self.ccd['amp_coord'] = u.make_ccd_2d_array(fits_file, gains=gains_dict)
    self.ccd['xpos_ccob'] = self.config_data['xpos']
    self.ccd['ypos_ccob'] = self.config_data['ypos']

    return self.ccd

def compute_QE(self):
    """
    Computes the mosaicked QE image from the reconstructed beam image (self.beam) and the 
    CCOB-illuminated sensor image (self.ccd).

    First, the beam model image is matched to the data, given the position of the CCOB when the data
    were taken. Then the ratio data/beam produced the CCOB flat field from which the relative QE may be measured.
    """

    # First need to match a beam image to the ccd data, given the position of the CCOB
    ref_pix = u.pix_coord_in_mosaic(self.ccd['amp_coord'], self.beam.properties['ref_amp'], \
                                    self.beam.properties['ref_pix_y'],
                                    self.beam.properties['ref_pix_x'])

    delta_x = float(self.ccd['ypos_ccob']) - self.beam.properties['max_yccob']
    delta_y = float(self.ccd['xpos_ccob']) - self.beam.properties['max_xccob']
    delta_x_pix = int(delta_x/0.01)
    delta_y_pix = int(delta_y/0.01)

    print('deplacement in mm: dx=%0.2f dy=%0.2f'%(delta_x, delta_y))
    print('deplacement in pixels: dx=%i dy=%i'%(delta_x_pix, delta_y_pix))

    geom_center_pos=(ref_pix[0]+delta_x_pix, ref_pix[1]+delta_y_pix)

    # distance from beam center to ccd edges in mm
    dist_to_left = geom_center_pos[0]*0.01
    dist_to_bottom = geom_center_pos[1]*0.01
    dist_to_right = (self.ccd['mosaic'].shape[1]-geom_center_pos[0])*0.01
    dist_to_top = (self.ccd['mosaic'].shape[0]-geom_center_pos[1])*0.01

    # bounding box to use to reconstruct the CCOB with dimensions matching ccd
    bbox = (self.beam.properties['max_xccob'] - dist_to_left,
    self.beam.properties['max_xccob'] + dist_to_right,
    self.beam.properties['max_yccob'] - dist_to_bottom,
    self.beam.properties['max_yccob'] + dist_to_top)

    xarr = np.linspace(bbox[0],bbox[1],self.ccd['mosaic'].shape[1])
    yarr = np.linspace(bbox[2],bbox[3],self.ccd['mosaic'].shape[0])
    self.beam_image = self.beam.beam_image['f_interp'](xarr, yarr)

    # Raw QE (flat) is simply the ccd-to-beam ratio
    self.QE = self.ccd['mosaic']/self.beam_image
    return self.QE

def plot_QE(self):
    plt.imshow(self.QE.QE, vmin=0.695, vmax=0.715, cmap='hot')
    plt.colorbar()
    plt.show()

def make_fits(self, outfile, template_file):
    """
    Saves the mosaicked QE image into a fits file, using the default format 

    Parameters
    ----------
        outfile: string
            Name of the output fits file into which the data will be saved
        template_file: string
            Name of the file to be used as template in the writeFits function

    """
    amp_dict = {}

    for i,amp_pos in enumerate(self.ccd['amp_coord']):

        amp = self.ccd['amp_coord'][amp_pos]['amp']
        xmin = self.ccd['amp_coord'][amp_pos]['xmin']
        xmax = self.ccd['amp_coord'][amp_pos]['xmax']
        ymin = self.ccd['amp_coord'][amp_pos]['ymin']
        ymax = self.ccd['amp_coord'][amp_pos]['ymax']
        flipx = self.ccd['amp_coord'][amp_pos]['flipx']
        flipy = self.ccd['amp_coord'][amp_pos]['flipy']
        detsec = self.ccd['amp_coord'][amp_pos]['detsec']
        datasec = self.ccd['amp_coord'][amp_pos]['datasec']

        foo = self.QE[ymin:ymax,xmin:xmax]

        if flipx:
            amp_dict[amp] = foo[:,::-1]
        elif flipy:
            amp_dict[amp] = foo[::-1,:]
        else:
            amp_dict[amp] = foo

    ccd_dict = sensorTest.MaskedCCD(template_file)
    shape = ccd_dict[1].getImage().getArray().shape        

    amp_dict_w_overscan = {}
    for amp in amp_dict:
        arr = np.zeros(shape)
        arr[datasec['ymin']-1:datasec['ymax'],datasec['xmin']-1:datasec['xmax']] = amp_dict[amp]
        amp_dict_w_overscan[amp] = arr

    u.writeFits_from_dict(amp_dict_w_overscan, outfile, template_file, bitpix=-32)
    
    
    
if __name__ == '__main__':

#    led_names = ['nm850', 'nm750', 'red']
#    config_file_beam = 'ccob_beam_recons_config.yaml'
#    config_file_data = 'ccob_config_RTM-006.yaml'

#    xpos = [253.0, 253.0, 253.0, 295.0, 295.0, 295.0, 337.0, 337.0, 337.0]
#    ypos = [237.0, 195.0, 153.0, 237.0, 195.0, 153.0, 237.0, 195.0, 153.0]
#    slot_names=['00','01','02','10','11','12','20','21','22']
#    ccd_pos_dict={}
#    u.define_ccd_pos(ccd_pos_dict, 'RTM-006', slot_names, xpos, ypos)        


    led_names = ['red']
    config_file_beam = 'ccob_beam_recons_config.yaml'
    config_file_data = 'ccob_config_RTM-006.yaml'

    xpos = [253.0]
    ypos = [237.0]
    slot_names=['00']
    ccd_pos_dict={}
    u.define_ccd_pos(ccd_pos_dict, 'RTM-006', slot_names, xpos, ypos)        
  
    
    for led in led_names:

        for slot in ccd_pos_dict['RTM-006']:
            print(slot)
            QE = qe.CcobQE(config_file_beam, config_file_data)
            QE.make_ccob_beam(led_name=led, ref_amp=13, ref_slot=slot, ref_pix_x=1000,ref_pix_y=256)
            print(QE.beam.properties)

            QE.config_data['xpos'] = str(ccd_pos_dict['RTM-006'][slot][0])
            QE.config_data['ypos'] = str(ccd_pos_dict['RTM-006'][slot][1])
            ccd = QE.load_ccd(led_name=led)

            QE.compute_QE()
            template_path = '/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive-test/LCA-11021_RTM/LCA-11021_RTM-006-Dev/5867D/qe_raft_acq/v0/38892/S'+slot
            template_file = glob.glob(os.path.join(os.path.join(QE.config_data['path'],led),slot+'*'))[0]
            outfile = os.path.join(QE.config_data['tmp_dir'],'QE_S' + slot + '_' + led
                                   + '_amp'+str(QE.beam.properties['ref_amp'])
                                   + '_refx'+str(QE.beam.properties['ref_pix_x'])
                                   + '_refy'+str(QE.beam.properties['ref_pix_y'])+'.fits')
            QE.make_fits(outfile, template_file)        

