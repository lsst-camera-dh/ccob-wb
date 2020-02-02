import os
import glob
import pickle as pkl
import scipy
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import numpy as np
import ccob_utils as u
import ccob_beam as b
import ccob_qe_data as qe_data
import pickle
import matplotlib.pyplot as plt
import pdb 
from lsst.obs.lsst.cameraTransforms import LsstCameraTransforms
from lsst.obs.lsst import LsstCamMapper as camMapper


def load_beam_model(path_to_beam, led_name='red', ref_amp=5, ccdid='R22_S11'):
    """
    Load a CCOB beam object from a scan, find the maximum

    Parameters
    ----------
        led_name: string
            Choice of CCOB LED
        ref_amp: int
            Amplifier where the bunch of pixels is located
        ccdid: string
            RXX_SYY - Raft/Sensor where the bunch of pixels is located
    """

    fl = glob.glob(os.path.join(path_to_beam,'*'+ccdid+'_'+led_name+'*'))
    beam_file = fl[0]
    b = pkl.load(open(beam_file,'rb'))
    b.interp_beam_BOT(amp=ref_amp, pd_corr=True)
    b.make_image_BOT()
    b.find_max_from_avg()
    return b

def compute_offsets(beam, lct, ccdid='R22_S11', ref_pix_x=2304, ref_pix_y=3003):
    '''    
    beam: beam object
    ccdid: raft_sensor from which the beam object was computed
    ref_pix_x, ref_pix_y: pixel in raft_sensor from which the beam was reconstructed. These are given in FP coordinates/convention - (0,0) is the bottom left corner of the sensor
    '''    

    det = lct.getDetector(ccdid)
    # FP coordinates of the reconstructed beam maximum
    tmp_ccdid, xmax, ymax = lct.focalMmToCcdPixel(beam.properties['max_yccob'], beam.properties['max_xccob'])
        
    delta_x = ref_pix_x - xmax
    delta_y = ref_pix_y - ymax
    return delta_x, delta_y # in FP coordinates


def define_model_bbox(beam_model, mosaic, lct, pos, delta_x, delta_y):
    
    xccob = float(pos.split('_')[0]) # CCS
    yccob = float(pos.split('_')[1]) # CCS
    ccdid, x_from_ccob_pos, y_from_ccob_pos = lct.focalMmToCcdPixel(yccob, xccob) # FP coord

    geom_center_pos = (np.shape(mosaic)[1]-(x_from_ccob_pos+delta_x),y_from_ccob_pos+delta_y) # EOtest coord

    dist_to_left = geom_center_pos[0]*0.01
    dist_to_top = geom_center_pos[1]*0.01
    dist_to_right = (mosaic.shape[1] - geom_center_pos[0])*0.01
    dist_to_bottom = (mosaic.shape[0] - geom_center_pos[1])*0.01
    print(dist_to_left,dist_to_bottom,dist_to_right,dist_to_top)
    print(dist_to_left+dist_to_right, dist_to_bottom+dist_to_top)

    bbox=(beam_model.properties['max_xccob']-dist_to_bottom,
          beam_model.properties['max_xccob']+dist_to_top,
          beam_model.properties['max_yccob']-dist_to_left,
          beam_model.properties['max_yccob']+dist_to_right)
    print(np.shape(mosaic))
    print(-bbox[0]+bbox[1])
    print(-bbox[2]+bbox[3])
    
    return bbox

def plot_results(qe, model, mosaic, data_ccdid, lct, pos, delta_x, delta_y, filename):

    xccob = float(pos.split('_')[0]) # CCS
    yccob = float(pos.split('_')[1]) # CCS
    ccdid, x_from_ccob_pos, y_from_ccob_pos = lct.focalMmToCcdPixel(yccob, xccob) # FP coord
    geom_center_pos = (np.shape(mosaic)[1]-(x_from_ccob_pos+delta_x),y_from_ccob_pos+delta_y) # EOtest coord

    fig, axes = plt.subplots(ncols=3, nrows=1,  figsize=(12, 5))
    vmin = np.median(mosaic.flatten())*0.95
    vmax = np.median(mosaic.flatten())*1.02
    im0 = axes[0].imshow(mosaic, vmin=vmin, vmax=vmax)
    axes[0].plot([np.shape(mosaic)[1]-x_from_ccob_pos],[y_from_ccob_pos], marker='+',color='red')
    axes[0].plot([geom_center_pos[0]],[geom_center_pos[1]], marker='o',color='yellow')
    im1 = axes[1].imshow(model)

    qe_sm = scipy.ndimage.filters.gaussian_filter(qe, 10, mode='constant')
    vmin = np.median(qe_sm.flatten())*0.99
    vmax = np.median(qe_sm.flatten())*1.01
#    vmin = np.max(mosaic.flatten())*0.99
#    vmax = np.max(mosaic.flatten())*1.01
#    vmin = 0.99
#    vmax = 1.01
    im2 = axes[2].imshow(qe_sm, vmin=vmin, vmax=vmax)

    axes[0].set_title(data_ccdid+', '+pos)
    axes[1].set_title('Beam model')
    axes[2].set_title('Flat = data/model')

    fig.colorbar(im0, ax=axes[0],fraction=0.046, pad=0.04)
    fig.colorbar(im1, ax=axes[1],fraction=0.046, pad=0.04)
    fig.colorbar(im2, ax=axes[2],fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close()



def make_fits(QE_map, amp_coord, outfile, template_file):
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

    for amp_pos in amp_coord.keys():

        amp = amp_coord[amp_pos]['amp']
        xmin = amp_coord[amp_pos]['xmin']
        xmax = amp_coord[amp_pos]['xmax']
        ymin = amp_coord[amp_pos]['ymin']
        ymax = amp_coord[amp_pos]['ymax']
        flipx = amp_coord[amp_pos]['flipx']
        flipy = amp_coord[amp_pos]['flipy']
        detsec = amp_coord[amp_pos]['detsec']
        datasec = amp_coord[amp_pos]['datasec']

        foo = QE_map[ymin:ymax,xmin:xmax]

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

    u.writeFits_from_dict(amp_dict_w_overscan, outfile, template_file)#,bitpix=-32)
    
    
if __name__ == '__main__':

    camera = camMapper._makeCamera()
    lct = LsstCameraTransforms(camera)

    
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
    path_to_beam = '/home/combet/tmp_9rafts/60x60/'
    outdir = '/home/combet/tmp_9rafts/QE_results/'

    led = 'red'
    
    model_ccdid = 'R22_S11'
    beam_model = load_beam_model(path_to_beam, led_name=led, ref_amp=5, ccdid=model_ccdid)
    delta_x, delta_y = compute_offsets(beam_model, lct, ccdid=model_ccdid, ref_pix_x=2304, ref_pixy=3003)
    
    raft_list = ['R22', 'R30', 'R01']
    sensor_list = ['S00','S01','S02','S10','S11','S12','S20','S21','S22']
    ccdid_list = sorted([raft+'_'+sensor for sensor in sensor_list for raft in raft_list])
    
    for data_ccdid in ccdid_list:
        gainfile = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_stage-test/LCA-10134_Cryostat/LCA-10134_Cryostat-0001/6801D/fe55_analysis_BOT/v0/47706/'+data_ccdid+'_6801D_eotest_results.fits'
            
        raft = data_ccdid.split('_')[0]
        data = qe_data.CcobQeData(data_ccdid, led, path_to_data[raft], gainfile, biasfile=None)
        data.find_dir()

        for pos in data.pos_list:
            data.make_avg_mosaic_at_pos(pos, '/home/combet/tmp_9rafts/')
            mosaic = data.data[pos]['mosaic']
            bbox = define_model_bbox(beam_model, mosaic, data.lct, pos, delta_x, delta_y)

            xarr = np.linspace(bbox[0],bbox[1],mosaic.shape[0])
            yarr = np.linspace(bbox[2],bbox[3],mosaic.shape[1])
            tmp = beam_model.beam_image['f_interp'](xarr, yarr) # interp works in cam coordinates
            tmp = np.flip(np.flip(tmp,axis=0),axis=1) # invert the model in x and in y 
            model_eotestDef = np.flip(tmp.T, axis=1) # beam model followinf EOTest convention
            model_normalised = model_eotestDef/np.max(model_eotestDef.flatten())
            qe = mosaic/model_normalised
      
            outfile = outdir+'fits/QE_'+data_ccdid+'_'+led+'_'+pos+'.fits'
            make_fits(qe, outfile, data.data[pos]['amp_coord'], data.template_file)

            figfile = outdir+'fits/QE_'+data_ccdid+'_'+led+'_'+pos+'.png'
            plot_results(qe, model, mosaic, data_ccdid, data.lct, pos, figfile)

