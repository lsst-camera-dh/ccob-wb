import os
import glob
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

class CcobQeData: 
    """
    A class to compute and save the relative QE measurements perfomed using CCOB data
    
    Attributes
    ----------
        ccdid: string
            Path to the file containing information necessary to the CCOB beam reconstruction
        led: string
            Path to the file containing information to read in the sensor data used to compute the CCOB QE
        path_to_data: CcobBeam object
            Contains all the properties of the CCOB beam reconstructed from the config_file_beam configuration file
        gainfile: 2D array
            Mosaic image of the CCOB-illuminated sensor to be used in the QE calculation
        biasfile: 2D array
            Mosaic image of the QE obtained from ccd/beam
    
    """
    def __init__(self, ccdid, led, path_to_data, gainfile, biasfile=None):
        self.ccdid = ccdid
        self.gains = u.gains(gainfile)
        self.biasfile= biasfile
        self.led = led
        self.path_to_data = path_to_data
        self.data = {}
        camera = camMapper._makeCamera()
        self.lct = LsstCameraTransforms(camera)

        
    def find_dir(self):
        allfiles = glob.glob(os.path.join(self.path_to_data,'*'+self.led+'*'))
        filenames=[os.path.basename(f) for f in allfiles]
        fnames_comp=[f.split('_') for f in filenames]
        x_ccob = np.array(fnames_comp).T[2].astype(np.float)
        y_ccob = np.array(fnames_comp).T[3].astype(np.float)

        tmp = [self.lct.focalMmToCcdPixel(t[0],t[1]) for t in zip(y_ccob,x_ccob)]
        ccdid_list=np.array(tmp).T[0]
        flag = (ccdid_list == self.ccdid)
        
        self.dir_list = np.sort(np.array(allfiles)[flag])

        dirnames=[os.path.basename(d) for d in self.dir_list]
        dnames_comp=[d.split('_') for d in dirnames]
        x_ccob = np.array(dnames_comp).T[2].astype(np.float)
        y_ccob = np.array(dnames_comp).T[3].astype(np.float)
        pos_str = [str(x_ccob[i])+'_'+str(y_ccob[i]) for i in np.arange(len(x_ccob))]
        self.pos_list = np.unique(pos_str)


    def make_avg_mosaic_at_pos(self, pos, outdir):
        """
        Create a mosaic image from a given sensor illuminated 
        by the CCOB. The list of fileThe path to the data is provided in the self.config_file_data file.

        Parameters
        ----------
            led_name: choice of the CCOB LED. Either one of ['nm960','nm850','nm750,'red','blue,'uv'] 
        """

        dirlist = [d for d in self.dir_list if pos in d]
        flist = np.reshape(np.array([glob.glob(os.path.join(d, '*'+self.ccdid+'*')) for d in dirlist]),(len(dirlist)))
        print('Averaging the following files:\n')
        print(flist)
        
        tmp_file = os.path.join(outdir,'tmp.fits')
        imutils.fits_mean_file(flist, os.path.join(outdir,tmp_file))
        
        self.template_file = flist[0]
        
        mosaic, amp_coord = u.make_ccd_2d_array(tmp_file, gains=self.gains, biasfile=self.biasfile)
        self.data[pos] = {'mosaic':mosaic, 
                          'amp_coord':amp_coord}
        
    def plot_mosaic(self, pos):

        mosaic = self.data[pos]['mosaic']
        xccob = float(pos.split('_')[0]) # CCS
        yccob = float(pos.split('_')[1]) # CCS

        ccdid, x_from_ccob_pos, y_from_ccob_pos = self.lct.focalMmToCcdPixel(yccob, xccob) # FP coord

        vmin = np.median(mosaic.flatten())*0.9
        vmax = np.median(mosaic.flatten())*1.02
        
        plt.imshow(mosaic, vmin=vmin, vmax=vmax)
        plt.plot([np.shape(mosaic)[1]-x_from_ccob_pos],[y_from_ccob_pos], marker='+',color='red')
        plt.colorbar()
        plt.show()

