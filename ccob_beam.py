import os
import glob
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import lsst.eotest.raft as raft
import matplotlib.pyplot as plt
import numpy as np
import yaml
import sys
from scipy import interpolate
from astropy.table import Table
from numpy import unravel_index
import ccob_utils as u

class CcobBeam:
    
    def __init__(self, config_file, npix_x=300, npix_y=300):
        self.config = u.load_ccob_config(config_file)
        self.npix_x = npix_x
        self.npix_y = npix_y
    

    def recons(self, ref_slot='11', ref_amp=4, ref_pix_x=1000,
               ref_pix_y=256, npix_for_avg=30):

        self.ref_slot = ref_slot
        self.ref_amp = ref_amp
        self.ref_pix_x = int(ref_pix_x)
        self.ref_pix_y = int(ref_pix_y)
        self.npix_for_avg = int(npix_for_avg)
        
        recons = {}
        led = self.config['led_name']
        self.config['path'] = os.path.join(self.config['path'],led)
        print(led)

        flist = sorted(u.find_files(self.config, slot=ref_slot))
        nodes = {}
        nodes['xarr'] = []
        nodes['yarr'] = []
        nodes['val'] = []
        for i,f in enumerate(sorted(flist)):
            nodes['xarr'].append(float(os.path.basename(f).split('_')[3].split('x')[1]))
            nodes['yarr'].append(float(os.path.basename(f).split('_')[4].split('y')[1]))
            ccd_dict = sensorTest.MaskedCCD(f)
            image = ccd_dict.unbiased_and_trimmed_image(ref_amp)
            arr = image.getArrays()
            nodes['val'].append(np.mean((arr[0][int(ref_pix_x-npix_for_avg/2):int(ref_pix_x+npix_for_avg/2),\
                                                int(ref_pix_y-npix_for_avg/2):int(ref_pix_y+npix_for_avg/2)])))

        f_interp = interpolate.interp2d(np.unique(nodes['xarr']), np.unique(nodes['yarr']), nodes['val'], kind='cubic')
        xarr2 = np.linspace(min(nodes['xarr']),max(nodes['xarr']),self.npix_x)
        yarr2 = np.linspace(min(nodes['yarr']),max(nodes['yarr']),self.npix_y)
        tmp = f_interp(xarr2, yarr2)
        self.beam_image={}
        self.beam_image['nodes'] = nodes
        self.beam_image['xarr'] = xarr2
        self.beam_image['yarr'] = yarr2
        self.beam_image['beam'] = tmp/max(tmp.flatten())
        self.beam_image['f_interp'] = f_interp
        return self.beam_image
                
        
    def plot(self):
        extent = [min(self.beam_image['nodes']['xarr']),
                  max(self.beam_image['nodes']['xarr']),
                  min(self.beam_image['nodes']['yarr']), 
                  max(self.beam_image['nodes']['yarr'])]
        plt.imshow(self.beam_image['beam'], extent=extent, origin='lower')
        plt.colorbar()
        plt.scatter(self.beam_image['nodes']['xarr'],self.beam_image['nodes']['yarr'], marker='+', color='blue')
        plt.show()
    
    def find_max(self):
        xarg, yarg = unravel_index(self.beam_image['beam'].argmax(), self.beam_image['beam'].shape)
        self.max_xccob = self.beam_image['xarr'][xarg]
        self.max_yccob = self.beam_image['yarr'][yarg]
        return self.max_xccob, self.max_yccob, xarg, yarg
    
    def write_to_file(self, filename):
        my_table=Table([self.beam_image['beam'],
               self.beam_image['xarr'],
               self.beam_image['yarr']])
        my_table.write(filename)
