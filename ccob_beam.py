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
import pickle 

class CcobBeam:
    
    def __init__(self, config):
        self.config = config
        self.properties = {}
        self.beam_image={}


    def recons(self, ref_slot='11', ref_amp=4, ref_pix_x=1000,
               ref_pix_y=256, npix_for_avg=30):

        self.properties["ref_slot"] = ref_slot
        self.properties["ref_amp"] = ref_amp
        self.properties["ref_pix_x"] = int(ref_pix_x)
        self.properties["ref_pix_y"] = int(ref_pix_y)
        self.properties["npix_for_avg"] = int(npix_for_avg)
        
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
        self.beam_image['nodes'] = nodes
        self.beam_image['f_interp'] = f_interp
#        return self.beam_image
                
        
    def make_image(self, ncols=300, nrows=300):
        
        self.properties['ncols'] = ncols
        self.properties['nrows'] = nrows
        
        extent = [min(self.beam_image['nodes']['xarr']),
                  max(self.beam_image['nodes']['xarr']),
                  min(self.beam_image['nodes']['yarr']), 
                  max(self.beam_image['nodes']['yarr'])]
 
        xarr = np.linspace(extent[0],extent[1],ncols)
        yarr = np.linspace(extent[2],extent[3],nrows)
        self.beam_image['xarr'] = xarr
        self.beam_image['yarr'] = yarr
        self.beam_image['beam'] = self.beam_image['f_interp'](xarr, yarr)
        return self.beam_image['beam']
    

    
    def find_max(self):
        yarg,xarg = unravel_index(self.beam_image['beam'].argmax(), self.beam_image['beam'].shape)
        self.properties["max_xccob"] = self.beam_image['xarr'][xarg]
        self.properties["max_yccob"] = self.beam_image['yarr'][yarg]
        self.properties["max_xarg"] = xarg
        self.properties["max_yarg"] = yarg
#        return self.properties["max_xccob"], self.propperties["self.max_yccob"], xarg, yarg
    
    def plot(self):        
        extent = [min(self.beam_image['nodes']['xarr']),
                  max(self.beam_image['nodes']['xarr']),
                  min(self.beam_image['nodes']['yarr']), 
                  max(self.beam_image['nodes']['yarr'])]
        plt.imshow(self.beam_image['beam'], extent=extent, origin='lower', aspect='equal')
        plt.colorbar()
        plt.scatter(self.beam_image['nodes']['xarr'],self.beam_image['nodes']['yarr'], marker='+', color='blue')
        if 'max_xccob' in self.properties:
            plt.plot([self.properties['max_xccob']],[self.properties['max_yccob']], marker='x', color='red', markersize='6')
        plt.show()

    
    def save(self, filename):
        with open(filename, 'wb') as f:  # Overwrites any existing file.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
            
        
if __name__ == '__main__':
    config_file = '/home/combet/ccob-wb/ccob_beam_recons_config.yaml'
    beam = b.CcobBeam(config_file=config_file)
    beam.recons()
    beam.make_image()
    beam.find_max()
    beam.save(os.path.join(beam.config['tmp_dir'],beam.config['led_name']+'_beam.pkl'))
