"""
All that is needed to reconstruct the CCOB beam and generate the beam model
used to produce the synthetic flat fields.
"""

import os
import sys
import glob
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.ndimage.filters import median_filter, gaussian_filter
from astropy.io import fits
import lsst.eotest.sensor as sensorTest
import ccob_utils as u

class CcobBeam:
    """Object that contains all relevant data and information for the beam model.
    The latter is obtained by scanning a bunch of reference pixels with the CCOB-WB projector.

    Attributes
    ----------
    config : dict
        Configuration for the beam reconstruction, containing paths to data, led, etc.
    properties : dict
        Contains properties that identify how the beam object was created, e.g., from which
        sensor, amplifier, reference pixels. Also contains location of the beam maximum.
    beam_image : dict
        Contains the information required to evaluate the intensity of the beam at any position.
    raw_data : dict
        Contains the raw data () at the location of the reference pixels throughout the scan.
        This records position of the CCOB and value of the control photodiode.
    """

    def __init__(self, config):
        self.config = config
        self.properties = {}
        self.beam_image = {}
        self.raw_data = {}

        
    def read_multibunch(self, config, dirlist=None, silent=False):
        """ Reads the data from a bunch of reference pixels after a CCOB scan and fills in self.raw_data 
        and self.properties
        
        Parameters
        ----------
        config : dict
            Contains all required information to reconstruct the beam
        dirlist : list
            List of directories containing the data
        silent : boolean
            If True, track the progress of the beam reconstruction
        """
        
        self.properties["ref_raft"] = ref_raft = 'R22' if ('ref_raft' not in config) else config['ref_raft']
        self.properties["ref_slot"] = ref_slot = 'S11' if ('ref_slot' not in config) else config['ref_slot']
        self.properties["ref_amp"] = ref_amps = [5] if ('ref_amps' not in config) else config['ref_amps']
        self.properties["ref_pix_x"] = ref_pix_x = 1000 if ('ref_pix_x' not in config) else int(config['ref_pix_x'])
        self.properties["ref_pix_y"] = ref_pix_y = 256 if ('ref_pix_y' not in config) else int(config['ref_pix_y'])
        self.properties["npix_for_avg"] = npix_for_avg = 30 if ('npix_for_avg' not in config) else int(config['npix_for_avg'])
        biasfile = None if ('biasfile' not in config) else config['biasfile']
        outdir = './' if ('tmpdir' not in config) else config['tmpdir']
               
        recons = {}
        led = self.config['led_name']

        dirlist = sorted(dirlist)
        if 'xarr' not in self.raw_data: self.raw_data['xarr']=[]
        if 'yarr' not in self.raw_data: self.raw_data['yarr']=[]
        if 'val' not in self.raw_data: self.raw_data['val']={}
        if 'pd_value' not in self.raw_data: self.raw_data['pd_value']=[]

        for i,d in enumerate(dirlist):
            dd = os.path.basename(d)
            xcurr = float(dd.split('_')[2])
            ycurr = float(dd.split('_')[3])
            l = list(zip(self.raw_data['xarr'], self.raw_data['yarr']))
#            print(i)

            if (xcurr,ycurr) in l:
                print('here')
                continue

            self.raw_data['xarr'].append(xcurr)
            self.raw_data['yarr'].append(ycurr)

            f = glob.glob(os.path.join(d, "*"+ref_raft+"*"+ref_slot+'*'))
            bb = fits.open(f[0])
            xh = np.round(bb[0].header['BOTXCAM'],1)
            yh = np.round(bb[0].header['BOTYCAM'],1)
#            self.raw_data['xarr'].append(np.round(bb[0].header['BOTX'],1))
#            self.raw_data['yarr'].append(np.round(bb[0].header['BOTY'],1))

            self.raw_data['pd_value'].append(bb[0].header['CCOBADC'])
            bb.close()
            ccd_dict = sensorTest.MaskedCCD(f[0], bias_frame=biasfile)
            for amp in ref_amps:
                image = ccd_dict.unbiased_and_trimmed_image(amp)
                arr = image.getArrays()
                val = np.mean((arr[0][int(ref_pix_x-npix_for_avg/2):int(ref_pix_x+npix_for_avg/2),\
                                                int(ref_pix_y-npix_for_avg/2):int(ref_pix_y+npix_for_avg/2)]))
                if amp in self.raw_data['val']:
                    self.raw_data['val'][amp].append(val)
                else:
                    self.raw_data['val'][amp]=[val]
            if not silent:
                #print(os.path.basename(f[0]))
                print(amp, self.raw_data['xarr'][-1], xh, self.raw_data['yarr'][-1], yh, val)
                #print(self.raw_data['val'])

        # Reordering the data
        x = self.raw_data['xarr']
        y = self.raw_data['yarr']
        val = self.raw_data['val']
        pd = self.raw_data['pd_value']
        newx = [x for x,_,_,_ in sorted(zip(x,y,val[ref_amps[0]],pd))]
        newy = [y for _,y,_,_ in sorted(zip(x,y,val[ref_amps[0]],pd))]
        newpd = [pd for _,_,_,pd in sorted(zip(x,y,val[ref_amps[0]],pd))]
        newval = {}
        for amp in ref_amps:
            newval[amp] = [val for _,_,val,_ in sorted(zip(x,y,val[amp],pd))]

        self.raw_data['xarr'] = newx
        self.raw_data['yarr'] = newy
        self.raw_data['val'] = newval
        self.raw_data['pd_value'] = newpd

        if outdir is not None:
            for amp in ref_amps:
                outf =os.path.join(outdir,'beam_raw_'+led+'_'+ref_raft+'_'+ref_slot+'_'
                                   +str(amp)+'_'+str(ref_pix_x)+'_'+str(ref_pix_y)+'.txt')
                np.savetxt(outf, np.array([newx, newy, newval[amp], newpd]).T, delimiter='   ')
 
        
  

    def interp_beam_BOT(self, xrange=None, yrange=None, step=1, pd_corr=False, amp=1, use_filt = False):

        """ Given the raw data, creates the corresponding beam model from cubic 
        spline interpolation, using the raw data as nodes.
        
        Parameters
        ----------
        xrange : tuple
            (xmin,xmax) for the interpolation
        yrange : tuple
            (ymin,ymax) for the interpolation
        step : int
            Defines the nodes for the interpolation. If step = 1, all positions in the raw data are used
            for the interpolation. step == 2, every second raw data is used, etc.
        pd_corr : boolean
            If True, interpolation is performed on the values corrected by the control photodiode value.
        amp : int
            Define the amp from which the data are to be considered for the interpolation
        use_filt : boolean
            If True, the data is first smoothed before being interpolated
            
        Returns
        -------
        The interpolation function is returned as a new attribute of the beam object, 
        self.beam_image['f_interp'] 
        """

        self.properties['analysis_amp'] = amp

        if not pd_corr:
            val_arr = np.array(self.raw_data['val'][amp])
        else:
            val_arr = np.array(self.raw_data['val'][amp])/np.array(self.raw_data['pd_value'])


        raw_image = np.reshape(val_arr,
                               (len(np.unique(self.raw_data['yarr'])), # number of lines
                                len(np.unique(self.raw_data['xarr'])))) # number of columns
#        self.raw_data['filtered'] = median_filter(raw_image, size=(5,5), mode = 'nearest').flatten()
        self.raw_data['filtered'] = gaussian_filter(raw_image, sigma=(1,1)).flatten()

        
        if use_filt: val_arr = np.array(self.raw_data['filtered'])

        if xrange is None:
            xrange = (min(self.raw_data['xarr']),max(self.raw_data['xarr']))

        if yrange is None:
            yrange= (min(self.raw_data['yarr']),max(self.raw_data['yarr']))

        xmin,xmax = xrange
        ymin,ymax = yrange

        filtx = ((np.array(self.raw_data['xarr']) >= xmin)*(np.array(self.raw_data['xarr']) <= xmax)) 

        tmp_x = np.array(self.raw_data['xarr'])[filtx]
        good_x = np.unique(tmp_x)[::step]
        filtx2 = [e in good_x for e in tmp_x]
        tmp_x = tmp_x[filtx2]

        tmp_y = np.array(self.raw_data['yarr'])[filtx][filtx2]
        tmp_val = val_arr[filtx][filtx2]

        filty = ((tmp_y >= ymin)*(tmp_y <= ymax))

        nodes = {}
        nodes['xarr'] = tmp_x[filty][::step]
        nodes['yarr'] = tmp_y[filty][::step]
        nodes['val'] = tmp_val[filty][::step]

        self.beam_image['nodes'] = nodes
        raw_image = np.reshape(nodes['val'],(len(np.unique(nodes['yarr'])),len(np.unique(nodes['xarr']))))
        f_interp = interpolate.interp2d(np.unique(nodes['xarr']), 
                                        np.unique(nodes['yarr']), 
                                        raw_image.T,
                                        kind='cubic')

        self.beam_image['f_interp'] = f_interp


    def make_image_BOT(self, ncols=300, nrows=300):

        """ Given the interpolation function, create a 2D array of the beam image.

        Parameters
        ----------
        ncols : int
            Number of columns of the beam image
        nrows : int
            Number of rows of the beam image
        
        Returns
        -------
        Beam model image as a 2D array
        """

        self.properties['ncols'] = ncols
        self.properties['nrows'] = nrows
        
#        extent = [min(self.beam_image['nodes']['xarr']),
#                  max(self.beam_image['nodes']['xarr']),
#                  min(self.beam_image['nodes']['yarr']), 
#                  max(self.beam_image['nodes']['yarr'])]

        extent = [min(self.raw_data['xarr']),
                  max(self.raw_data['xarr']),
                  min(self.raw_data['yarr']), 
                  max(self.raw_data['yarr'])]


        xarr = np.linspace(extent[0],extent[1],ncols)
        yarr = np.linspace(extent[2],extent[3],nrows)
        self.beam_image['xarr'] = xarr 
        self.beam_image['yarr'] = yarr
        self.beam_image['beam'] = self.beam_image['f_interp'](xarr, yarr)
        return self.beam_image['beam']
 

    def find_max(self):
        """ Finds the location of the beam maximum and save it in self.properties"""

        yarg,xarg = np.unravel_index(self.beam_image['beam'].argmax(), self.beam_image['beam'].shape)
        self.properties["max_xccob"] = self.beam_image['xarr'][xarg]
        self.properties["max_yccob"] = self.beam_image['yarr'][yarg]
        self.properties["max_xarg"] = xarg
        self.properties["max_yarg"] = yarg
        
   
    def find_max_from_avg(self):
        """ Same as find_max() but averages the position of the maximum on each row
        and column to define the beam maximum. This is more stable that using find_max()."""

        im_sm = gaussian_filter(self.beam_image['beam'], 5, mode='constant')

        xarg = np.mean([np.argmax(im_sm[i]) for i in np.arange(np.shape(im_sm)[0])])
        yarg = np.mean([np.argmax(im_sm[:,i]) for i in np.arange(np.shape(im_sm)[1])])
        
        
        self.properties["max_xccob"] = self.beam_image['xarr'][int(np.round(xarg))]
        self.properties["max_yccob"] = self.beam_image['yarr'][int(np.round(yarg))]
        self.properties["max_xarg"] = int(np.round(xarg))
        self.properties["max_yarg"] = int(np.round(yarg))
               
    def plot_BOT(self, aspect=None, outfile=None):        
        """
        Plots the beam and the location of its maximum if the information is available.
        """
#        extent = [min(self.beam_image['nodes']['yarr']),
#                  max(self.beam_image['nodes']['yarr']),
#                  min(self.beam_image['nodes']['xarr']), 
#                  max(self.beam_image['nodes']['xarr'])]
        
        colsize = (max(self.beam_image['nodes']['xarr'])-min(self.beam_image['nodes']['xarr']))/(self.properties['ncols'])
        rowsize = (max(self.beam_image['nodes']['yarr'])-min(self.beam_image['nodes']['yarr']))/(self.properties['nrows'])
        extent = [max(self.beam_image['nodes']['xarr'])+colsize/2,
                  min(self.beam_image['nodes']['xarr'])-colsize/2,
                  min(self.beam_image['nodes']['yarr'])-rowsize/2, 
                  max(self.beam_image['nodes']['yarr'])+rowsize/2]

        im = np.flip(self.beam_image['beam'],axis=1)

        plt.imshow(im/np.max(im.flatten()), extent=extent, aspect=aspect, origin='lower', vmin=0.7, vmax=1)
        plt.colorbar()
#        plt.scatter(self.beam_image['nodes']['xarr'],self.beam_image['nodes']['yarr'], marker='+', color='blue')
        if 'max_xccob' in self.properties:
            plt.plot([self.properties['max_xccob']],[self.properties['max_yccob']], marker='x', color='red', markersize='6')
        if outfile is None:
            plt.show()
        else:
            plt.savefig(outfile)

    
    def save(self, filename):
        """
        Saves the beam object as a pickle file for use at a later time.
        
        Parameters
        ----------        
        filename : str
            Name of the file into which to save the pickled beam object.
        """
        
        with open(filename, 'wb') as f:  # Overwrites any existing file.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)                

            
def main():
    """
    First, finds the list of files relevant to a CCOB scan, given the paths available in the configuration file. 
    Then,generates a CcobBeam object containing the raw data, the interpolated model of the beam, the 2D beam array,
    and the location of the beam maximum.
    """
    configfile = sys.argv[1]
    config = u.load_ccob_config(configfile)
    dirlist = []
    for d in config['rundir']:
        dirlist += glob.glob(config['rootdir']+d+'ccob_'+config['led_name']+'*')
    
    assert (len(dirlist)/config['scan_size']).is_integer, f'{len(dirlsit)} = Wrong number of scan locations'
        
    print(f'Scan over {len(dirlist)} locations')
    
    b = CcobBeam(config)
    for i in np.arange(config['scan_size']): 
        # that's just to allow saving the data at every line of the scan in case
        # the job is interrupted.
        print(f'Loading row {i}')
        start = i*config['scan_size']
        end = (i+1)*config['scan_size']
        b.read_multibunch(config, dirlist=dirlist[start:end], silent=True)
        b.save(os.path.join(config['tmpdir'],
                            'beam_object_'+config['ref_raft']+'_'+config['ref_slot']+'_'+config['led_name']+'.pkl'))   

    b.interp_beam_BOT(amp=config['ref_amp'], pd_corr=True)
    im = b.make_image_BOT()
    b.find_max_from_avg()
    b.save(os.path.join(config['tmpdir'],
                        'beam_object_'+config['ref_raft']+'_'+config['ref_slot']+'_'+config['led_name']+'.pkl'))   
      
if __name__ == '__main__':
    main()
 