import os
import lsst.eotest.sensor as sensorTest
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.ndimage.filters import median_filter, gaussian_filter
from numpy import unravel_index
import ccob_utils as u
import pickle 
import glob
import pdb
from astropy.io import fits as fits

class CcobBeam:
    
    def __init__(self, config):
        self.config = config
        self.properties = {}
        self.beam_image = {}
        self.raw_data = {}

    def read_multibunch(self, ref_raft='R22', ref_slot='S11', ref_amps=np.arange(1,17), ref_pix_x=1000,
               ref_pix_y=256, npix_for_avg=30, biasfile = None, dirlist=None, outdir=None):
        
        self.properties["ref_raft"] = ref_raft
        self.properties["ref_slot"] = ref_slot
        self.properties["ref_amp"] = ref_amps
        self.properties["ref_pix_x"] = int(ref_pix_x)
        self.properties["ref_pix_y"] = int(ref_pix_y)
        self.properties["npix_for_avg"] = int(npix_for_avg)
        
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
            ccd_dict = sensorTest.MaskedCCD(f[0])
            for amp in ref_amps:
                image = ccd_dict.unbiased_and_trimmed_image(amp)
                arr = image.getArrays()
                val = np.mean((arr[0][int(ref_pix_x-npix_for_avg/2):int(ref_pix_x+npix_for_avg/2),\
                                                int(ref_pix_y-npix_for_avg/2):int(ref_pix_y+npix_for_avg/2)]))
                if amp in self.raw_data['val']:
                    self.raw_data['val'][amp].append(val)
                else :
                    self.raw_data['val'][amp]=[val]
            #print(os.path.basename(f[0]))
            print(amp, self.raw_data['xarr'][-1], xh, self.raw_data['yarr'][-1], yh, val)
                #print(self.raw_data['val'])

        # Reordering the data
        x = self.raw_data['xarr']
        y = self.raw_data['yarr']
        val = self.raw_data['val']
        pd = self.raw_data['pd_value']
        newx = [x for x,_,_,_ in sorted(zip(x,y,val[1],pd))]
        newy = [y for _,y,_,_ in sorted(zip(x,y,val[1],pd))]
        newpd = [pd for _,_,_,pd in sorted(zip(x,y,val[1],pd))]
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
 
    def make_image_BOT(self, ncols=300, nrows=300):
        
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
        yarg,xarg = unravel_index(self.beam_image['beam'].argmax(), self.beam_image['beam'].shape)
        self.properties["max_xccob"] = self.beam_image['xarr'][xarg]
        self.properties["max_yccob"] = self.beam_image['yarr'][yarg]
        self.properties["max_xarg"] = xarg
        self.properties["max_yarg"] = yarg
#        return self.properties["max_xccob"], self.propperties["self.max_yccob"], xarg, yarg
 
    def find_max_from_avg(self):
        im_sm = gaussian_filter(self.beam_image['beam'], 5, mode='constant')

        xarg = np.mean([np.argmax(im_sm[i]) for i in np.arange(np.shape(im_sm)[0])])
        yarg = np.mean([np.argmax(im_sm[:,i]) for i in np.arange(np.shape(im_sm)[1])])
        
        
        self.properties["max_xccob"] = self.beam_image['xarr'][int(np.round(xarg))]
        self.properties["max_yccob"] = self.beam_image['yarr'][int(np.round(yarg))]
        self.properties["max_xarg"] = int(np.round(xarg))
        self.properties["max_yarg"] = int(np.round(yarg))
               
    def plot_BOT(self, aspect=None):        
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

        plt.imshow(im, extent=extent, aspect=aspect, origin='lower')
        plt.colorbar()
#        plt.scatter(self.beam_image['nodes']['xarr'],self.beam_image['nodes']['yarr'], marker='+', color='blue')
        if 'max_xccob' in self.properties:
            plt.plot([self.properties['max_xccob']],[self.properties['max_yccob']], marker='x', color='red', markersize='6')
        plt.show()

    
    def save(self, filename):
        with open(filename, 'wb') as f:  # Overwrites any existing file.
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
            
            
#     def read_in_beamdata_BOT(self, ref_raft='R10', ref_slot='S11', ref_amp=4, ref_pix_x=1000,
#                ref_pix_y=256, npix_for_avg=30, dirlist=None, outfile=None):
        
#         self.properties["ref_raft"] = ref_raft
#         self.properties["ref_slot"] = ref_slot
#         self.properties["ref_amp"] = ref_amp
#         self.properties["ref_pix_x"] = int(ref_pix_x)
#         self.properties["ref_pix_y"] = int(ref_pix_y)
#         self.properties["npix_for_avg"] = int(npix_for_avg)
        
#         recons = {}
#         led = self.config['led_name']
#         print(led)

#         dirlist = sorted(dirlist)
#         if 'xarr' not in self.raw_data: self.raw_data['xarr']=[]
#         if 'yarr' not in self.raw_data: self.raw_data['yarr']=[]
#         if 'val' not in self.raw_data: self.raw_data['val']=[]
#         if 'pd_value' not in self.raw_data: self.raw_data['pd_value']=[]
        
#         for i,d in enumerate(dirlist):
#             dd = os.path.basename(d)
# #            print(d, dd)
#             self.raw_data['xarr'].append(float(dd.split('_')[2]))
#             self.raw_data['yarr'].append(float(dd.split('_')[3]))

#             f = glob.glob(os.path.join(d, "*"+ref_raft+"*"+ref_slot+'*'))
#             bb = fits.open(f[0])
#             xh = np.round(bb[0].header['BOTXCAM'],1)
#             yh = np.round(bb[0].header['BOTYCAM'],1)
#             self.raw_data['pd_value'].append(bb[0].header['CCOBADC'])
#             bb.close()
#             ccd_dict = sensorTest.MaskedCCD(f[0])
#             image = ccd_dict.unbiased_and_trimmed_image(ref_amp)
#             arr = image.getArrays()
#             val = np.mean((arr[0][int(ref_pix_x-npix_for_avg/2):int(ref_pix_x+npix_for_avg/2),\
#                                                 int(ref_pix_y-npix_for_avg/2):int(ref_pix_y+npix_for_avg/2)]))
#             self.raw_data['val'].append(val)
# #            print(os.path.basename(f[0]))
#             print(self.raw_data['xarr'][-1], xh, self.raw_data['yarr'][-1], yh, val)

#         x = self.raw_data['xarr']
#         y = self.raw_data['yarr']
#         val = self.raw_data['val']
#         pd = self.raw_data['pd_value']
        
#         newx = [x for x,_,_,_ in sorted(zip(x,y,val,pd))]
#         newy = [y for _,y,_,_ in sorted(zip(x,y,val,pd))]
#         newval = [val for _,_,val,_ in sorted(zip(x,y,val,pd))]
#         newpd = [pd for _,_,_,pd in sorted(zip(x,y,val,pd))]

#         self.raw_data['xarr'] = newx
#         self.raw_data['yarr'] = newy
#         self.raw_data['val'] = newval
#         self.raw_data['pd_value'] = newpd


#         if outfile is not None:
#             np.savetxt(outfile, np.array([newx, newy, newval, newpd]).T, delimiter='   ')
     
    
#     def recons(self, ref_slot='11', ref_amp=4, ref_pix_x=1000,
#                ref_pix_y=256, npix_for_avg=30):

#         self.properties["ref_slot"] = ref_slot
#         self.properties["ref_amp"] = ref_amp
#         self.properties["ref_pix_x"] = int(ref_pix_x)
#         self.properties["ref_pix_y"] = int(ref_pix_y)
#         self.properties["npix_for_avg"] = int(npix_for_avg)
        
#         recons = {}
#         led = self.config['led_name']
#         print(led)

#         flist = sorted(u.find_files(self.config, slot=ref_slot))
#         nodes = {}
#         nodes['xarr'] = []
#         nodes['yarr'] = []
#         nodes['val'] = []

#         for i,f in enumerate(sorted(flist)):
#             nodes['xarr'].append(float(os.path.basename(f).split('_')[3].split('x')[1]))
#             nodes['yarr'].append(float(os.path.basename(f).split('_')[4].split('y')[1]))
#             ccd_dict = sensorTest.MaskedCCD(f)
#             image = ccd_dict.unbiased_and_trimmed_image(ref_amp)
#             arr = image.getArrays()
#             nodes['val'].append(np.mean((arr[0][int(ref_pix_x-npix_for_avg/2):int(ref_pix_x+npix_for_avg/2),\
#                                                 int(ref_pix_y-npix_for_avg/2):int(ref_pix_y+npix_for_avg/2)])))
#         f_interp = interpolate.interp2d(np.unique(nodes['xarr']), np.unique(nodes['yarr']), nodes['val'], kind='cubic')
#         self.beam_image['nodes'] = nodes
#         self.beam_image['f_interp'] = f_interp
# #        return self.beam_image
                
                     
        
if __name__ == '__main__':
    config_file = '/home/combet/ccob-wb/ccob_beam_recons_config.yaml'
    beam = b.CcobBeam(config_file=config_file)
    beam.recons()
    beam.make_image()
    beam.find_max()
    beam.save(os.path.join(beam.config['tmp_dir'],beam.config['led_name']+'_beam.pkl'))
