import os
import glob
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import lsst.eotest.raft as raft
import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy import interpolate
import ccob_utils as u


led_names = ['uv', 'blue', 'red', 'nm750','nm850','nm960']


def beam_recons(led_names, ref_slot='11', ref_amp=4, ref_pix_x=1000,
               ref_pix_y=256, npix_for_avg=30, npix_beam_image=300):


    mean_bias_file = ref_slot+'_mean_bias_image_RTM-006.fits'
    recons = {}

    for led in led_names:
        config = u.load_ccob_config('/home/combet/ccob-wb/ccob_config_RTM-006.yaml')
        config['led_name'] = led
        config['path'] = os.path.join(config['path'],config['led_name'])

        flist = sorted(u.find_files(config, slot=ref_slot))
        ccd_dict = sensorTest.MaskedCCD(flist[0], bias_frame=os.path.join(config['tmp_dir'],mean_bias_file))
        eotest_results_file = os.path.join(config['eo_data_path'], '{}_eotest_results.fits'.format(ccd_dict.md('LSST_NUM')))
        gains_dict = u.gains(eotest_results_file)

        nodes = {}
        nodes['xarr'] = []
        nodes['yarr'] = []
        nodes['val'] = []
        for i,f in enumerate(sorted(flist)):
            nodes['xarr'].append(float(os.path.basename(f).split('_')[3].split('x')[1]))
            nodes['yarr'].append(float(os.path.basename(f).split('_')[4].split('y')[1]))
            ccd_dict = sensorTest.MaskedCCD(f, bias_frame=os.path.join(config['tmp_dir'],mean_bias_file))
            image = ccd_dict.unbiased_and_trimmed_image(ref_amp)
            image *= gains_dict[ref_amp]
            arr = image.getArrays()
            nodes['val'].append(np.mean((arr[0][ref_pix_x-npix_for_avg/2:ref_pix_x+npix_for_avg/2,
                                                ref_pix_y-npix_for_avg/2:ref_pix_y+npix_for_avg/2])))

        f_interp = interpolate.interp2d(np.unique(nodes['xarr']), np.unique(nodes['yarr']), nodes['val'], kind='cubic')
        xarr2 = np.linspace(min(nodes['xarr']),max(nodes['xarr']),npix_beam_image)
        yarr2 = np.linspace(min(nodes['yarr']),max(nodes['yarr']),npix_beam_image)
        tmp = f_inerp(xarr2, yarr2)
        recons[led] = tmp/max(tmp.flatten())

    return recons
    
    
    
extent = [min(nodes['xarr']), max(nodes['xarr']),min(nodes['yarr']), max(nodes['yarr'])]
fig, axes = plt.subplots(ncols=2, nrows=1)
axes[0].imshow(recons[led], extent=extent)
axes[0].colorbar()
axes[0].scatter(nodes['xarr'],nodes['yarr'], marker='+', color='blue')
            
            