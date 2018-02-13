import pyfits
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.ndimage as ndimage


def fill_seg_dict(d):
    seg_dict = {}
    for i in np.arange(1,17):
        datasec = [d[i].header['DATASEC'].strip('[]').split(',')[0].split(':')[0],
                  d[i].header['DATASEC'].strip('[]').split(',')[0].split(':')[1],
                  d[i].header['DATASEC'].strip('[]').split(',')[1].split(':')[0],
                  d[i].header['DATASEC'].strip('[]').split(',')[1].split(':')[1]]
        detsec = [d[i].header['DETSEC'].strip('[]').split(',')[0].split(':')[0],
                  d[i].header['DETSEC'].strip('[]').split(',')[0].split(':')[1],
                  d[i].header['DETSEC'].strip('[]').split(',')[1].split(':')[0],
                  d[i].header['DETSEC'].strip('[]').split(',')[1].split(':')[1]]
        detsize = [d[i].header['DETSIZE'].strip('[]').split(',')[0].split(':')[0],                      
                   d[i].header['DETSIZE'].strip('[]').split(',')[0].split(':')[1],
                  d[i].header['DETSIZE'].strip('[]').split(',')[1].split(':')[0],
                  d[i].header['DETSIZE'].strip('[]').split(',')[1].split(':')[1]]
     
        seg_dict[d[i].header['EXTNAME']] = {'datasec':list(map(int,datasec)),
                                           'detsec':list(map(int,detsec)),
                                           'detsize':list(map(int,detsize)),
                                           'data':d[i].data-d[i].header['AVGBIAS']}
    return seg_dict


def channel_intercal(seg_dict):
    
    seg_ref = list(seg_dict.keys())[0]
    xmin = seg_dict[seg_ref]['datasec'][0]
    xmax = seg_dict[seg_ref]['datasec'][1]
    ymin = seg_dict[seg_ref]['datasec'][2]
    ymax = seg_dict[seg_ref]['datasec'][3]

    filt=np.abs(seg_dict[seg_ref]['data'][ymin-1:ymax-1,xmin-1:xmax-1].flatten() - np.median(seg_dict[seg_ref]['data'][ymin-1:ymax-1,xmin-1:xmax-1]))<2*np.std(seg_dict[seg_ref]['data'][ymin-1:ymax-1,xmin-1:xmax-1])
    ref_channel_mean = np.median(seg_dict[seg_ref]['data'][ymin-1:ymax-1,xmin-1:xmax-1].flatten()[filt])
#    print(np.sum(filt),ref_channel_mean)
    norm=np.zeros(16)
    for i,seg in enumerate(list(seg_dict.keys())):
        xmin = seg_dict[seg]['datasec'][0]
        xmax = seg_dict[seg]['datasec'][1]
        ymin = seg_dict[seg]['datasec'][2]
        ymax = seg_dict[seg]['datasec'][3]
        
        filt=np.abs(seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1].flatten()- np.median(seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1]))<2*np.std(seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1])

        norm = ref_channel_mean / np.mean(seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1].flatten()[filt])    
        seg_dict[seg]['intercalib']=norm

def make_ccd_image(seg_dict, seg_dict_cal, intercal=False):
    img_tot = np.empty(shape=(seg_dict['Segment00']['detsize'][3],seg_dict['Segment00']['detsize'][1]))
    for seg in seg_dict.keys():
        xmin = seg_dict[seg]['datasec'][0]
        xmax = seg_dict[seg]['datasec'][1]
        ymin = seg_dict[seg]['datasec'][2]
        ymax = seg_dict[seg]['datasec'][3]
        xmin_detset = seg_dict[seg]['detsec'][0]
        xmax_detset = seg_dict[seg]['detsec'][1]
        ymin_detset = seg_dict[seg]['detsec'][2]
        ymax_detset = seg_dict[seg]['detsec'][3]
    
        if (xmin_detset < xmax_detset) & (ymin_detset < ymax_detset):
#            print('case1')
            if intercal:
                img_tot[ymin_detset-1:ymax_detset-1,xmin_detset-1:xmax_detset-1] =\
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1]*seg_dict_cal[seg]['intercalib']
            else:
                img_tot[ymin_detset-1:ymax_detset-1,xmin_detset-1:xmax_detset-1] = \
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1]

        if (xmin_detset < xmax_detset) & (ymin_detset > ymax_detset):
#            print('case2')
            if intercal:
                img_tot[ymax_detset-1:ymin_detset-1,xmin_detset-1:xmax_detset-1] =\
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1][::-1,:]*seg_dict_cal[seg]['intercalib']
            else:
                img_tot[ymax_detset-1:ymin_detset-1,xmin_detset-1:xmax_detset-1] =\
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1][::-1,:]

        if (xmin_detset > xmax_detset) & (ymin_detset < ymax_detset):
#            print('case3')
            if intercal:
                img_tot[ymin_detset-1:ymax_detset-1,xmax_detset-1:xmin_detset-1] =\
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1][:,::-1]*seg_dict_cal[seg]['intercalib']
            else:
                img_tot[ymin_detset-1:ymax_detset-1,xmax_detset-1:xmin_detset-1] =\
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1][:,::-1]

        if (xmin_detset > xmax_detset) & (ymin_detset > ymax_detset):
#            print('case4')
            if intercal:
                img_tot[ymax_detset-1:ymin_detset-1,xmax_detset-1:xmin_detset-1] =\
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1][::-1,::-1]*seg_dict_cal[seg]['intercalib']
            else:
                img_tot[ymax_detset-1:ymin_detset-1,xmax_detset-1:xmin_detset-1] =\
                  seg_dict[seg]['data'][ymin-1:ymax-1,xmin-1:xmax-1][::-1,::-1]

    return img_tot

def raw_rack_image(d, dcal, intercal=True):
    im=[]
    for ccd in sorted(d.keys()):
        print(ccd)
        im.append(make_ccd_image(d[ccd], dcal[ccd], intercal=intercal))
    return im

def filtered_rack_image(im):
    return ndimage.gaussian_filter(im, sigma=(9, 9), order=0)


def plot_rack_image(im):

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10,10))
    axes[0,2].imshow(np.flip(im[0],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[0,2].axis('off')
    img = axes[0,1].imshow(np.flip(im[1],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[0,1].axis('off')
    axes[0,0].imshow(np.flip(im[2],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[0,0].axis('off')
    axes[1,2].imshow(np.flip(im[3],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[1,2].axis('off')
    axes[1,1].imshow(np.flip(im[4],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[1,1].axis('off')
    axes[1,0].imshow(np.flip(im[5],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[1,0].axis('off')
    axes[2,2].imshow(np.flip(im[6],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[2,2].axis('off')
    axes[2,1].imshow(np.flip(im[7],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[2,1].axis('off')
    axes[2,0].imshow(np.flip(im[8],axis=1), vmin=10000, vmax=120000, origin='lower', cmap='plasma')
    axes[2,0].axis('off')
    fig.subplots_adjust(wspace=0, hspace=0)
    cbar_ax = fig.add_axes([0.95, 0.12, 0.03, 0.75])
    fig.colorbar(img, cax=cbar_ax)
    fig.show()




ledcurrent = '_0.01A'
exptime = '_0.1s'
led = 'red'

dirname = '/sps/lsst/data/ccombet/CCOB/data_slac_1711/171117/testccob/'
ccd_names = ['00','01','02','10','11','12','20','21','22']
ccd_pos1 = ['-42','-42','-42','0','0','0','42','42','42']
ccd_pos2 = ['42','0','-42','42','0','-42','42','0','-42']
pos = list(zip(ccd_pos1,ccd_pos2))
d_cal={}
for i,ccd in enumerate(ccd_names):
    filename_base = ccd + '_CCOB_' + led + ledcurrent + exptime
    l = glob.glob(dirname+'xy_'+pos[i][0]+'_'+pos[i][1]+'/'+filename_base+'*')
    myfile = l[0]
    dd = pyfits.open(myfile)
    print(ccd, myfile)
    d_cal[ccd]=fill_seg_dict(dd)
    channel_intercal(d_cal[ccd])




ledcurrent = '_0.01A'
exptime = '_0.1s'
led = 'red'

dirname = '/sps/lsst/data/ccombet/CCOB/data_slac_1711/171116/testccob/'
ccd_names=['00','01','02','10','11','12','20','21','22']
d={}
for ccd in ccd_names:
    filename_base = ccd + '_CCOB_' + led + ledcurrent + exptime
    l = glob.glob(dirname+filename_base+'*')
    myfile = l[0]
    dd = pyfits.open(myfile)
    print(ccd, myfile)
    d[ccd]=fill_seg_dict(dd)
    channel_intercal(d[ccd])





