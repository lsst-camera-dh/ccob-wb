import ccob_qe_analysis as qe
import ccob_utils as u
import os
import glob

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
