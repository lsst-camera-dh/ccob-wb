import ccob_qe_analysis as qe
import ccob_utils as u
import os
import glob

if __name__ == '__main__':

 
    slot_names=['11']
#    slot_names = ['00','01','02','10','11','12','20','21','22']
#    led_names = ['nm960','850nm', '750nm', 'red', 'blue', 'uv']
    led_names = ['uv','nm850', 'nm750', 'blue']

    config_file_beam = 'ccob_beam_recons_config.yaml'
    config_file_data = 'ccob_config_RTM-006.yaml'

    for led in led_names:
        for slot in slot_names:
            print("=============", led, slot)
            QE = qe.CcobQE(config_file_beam, config_file_data)
            QE.make_ccob_beam(led_name=led, ref_slot=slot)
            ccd = QE.load_ccd(led_name=led)
            QE.compute_QE()

            config = u.load_ccob_config(config_file_data)
            flist = glob.glob(os.path.join(os.path.join(config['path'],led),slot+'*'))
            template_file = flist[0]
            print(template_file)
            outfile = os.path.join(config['tmp_dir'], slot+'_CCOB_QE_'+ led + '.fits')
            QE.make_fits(outfile, template_file)
        #    QE.plot_QE()