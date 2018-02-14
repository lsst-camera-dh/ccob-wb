import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import ccob_utils

def main():
    config_file = 'ccob_config.yaml'
    config = ccob_utils.load_ccob_config(config_file)
    slot_names = ['00','01','02','10','11','12','20','21','22']
    led_names = ['uv', 'blue','nm750', 'nm850', 'red', 'nm960']
    currents = ['0.012A','0.012A','0.012A','0.012A','0.012A','0.012A']
    exp_times = ['0.6s', '0.03s', '0.06s','0.12s','0.06s','0.18s']

    xpos = ['258', '300', '342']
    ypos = ['234', '192', '150']
    
    for slot in slot_names:
        ccob_utils.build_mean_bias_frame(config, slot)  

    for i,led in enumerate(led_names):
        config['led_name'] = led
        config['current'] = currents[i]
        config['exp_time'] = exp_times[i]
        for x in xpos:
            for y in ypos:
                print "Working on led = "+led+", x = "+x+", y = "+y
                config['xpos'] = x
                config['ypos'] = y
                im_raw, im_corr = ccob_utils.make_image(config, slot_names)  
                fig_raw = im_raw.plot(nsig=2, title = config['led_name'])
                fig_corr = im_corr.plot(nsig=2, title = config['led_name'])
        
                filename = 'raft_image_raw_ccob_X'+config['xpos']+'_Y'+config['ypos']+'_'+config['led_name']+'.png'
                fig_raw.savefig(os.path.join(config['tmp_dir'],filename))
                filename = 'raft_image_corr_ccob_X'+config['xpos']+'_Y'+config['ypos']+'_'+config['led_name']+'.png'
                fig_corr.savefig(os.path.join(config['tmp_dir'],filename))
    
# this means that if this script is executed, then 
# main() will be executed
if __name__ == '__main__':
    main()