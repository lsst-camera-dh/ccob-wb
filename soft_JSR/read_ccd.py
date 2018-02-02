import sys
sys.path.append('/sps/lsst/CCOB/eotest/python/')
sys.path.append('/sps/lsst/CCOB/eotest/python/lsst/eotest/sensor/')

from MaskedCCD import MaskedCCD


image_file = "/sps/lsst/DataBE/lpnws5203/data/frames//20170616/10_CCD1_20170616091825.fz"

ccd = MaskedCCD(image_file)


