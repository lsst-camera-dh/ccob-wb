# ccob-wb

This repository contains the analysis software for the Camera Calibration Optical Bench - Wide Beam projector (CCOB-WB). The CCOB-WB aims at producing a composite flat field of the LSST focal plane (without the optics) at a precision of a few per mil. 

The CCOB-WB consists in 6 LED, emitting light at  wavelength characteristic of the 6 LSST filters, mounted on a integrating sphere. The outgoing beam, 3-4 cm wide, can illuminate any part fo the focal plane. The measurement is performed in two steps:

- First, the beam must be reconstructed at the per mil level. This is done by choosing a bunch of reference pixels and scanning those pixels with the beam over a fine grid. This allows to define a beam model. 
- Second, each CCD of the focal plane is illuminated by the beam, in each wavelength. The corresponding image divided by the beam model computed at step 1, will be provide a synthetic flat field of that CCD.

Stability during the scan is monitored/ensured by a control photodiode mounted on the sphere.

## Requirements

`ccob-wb` requires Python version 3.6 or later and depends on:

- [astropy](https://www.astropy.org/) 
- [eotest](https://github.com/lsst-camera-dh/eotest)
- [matplotlib](https://matplotlib.org/)
- [numpy](http://www.numpy.org/)
- [scipy](http://www.scipy.org/)
- [dm-stack](https://pipelines.lsst.io/) (LSST analysis pipeline)

## Running the analysis

At the moment, the analysis is not automated and the various steps need to be run by hand. Also, data may be spread over several directories, depending on how the acquisition went.

- `make_raw_beam.ipynb` [needs cleanup] in the notebook folder, examplifies how to reconstruct the beam model from the scan over the bunch of reference pixels. This produces a pickle file containing the beam object. There are several attributes to the beam object, including: