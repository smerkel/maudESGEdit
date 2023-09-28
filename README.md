# maudESGEdit
Utility to fix data in ESG files before Rietveld refinement in MAUD

This programs allows for
- Removing rubbish data points, due to detector gaps, for instance,
- Removing and correction your background. 

This program is written in python and should work on most platforms. Download maudESGEdit.py and run it with python. It is the only file you need.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

The source code is available online at https://github.com/smerkel/maudESGEdit/. 

## Extract from the User Manual (see help -> user manual for more)

This programs allows for
 * Removing rubbish data points, due to detector gaps, for instance,
 * Removing and correction your background.


### Rubbish data points

Rubbish data points occur when you  have either a parasite signal or gaps in your data collection. MAUD can exclude 2theta ranges, but will fail if the parasite signals occur at inconsistent 2theta values. This utility allows to browse through all azimuth in an ESG file and remove these rubbish data points.

How to proceed to remove rubbish data points?
 * Click on the zoom icon before you start, to activate the zoom mode,
 * Zoom-in on the data points you wish to remove,
 * Click on *Remove data points* or hit *Ctr-d*.


### Background

First, if your background can be treated in MAUD, treat them in MAUD. It is always a million times better to fit your signal within a single software with clear hypothesis for the modeling of experimental data. This utility is for hopeless cases, inconsistent detector chips inducing jumps in the overall background level, for instance.

How to subtract background?
 * Right click to generate background points (background interpolation is linear),
 * Select *Subtract background* when you have enough points.

### Baseline intensity

MAUD will fail with negative or zero intensity points. It actually is a feature, but can really throw you off with some datas. Your refinement will diverge. To fix this, You can
 * Add a fixed intensity to the data at the current azimuth by selecting *Shift dataset...*,
 * Add a fixed intensity to the data at all azimuth by selecting *Shift all datasets...*,
 * Set the minimum intensity at all azimuth by selecting *Set value for minimum intensity...*. For each azimuth, the program will evaluate the current intensity minimum and add the necessary shift to bring the minimum intensity to the value decided by the user. This operation is typically performed at the end, once you have edited data for all azimuth.

### Navigation and file management

You can navigate between spectra using the *Previous* or *Next* buttons or by using the *left* and *right* keyboard keys.

When you are done, save your new data in a new ESG file.

### Masks

As you remove rubbish data points, we record 2theta ranges along with the corresponding azimuths. You can save these ranges in a file, with a *msk* extension to reuse them later.

If you want to remove the same data ranges as in a previous processing, use the *Mask -> Load and apply mask*menu item.

### Final note

Is this data manipulation? If you use this sotware to remove actual data, it is. If you use this software to clean up spectra (due to gaps in your detectors, for instance), it is not.

Good luck with your data!
