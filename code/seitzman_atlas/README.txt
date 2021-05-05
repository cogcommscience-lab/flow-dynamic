There are two files called ROIs_300inVol_*_allInfo.txt, where * = MNI or Talairach. 
These files correspond to the set of 300 ROIs described in our manuscript found 
here:  https://www.biorxiv.org/content/early/2018/10/24/450452 .

There are 301 rows and 6 columns in each file. The first row describes each column, 
and the rest of the rows correspond to one ROI. The first 3 columns are the 
(x,y,z)-coordinate of each ROI (use * = MNI for volumetric MNI space, or 
* = Talairach for volumetric Talairach88 space). The fourth column is the RADIUS 
of the ROI in millimeters. The fifth column is a categorical label for each ROI. 
If using Connectome Workbench (https://www.humanconnectome.org/software/connectome-workbench),
these labels correspond to the colors in the power_surf color palette. The sixth 
column is the name of the functional network to which the ROI is assigned.

The cortical ROIs are listed first, and they are sorted by functional network. 
ROIs in the amygdala, ventral striatum, and hippocampus are grouped with their 
cortical functional networks, too. ROIs in the basal ganglia, thalamus, and 
cerebellum are listed at the end. Within each of these anatomical structures, 
the ROIs are sorted by functional network. These delineations are found in the 
file called ROIs_anatomicalLabels.txt

— — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — —
09/27/2019

After final revisions to the manuscript, a few consensus network labels changed, and each ROI 
in the basal ganglia, thalamus, and cerebellum has an “integrative” value. The integrative 
value is the percent of voxels contained within the ROI that have strong connectivity to two 
or more cortical functional networks. These changes have been added to the existing text files 
(so the 7th column is the integrative value of the ROIs).

— — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — — —
06/16/2020

A NIFTI image of the ROIs in MNI space has been added. Each voxel within an ROI has the same 
whole-number label (1-300, inclusive). These numbers correspond to the order and information 
found in ROIs_300inVol_MNI_allInfo.txt. The resolution of the image is 1x1x1mm voxels.