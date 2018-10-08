
# coding: utf-8

# In[ ]:


# Details:
# Study: Brain Network Organization and Behavior
# Note: This notebook merges single volume nii files into a 4d file. This is the 4d file used in main analysis


# Credits:
# Written by Richard Huskey


# Notes and dependencies
# Requires job scheduler that accepts .pbs files


# In[ ]:


# Makes a template file in the working directory

with open('fslmerge_template.txt', 'w') as readtemp:
    readtemp.write('#!/bin/bash\n#PBS -l walltime=240:00:00\n#PBS -l nodes=1:ppn=1:huskey,mem=16000mb\n#PBS -N fslmerge_sub{sub}_run{run}\n#PBS -j oe\n#PBS -m abe\n#PBS -m abe -M huskey.29@osu.edu\n\n#COMMANDS TO RUN\nmodule load fsl\nfslmerge -tr /fs/project/huskey.29/flow_dynamic/bids_nii_trimmed/merged/sub-{sub}_task-game_run-{run}_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc_merged.nii.gz /fs/project/huskey.29/flow_dynamic/bids_nii_trimmed/volumes/sub-{sub}_task-game_run-{run}_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc_{vols} 2.0')


# In[ ]:


# !/usr/bin/python

import os,subprocess,sys

# Define bids subject range
SUB = range(5,43)
RUN = range(1,4)
SUBDICT = [10, 27, 38]
VOLS = '{0010..0125}'
PBSTEMPLATE = '/fs/project/huskey.29/flow_dynamic/bids_nii_trimmed/pbs_fslsplit/fslmege_template.txt'

# This loop populates obs scripts that will merge each nii file
pbsscripts = []

for sub in SUB:
    for run in RUN:
            if sub in SUBDICT:
                continue
            with open(PBSTEMPLATE, 'r') as ptemplate:
                pstring = ptemplate.read()
                pbs = pstring.format(sub='%0.3d' % sub,run='%0.2d' %run)
                poutname = '/fs/project/huskey.29/flow_dynamic/bids_nii_trimmed/pbs_fslmerge/sub-{}_task-game_run-{}_bold.pbs'.format('%0.3d' % sub,'%0.2d' %run)
                with open(poutname, 'w') as poutscript:
                    poutscript.writelines(pbs)

