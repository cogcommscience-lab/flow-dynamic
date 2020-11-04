#
#SBATCH -J fmriprep
#SBATCH --time=120:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH -p high2 # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o /data/slurm_logs/fmriprep/%x-%A-%a.out
#SBATCH -e /data/slurm_logs/fmriprep/%x-%A-%a.err
# ------------------------------------------

STUDY_DIR="/data/"
DERIVS_DIR="bids_nii/derivatives/fmriprep-latest"

# Prepare some writeable bind-mount points.
TEMPLATEFLOW_HOST_HOME=$HOME/.cache/templateflow
FMRIPREP_HOST_CACHE=$HOME/.cache/fmriprep
mkdir -p ${TEMPLATEFLOW_HOST_HOME}
mkdir -p ${FMRIPREP_HOST_CACHE}

# Prepare derivatives folder
mkdir -p ${STUDY_DIR}/${DERIVS_DIR}
mkdir -p /data/fmriprepwork

# Make sure FS_LICENSE is defined in the container.
export SINGULARITYENV_FS_LICENSE=$HOME/license.txt

# Designate a templateflow bind-mount point
export SINGULARITYENV_TEMPLATEFLOW_HOME="/templateflow"

SINGULARITY_CMD="singularity run --cleanenv -B $STUDY_DIR:/data -B ${TEMPLATEFLOW_HOST_HOME}:${SINGULARITYENV_TEMPLATEFLOW_HOME} -B /group/rwhuskeygrp/fmriprepwork:/work /data/singularity_images/fmriprep-latest.simg"

# Parse the participants.tsv file and extract one subject ID from the line corresponding to this SLURM task.
subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${STUDY_DIR}/bids_nii/participants.tsv )

# Compose the command line
cmd="${SINGULARITY_CMD} /data/bids_nii /data/${DERIVS_DIR} participant --participant-label $subject -w /work/ -vv --nthreads $SLURM_CPUS_PER_TASK --mem_mb 55000 --output-spaces MNI152NLin6Asym:res-2 MNI152NLin2009cAsym --use-aroma --fs-license-file /data/license.txt"

# Setup done, run the command
echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
module load singularity
eval $cmd
exitcode=$?

# Output results to a table
echo "sub-$subject   ${SLURM_ARRAY_TASK_ID}    $exitcode" \
      >> ${SLURM_JOB_NAME}.${SLURM_ARRAY_JOB_ID}.tsv
echo Finished tasks ${SLURM_ARRAY_TASK_ID} with exit code $exitcode
exit $exitcode

