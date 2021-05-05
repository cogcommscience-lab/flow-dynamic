#!/bin/bash
# 
#SBATCH -J xcpengine
#SBATCH --time=120:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH -p high2 # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o /data/slurm_logs/xcpengine/%x-%A-%a.out
#SBATCH --array=1-139
# ------------------------------------------
DATA_DIR="/data/bids_nii/derivatives/fmriprep-latest"
FULL_COHORT="$DATA_DIR/fmriprep/cohort.csv"

SINGULARITY_IMG="/data/singularity_images/xcpEngine.simg"
LINE_NUM=$( expr ${SLURM_ARRAY_TASK_ID} + 1 )
LINE=$(awk "NR==$LINE_NUM" $FULL_COHORT)
mkdir -p ${DATA_DIR}/tempcohort
TEMP_COHORT=$DATA_DIR/tempcohort/${SLURM_ARRAY_TASK_ID}.csv
HEADER=$(head -n 1 $FULL_COHORT)
echo $HEADER > $TEMP_COHORT
echo $LINE >> $TEMP_COHORT

SINGULARITY_CMD="singularity run -B $DATA_DIR:/home/user/data $SINGULARITY_IMG"

# Compose the command line
cmd="${SINGULARITY_CMD} -d /home/user/data/fmriprep/fc-aroma-gsr-0407.dsn -c /home/user/data/tempcohort/${SLURM_ARRAY_TASK_ID}.csv -o /home/user/data/output_0407_gsr -r /home/user/data/ -i $TMPDIR"

# Setup done, run the command
echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
module load singularity
eval $cmd
exitcode=$?

# Output results to a table
echo " ${SLURM_ARRAY_TASK_ID}    $exitcode" \
      >> ${SLURM_JOB_NAME}.${SLURM_ARRAY_JOB_ID}.tsv
echo Finished tasks ${SLURM_ARRAY_TASK_ID} with exit code $exitcode
exit $exitcode
