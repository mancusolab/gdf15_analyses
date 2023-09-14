#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=16Gb
#SBATCH --array=1-4

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$(($1 - 1))
else
    idx=$(($SLURM_ARRAY_TASK_ID - 1))
fi

module load r

ddir=../merge/
rdir=../coloc/

OLINK=${ddir}/OLINK_GDF15_Chr19_subset.cis.1Mb.23andMe.merged.tsv
ROCHE=${ddir}/GS_GDF15_Chr19_subset.cis.1Mb.23andMe.merged.tsv

LD=(UKBB UKBB 1000G 1000G)
DATA=($OLINK $ROCHE $OLINK $ROCHE)
LABEL=(OLINK ROCHE OLINK ROCHE)

Rscript coloc.R \
    ${DATA[$idx]} \
    ${LD[$idx]} \
    ${rdir}/${LABEL[$idx]}.GDF15.23andMe.${LD[$idx]}LD
