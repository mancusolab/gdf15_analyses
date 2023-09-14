#!/bin/bash

module load r

ddir=../merge/
rdir=../mr/

OLINK=${ddir}/OLINK_GDF15_Chr19_subset.cis.1Mb.23andMe.merged.tsv
ROCHE=${ddir}/GS_GDF15_Chr19_subset.cis.1Mb.23andMe.merged.tsv

DATA=($OLINK $ROCHE $OLINK $ROCHE)
LD=(UKBB UKBB 1000G 1000G)
OUT=(OLINK ROCHE OLINK ROCHE)

for idx in {0..3} ; do
Rscript MR.R \
    ${DATA[$idx]} \
    ${LD[$idx]} \
    ${OUT[$idx]} \
    ${rdir}/${OUT[$idx]}.GDF15.23andMe.${LD[$idx]}LD

Rscript MR.cond.R \
    ${DATA[$idx]} \
    ${LD[$idx]} \
    ${OUT[$idx]} \
    ${rdir}/${OUT[$idx]}.GDF15.23andMe.${LD[$idx]}LD.cond
done
