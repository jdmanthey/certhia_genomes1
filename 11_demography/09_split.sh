#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N split1
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

source activate smcpp_mod

# run split
smc++ split --cores 24 -o split/ albescens_2/model.final.json americana_2/model.final.json --timepoints 5e3 2e6 *total.smc

