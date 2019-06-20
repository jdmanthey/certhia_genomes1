#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N training_1
#$ -q omni
#$ -pe sm 8
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

python ~/diploSHIC/diploSHIC.py train training_albescens/ training_albescens/ albescensModel
