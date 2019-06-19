#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N alb_h_1
#$ -q omni
#$ -pe sm 4
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

source activate selection

/home/jmanthey/discoal/discoal 7 1000 220000 \
-Pt 663.2616 2653.0464 \
-Pre 726 16588 \
-ws 0 \
-Pa 601.2815 60128.1500 \
-Pu 0 0.0021 \
-x 0.13636363636363635 \
-en 0.000632203 0 2.019974718 -en 0.000869762 0 2.972900594 -en 0.001390269 0 3.018899092 \
-en 0.001675123 0 2.670204345 -en 0.001977748 0 2.258370383 -en 0.003003693 0 2.031222469 \
-en 0.003389207 0 1.87704179 -en 0.003798772 0 2.031954987 -en 0.005187257 0 2.26223803 \
-en 0.00685217 0 3.152411573 -en 0.00814243 0 4.499390796 -en 0.010395677 0 6.300326634 \
-en 0.01214188 0 8.671650017 -en 0.01519136 0 11.63822226 -en 0.018847949 0 14.89414223 \
-en 0.037163017 0 16.59131765 -en 0.045193853 0 12.06869701 -en 0.054823515 0 1.085728395 \
-en 0.066370308 0 0.840997371 -en 0.080215905 0 1.121015277 \
> alb_hard_1.msOut
