#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N discoal_alb_hard
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-1500

source activate selection


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.045454545454545456 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.13636363636363635 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.22727272727272727 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.3181818181818182 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.4090909090909091 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.5 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.5909090909090909 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.6818181818181818 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.7727272727272727 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.8636363636363636 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut


/home/jmanthey/discoal/discoal 14 1 550000 \
-Pt 146.2941 1316.6473 \
-Pre 1815 41470 \
-ws 0 \
-Pa 79.574 7957.400 \
-Pu 0 0.01570865 \
-x 0.9545454545454546 \
-en 0.005550989 0 1.620863698 -en 0.006877349 0 4.967964833 -en 0.00978348 0 9.823375218 -en 0.014858601 0 13.98160443 \
-en 0.01676565 0 13.86850558 -en 0.020944106 0 10.22263277 -en 0.025660203 0 6.114622395 -en 0.033896155 0 4.288448104 \
-en 0.040278788 0 4.875597189 -en 0.051425095 0 7.003588505 -en 0.064790475 0 10.40588017 -en 0.075148271 0 15.17304764 \
-en 0.080816713 0 17.6515193 -en 0.0932366 0 22.15780897 -en 0.114926062 0 32.88720276 -en 0.131734789 0 49.00915295 \
-en 0.161088698 0 73.24380878 -en 0.183837158 0 111.3873055 -en 0.196286556 0 125.3681614 -en 0.308116121 0 40.81083602 \
-en 0.328319119 0 16.13400828 -en 0.349782556 0 7.719957136 -en 0.449888876 0 3.224188169 -en 0.542582783 0 4.371611812 \
>> output/${SGE_TASK_ID}-alb_hard_0.msOut





