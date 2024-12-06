#!/bin/bash
#SBATCH -t 60:0:0
#SBATCH -n 1
#SBATCH -p gpu
#SBATCH --gpus=a100:1
#SBATCH --mem=60gb
#SBATCH --mail-user=mwestbr2@umd.edu
#SBATCH --mail-type=ALL


eval "$(conda shell.bash hook)"

conda activate sleap

cd /scratch/zt1/project/juntti-lab/shared/sleap

sleap-train /scratch/zt1/project/juntti-lab/shared/sleap/trainingpackages/sleapmodel1_color_labels_07-26-2024.v002.slp.training_job/centroid.json /scratch/zt1/project/juntti-lab/shared/sleap/trainingpackages/sleapmodel1_color_labels_07-26-2024.v002.slp.training_job/sleapmodel1_color_labels.v002.pkg.slp

