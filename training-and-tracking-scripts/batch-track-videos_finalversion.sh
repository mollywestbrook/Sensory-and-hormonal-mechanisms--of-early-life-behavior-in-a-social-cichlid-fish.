#!/bin/bash
#SBATCH -t 100:0:0
#SBATCH -n 1
#SBATCH -p gpu
#SBATCH --gpus=a100:1
#SBATCH --mem=60gb
#SBATCH --mail-user=mwestbr2@umd.edu
#SBATCH --mail-type=ALL

PATH="/scratch/zt1/project/juntti-lab/shared/miniconda3/bin:$PATH"
eval "$(/scratch/zt1/project/juntti-lab/shared/miniconda3/bin/conda shell.bash hook)"
conda activate sleap

# Define the path to the SLEAP executable and model paths
sleap_executable="sleap-track"
model_paths="-m /scratch/zt1/project/juntti-lab/shared/sleap/models/240726_171343.centered_instance/training_config.json -m /scratch/zt1/project/juntti-lab/shared/sleap/models/240726_171343.centroid/training_config.json"
# Define tracking parameters
tracking_params="--tracking.tracker simplemaxtracks --tracking.max_tracking 1 --tracking.max_tracks 5 --tracking.pre_cull_iou_threshold 0.8 --tracking.similarity iou --tracking.match hungarian --tracking.track_window 5"

# Get the list of video files in Google Drive folder
video_folder="/scratch/zt1/project/juntti-lab/shared/sleap/MW-videos"
video_files=$(ls "$video_folder")

# Iterate through the video files and run SLEAP commands
for video_file in $video_files; do
    if [[ $video_file == *.mp4 ]]; then
        input_video_path="$video_folder/$video_file"
        output_sleap_path="/scratch/zt1/project/juntti-lab/shared/sleap/predictions/predictions_MW_$video_file.slp"
        command="$sleap_executable $input_video_path $tracking_params -o $output_sleap_path $model_paths"
        echo "Processing $input_video_path..."
        $command
    fi
done
