#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:20:00
#SBATCH --job-name="BDSKY"
#SBATCH --array=1-200
#SBATCH --mem-per-cpu=4096

i=$SLURM_ARRAY_TASK_ID

java --module-path=/cluster/software/stacks/2024-06/spack/opt/spack/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/javafx-20.0.1-3yh2tfmjuuvehv2uhgfvswuwzcushdow/lib \
     --add-modules javafx.controls,javafx.fxml -jar calibratedcpp-beast-0.0.1.jar -seed $i \
     -statefile mcmc_inference_$i.state \
     -D Index=$i \
     -overwrite mcmc_inference.xml
