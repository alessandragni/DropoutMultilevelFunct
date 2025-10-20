#!/bin/bash

# Parse input parameters

Is=(10 30 50)
Js=(4 10 50) 
balanceds=("TRUE" "FALSE")
grids=(200 1000 5000) 

mkdir -p logs

for I in "${Is[@]}"; do
  for J in "${Js[@]}"; do
    for balanced in "${balanceds[@]}"; do
      for grid in "${grids[@]}"; do
          echo "Submitting job with I=$I, J=$J, grid=$grid and balanced=$balanced"
          output_file="logs/output_I${I}_J${J}_grid${grid}_balanced${balanced}.txt"
          error_file="logs/error_I${I}_J${J}_grid${grid}_balanced${balanced}.txt"
          qsub -v I=$I,J=$J,balanced=$balanced,grid=$grid \
               -o "$output_file" \
               -e "$error_file" \
               sim.job &
      done
    done
  done
done
