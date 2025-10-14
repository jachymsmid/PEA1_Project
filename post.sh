#!/bin/bash

# Input data file
DATAFILE="output_t_0"

# Folder with data files
FOLDER="sim/"

# Run gnuplot commands
gnuplot -persist <<EOF
  set terminal pngcairo size 800,600 enhanced font 'Arial,12'
  set output '${FOLDER}${DATAFILE}.png'

  # Set lables
  set title "Surface Projection (Heatmap)"
  set xlabel "X"
  set ylabel "Y"

  # Top-down view
  set view map
  set pm3d map
  set palette rgb 33,13,10
  set cbrange [0:100]

  # Optional: axis ranges
  set xrange [-1.5:1.5]
  set yrange [-1.5:1.5]

  # Plot matrix data as heatmap
  splot '${FOLDER}${DATAFILE}.dat' matrix with pm3d

  unset output
EOF

echo "Plot saved as ${FOLDER}${DATAFILE}"
