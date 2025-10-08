#!/bin/bash

# Input data file
DATAFILE="data.dat"

# Output image
OUTPUT="contour.png"

# Run gnuplot commands
gnuplot <<EOF
  set terminal pngcairo size 800,600 enhanced font 'Arial,12'
  set output '${OUTPUT}'

  set title "Surface Projection (Heatmap)"
  set xlabel "X"
  set ylabel "Y"

  # Top-down view
  set view map
  set pm3d map
  set palette rgb 33,13,10

  # Optional: axis ranges
  # set xrange [0:50]
  # set yrange [0:50]

  # Plot matrix data as heatmap
  splot '${DATAFILE}' matrix with pm3d

  unset output
EOF

echo "Plot saved as ${OUTPUT}"
