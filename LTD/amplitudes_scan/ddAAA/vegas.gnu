#!/usr/bin/gnuplot

reset

# epslatex
set terminal epslatex size 9cm,5cm color colortext 8 standalone header \
"\\newcommand{\\ft}[0]{\\footnotesize} \\newcommand{\\ang}[1]{\\theta=\\frac{#1\pi}{10}}"  
set output 'sample_plot.tex'

# color definitions
set border linewidth 1.5
set style line 1 lc rgb '#800000' lt 2 lw 1 #dard red
set style line 2 lc rgb '#ff0000' lt 2 lw 1 #red
set style line 3 lc rgb '#ff4500' lt 2 lw 1 #orange
set style line 4 lc rgb '#ffa500' lt 2 lw 1 #yellow
set style line 5 lc rgb '#006400' lt 2 lw 1 #green
set style line 6 lc rgb '#0000ff' lt 2 lw 1 #blue
set style line 7 lc rgb '#9400d3' lt 2 lw 1 #purple

#dot size
set pointsize .5

#unset key
set key top right box 

# Axes
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

set mxtics 2
set mytics 10

set title 'Vegas 10x10M'
set xrange [-0.5:9.5]
#set xtics 0.1,0.2,1.0
set xlabel 'seed nr.' offset 0,0 
set ylabel 'rust/MG5' offset 0,0 
#set logscale y
set yrange [0.:2]
#set format '$%g$'
#set format y '$%T$'

set datafile separator ","
plot  'data/ddAAA_vegas_100M.csv'  u 1:2:3 w errorbar ls 2 title 'real',\
      'data/ddAAA_vegas_100M.csv'  u 1:4:5 w errorbar ls 6 title 'imag',

set output # finish the current output file
system('dvilualatex sample_plot.tex && dvips sample_plot.dvi && ps2eps sample_plot.ps')
system('rm sample_plot.dvi && rm sample_plot.log && rm sample_plot.aux')
system('rm sample_plot.tex && rm sample_plot.ps && rm sample_plot-inc*')

