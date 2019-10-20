#!/usr/bin/gnuplot

reset

# epslatex
set terminal epslatex size 15cm,15cm color colortext 8 standalone header \
"\\newcommand{\\ft}[0]{\\footnotesize} \\newcommand{\\ang}[1]{\\theta=\\frac{#1\pi}{10}}"  
set output 'sample_plot.tex'

# color definitions
set border linewidth 1.5
set style line 1 lc rgb '#800000' lt 1 lw 2 #dark red
set style line 2 lc rgb '#ff0000' lt 1 lw 2 #red
set style line 3 lc rgb '#ff4500' lt 1 lw 2 #orange
set style line 4 lc rgb '#ffa500' lt 1 lw 2 #yellow
set style line 5 lc rgb '#006400' lt 1 lw 2 #green
set style line 6 lc rgb '#0000ff' lt 1 lw 2 #blue
set style line 7 lc rgb '#9400d3' lt 1 lw 2 #purple

#dot size
#set pointsize .5

#unset key
set key bottom right box


# Axes
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

set mxtics 2
set mytics 10

set title 'TITLE_NAME'
set autoscale
#set xrange [0.5:9.5]
#set xtics 0.1,0.2,1.0
set xlabel '$\delta$' 
set ylabel 'integrand scaling'
#offset 1.5,-0.5 lc lt 1 
#set logscale y
#set yrange [0.9:1.1]
#set format '$%g$'
#set format y '$%T$'
set logscale 

set datafile separator ","
plot 'data/LIMIT_FILE'  u 1:20 w l ls 6 title "Sum",\
      for [col=3:13:2] 'data/LIMIT_FILE' u 1:col title  sprintf("cut %d", col/2),\

set output # finish the current output file
system('dvilualatex sample_plot.tex && dvips sample_plot.dvi && ps2eps sample_plot.ps')
system('rm sample_plot.dvi && rm sample_plot.log && rm sample_plot.aux')
system('rm sample_plot.tex && rm sample_plot.ps && rm sample_plot-inc*')

