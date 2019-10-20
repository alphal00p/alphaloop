#!/usr/bin/gnuplot
#
# Plot some photon flux density with the epslatex terminal and change the font
# size of the labels
#
# AUTHOR: Hagen Wierstorf
# gnuplot 4.6 patchlevel 0

reset

# epslatex
set terminal epslatex size 9cm,7cm color colortext 8 standalone header \
"\\newcommand{\\ft}[0]{\\footnotesize} \\newcommand{\\ang}[1]{\\theta=\\frac{#1\pi}{10}}"  
set output 'sample_plot.tex'

# color definitions
set border linewidth 1.5
set style line 1 lc rgb '#800000' lt 1 lw 2
#set style line 2 lc rgb '#ff0000' lt 1 lw 2
#set style line 3 lc rgb '#ff4500' lt 1 lw 2
#set style line 4 lc rgb '#ffa500' lt 1 lw 2
#set style line 5 lc rgb '#006400' lt 1 lw 2
#set style line 6 lc rgb '#0000ff' lt 1 lw 2
#set style line 7 lc rgb '#9400d3' lt 1 lw 2

unset key

# Axes
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

set mxtics 2
set mytics 10

set title 'Shear Force'
set xrange [0.5:9.5]
#set xtics 0.1,0.2,1.0
set xlabel '$\frac{r}{H}$' 
set ylabel '$\frac{S\left(\frac{r}{H}, n\frac{\pi}{20}\right)}{mg}$'
#offset 1.5,-0.5 lc lt 1 
#set logscale y
set yrange [0.9:1.1]
#set format '$%g$'
#set format y '$%T$'

set datafile separator ","

# getting slope for text placing
set label 2 '\ft $NLO$'       at 0.070,0.010    rotate by  -0.0 center tc ls 1
#set label 3 '\ft $n=1$'       at 0.070,0.040    rotate by -17.0 center tc ls 2
#set label 4 '\ft $n=2$'       at 0.070,0.070    rotate by -30.0 center tc ls 3
#set label 5 '\ft $n=3$'       at 0.070,0.095    rotate by -40.0 center tc ls 4
#set label 6 '\ft $n=4$'       at 0.070,0.125    rotate by -45.0 center tc ls 5
#set label 7 '\ft $n=6$'       at 0.075,0.160    rotate by -55.0 center tc ls 6
#set label 8 '\ft $n=10$'      at 0.080,0.195    rotate by -60.0 center tc ls 7

 
plot  'ddAAA_cuhre_2M.csv'  u 1:2:3 w errorbar,
#     'Shear_angle_10.csv'  u 1:3 w l ls 2, \
#     'Shear_angle_10.csv'  u 1:4 w l ls 3, \
#     'Shear_angle_10.csv'  u 1:5 w l ls 4, \
#     'Shear_angle_10.csv'  u 1:6 w l ls 5, \
#     'Shear_angle_10.csv'  u 1:8 w l ls 6, \
#     'Shear_angle_10.csv'  u 1:12 w l ls 7

set output # finish the current output file
system('dvilualatex sample_plot.tex && dvips sample_plot.dvi && ps2eps sample_plot.ps')
system('rm sample_plot.dvi && rm sample_plot.log && rm sample_plot.aux')
system('rm sample_plot.tex && rm sample_plot.ps && rm sample_plot-inc*')

