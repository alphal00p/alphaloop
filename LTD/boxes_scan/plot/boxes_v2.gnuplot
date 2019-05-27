################################################################################
# Gnuplot file
################################################################################


reset
set encoding iso_8859_1

set lmargin 10
set rmargin 0

#set terminal postscript portrait enhanced mono dashed lw 1.0 "Helvetica" 14 
set terminal postscript portrait enhanced color dashed lw 1.0 "Helvetica" 14 
#set terminal postscript eps enhanced color font 'Helvetica,10

box = 'VEGAS_box'
doublebox = 'QP_doublebox'
triplebox = 'NEW_1B_triplebox'

#set size ratio 0.75 
set key font ",9"
set key samplen "2"
set output "boxesv2.ps"
set style line 100 lt 2 dt 3 lc rgb "black" lw 2

set style line 1 lt 1 lc rgb "black" lw 1.8 pt 1 ps 0
set style line 2 lt 1 lc rgb "red" lw 1.8 pt 1 ps 0
set style line 3 lt 1 lc rgb "forest-green" lw 1.8 pt 1 ps 0
set style line 4 lt 1 lc rgb "blue" lw 1.8 pt 1 ps 0
set style line 5 lt 1 lc rgb "goldenrod" lw 1.8 pt 1 ps 0
set style line 6 lt 1 lc rgb "dark-magenta" lw 1.8 pt 1 ps 0

#set style line 11 lt 2 lc rgb "black" lw 1.8
#set style line 12 lt 2 lc rgb "#530FAD" lw 1.8
#set style line 13 lt 2 lc rgb "#FD0006" lw 1.8
#set style line 14 lt 2 lc rgb "#01939A" lw 1.8

set style line 11 lt 2 lc rgb "black" lw 1.8 pt 1 ps 0
set style line 12 lt 2 lc rgb "red" lw 1.8 pt 1 ps 0
set style line 13 lt 2 lc rgb "forest-green" lw 1.8 pt 1 ps 0
set style line 14 lt 2 lc rgb "blue" lw 1.8 pt 1 ps 0
set style line 15 lt 2 lc rgb "goldenrod" lw 1.8 pt 1 ps 0
set style line 16 lt 2 lc rgb "dark-magenta" lw 1.8 pt 1 ps 0

set style line 21 lt 5 lc rgb "black" lw 1.8 pt 1 ps 0
set style line 22 lt 5 lc rgb "red" lw 1.8 pt 1 ps 0
set style line 23 lt 5 lc rgb "forest-green" lw 1.8 pt 1 ps 0
set style line 24 lt 5 lc rgb "blue" lw 1.8 pt 1 ps 0
set style line 25 lt 5 lc rgb "goldenrod" lw 1.8 pt 1 ps 0
set style line 26 lt 5 lc rgb "dark-magenta" lw 1.8 pt 1 ps 0

set style data linespoints

set multiplot
set tics front

# Fraction of points with accuracy less than target {/Symbol D} [-]
# Comparison of {/Arial Ninja} reduction accuracy for various processes

# ==========================================================================
# Main plot
# ==========================================================================

set label "Results for 4-points multi-loop scalar integrals" font ",10" at graph 0.0, graph 0.92 offset 7.0,-0.2
set xrange [-7.:100.]
set yrange [0.7:6]
set origin 0.00, 0.5
set size 1, 0.23
set bmargin 0 
set tmargin 0
set key spacing 1.5
#set xtics 5 nomirror
unset xtics
set mxtics 4
set mytics 5
set ylabel "Integral [GeV^{-2(l+1)}]" offset 1.5, 0 font ",10"
set ytics 1 font ",10"
#set xtics nomirror
#set logscale y
#set logscale x
#set format y '10^{%T}'
#set format x '10^{%T}'
# Below is the position with no letter for the label prefixes
#set key at graph 0.94, graph 0.99 noautotitles spacing 2.4 horizontal font "Arial,10" width -5.0
#set label front 'MadGraph5\_aMC\@NLO' font "Courier,11" rotate by 90 at graph 1.02, graph 0.04

#set xlabel 't [GeV]' font ",12"

#set arrow from 1.0e-3,5e-3 to 1.0e-3,1e0 nohead lc rgb 'gray'

#'{/Symbol a}_s^2{/Symbol a}', \
#'{/Symbol a}_s^2{/Symbol a} + {/Symbol a}_s^3{/Symbol a}'
#'{/Symbol a}_s^2{/Symbol a} + {/Symbol a}_s^3{/Symbol a} + {/Symbol a}_s^2{/Symbol a}^2'
plot \
\
'../'.box.'/analytic_line_result.dat' using 1:($2*1.e3) ls 2 title '1-loop, analytic (x 10{^3})',\
'../'.box.'/ltd_results.dat' using 1:($2*1.e3):($3*1.e3) with yerrorbars ls 12 title '1-loop, LTD (x 10{^3})',\
'../'.doublebox.'/analytic_line_result.dat' using 1:($2*5.e4) ls 3 title '2-loop, analytic (x 5 10{^4})',\
'../'.doublebox.'/ltd_results.dat' using 1:($2*5.e4):($3*5.e4) with yerrorbars ls 13 title '2-loop, LTD (x 5 10{^4})',\
'../'.triplebox.'/analytic_line_result.dat' using 1:(($2)*1.e6) ls 4 title '3-loop, analytic (x -10{^7})',\
'../'.triplebox.'/ltd_results.dat' using 1:((-$2)*1.e6):((-$3)*1.e6) with yerrorbars ls 14 title '3-loop, LTD (x -10{^7})'

# ==========================================================================
# Second subplot
# ==========================================================================

unset label
unset format
#set label "Normalised deviation" font ",12" at graph 0.05, graph 0.88
set label "Relative deviation" font ",10" at graph 0.02, graph 0.2
set origin 0.0, 0.4
set size 1, 0.10
set yrange [-0.5+0.05:0.5-0.05]
set ytics 0.5 font ",8"
set ytics add (".25" 0.25)
set ytics add ("-.25" -0.25)
set mytics 2
set ylabel "Deviation [%]" offset 3.0, 0 font ",10"
#set bmargin 0 
#set tmargin 0
#set mytics 5
##set xtics nomirror
##unset logscale y
#unset ylabel
##set format y '10^{%T}'
#set nologscale y
#set format x '10^{%T}'

unset key

plot \
\
0 ls 1,\
'../'.box.'/ltd_results.dat' using 1:((($5-$2)/abs($5))*100.) ls 2,\
'../'.doublebox.'/ltd_results.dat' using 1:((($5-$2)/abs($5))*100.) ls 3,\
'../'.triplebox.'/ltd_results.dat' using 1:((($5-$2)/abs($5))*100.) ls 4

# ==========================================================================
# Third subplot
# ==========================================================================

unset logscale y
unset format y
unset label
unset format
set label "Deviation significance" font ",10" at graph 0.02, graph 0.2
set origin 0.0, 0.3
set size 1, 0.10
set yrange [-7.0:7.0]
set ytics 3.0 font ",10"
set mytics 3
set xtics (-7, 0, 25, 50, 75, 100) nomirror font ",10"
set ylabel "Deviation [{/Symbol s}]" offset 1.5, 0 font ",10"
#set bmargin 0 
#set tmargin 0
#set mytics 5
##set xtics nomirror
##unset logscale y
#unset ylabel
##set format y '10^{%T}'
#set nologscale y
#set format x '10^{%T}'
#set xlabel 't [GeV]' font ",10"
set label front 't [GeV]' font ",10" rotate by 90 at graph 1.02, graph -0.05

unset key

plot \
\
0 ls 1,\
3 ls 100,\
-3 ls 100,\
'../'.box.'/ltd_results.dat' using 1:(($5-$2)/abs($3)) ls 2,\
'../'.doublebox.'/ltd_results.dat' using 1:((($5-$2)/abs($3))) ls 3,\
'../'.triplebox.'/ltd_results.dat' using 1:((($5-$2)/abs($3))) ls 4

!ps2pdf "boxesv2.ps"
!open "boxesv2.ps"
