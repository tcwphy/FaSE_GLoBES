#!/usr/bin/gnuplot
reset
set terminal context color standalone size 12cm,9cm font ",17"
set output 'theta_phi.tex'
set colorsequence classic


set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot
set cntrparam levels discrete 11.83

set tmargin at screen 0.9; set bmargin at screen 0.2; set lmargin at screen 0.2; set rmargin at screen 0.8
set ytics offset -30.
unset ztics
set xtics("0.16" 0.16, "" 0.17, "0.18" 0.18, "" 0.19, "0.2" 0.2, "" 0.21, "0.22" 0.22, "" 0.23, "0.24" 0.24) rotate
set xlabel "$\\theta_\\nu/\\pi$" offset 0,0.6
set xtics offset 0,0.5
set label "$\\phi_\\nu/\\pi$" at graph -0.25,0.4 rotate by 90

set label "$3\\sigma$, WFS" at .165,.88 left
#set arrow from .233,.98 to .245,.98 linewidth 2 dt 4 lc rgb "black" nohead
#set label "DUNE" at .23,.98 right
set arrow from .233,1. to .245,1. linewidth 2 lc rgb "black" nohead
set label "MOMENT" at .23,1. right


set label at 0.209269,0.914924 " " point pointtype 7 lc 7 pointsize 2

splot[0.16:0.25][0.85:1.1]\
"theta_phi(MOMENT).dat" u 1:2:3 w l lt 6 lw 2 notitle
