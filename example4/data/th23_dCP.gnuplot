###############################################################
#use this gnuplot script with commend "gnuplot th23_dCP.gnuplot"
#then context the output tex file by "context th23_dCP.tex"
###############################################################
reset
set terminal context color standalone size 12cm,9cm font ",17"
set output 'S4_th23_dCP.tex'
set colorsequence classic

set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot

set cntrparam levels discrete 9.21  # Draw contours at 1,2,3,4,5 sigma
set tmargin at screen 0.9; set bmargin at screen 0.2; set lmargin at screen 0.2; set rmargin at screen 0.8
set xtics offset 0,0.8
set ytics offset -29.5
set xlabel "True $\\theta_{23}$ [$^\\circ$]" offset 0,.6
set label "True $\\delta$ [$^\\circ$]" at graph -0.26, 0.25 rotate by 90
unset ztics


set label at 49.6,215 " " point pointtype 3 lc 0 pointsize 2.
set label "G.F. " at 50,215 left

set arrow from 47.,198 to 46,198 linewidth 2 lc rgb "black" nohead
set label "MOMENT" at 45.5,198 right
#set arrow from 47.,170 to 46,170 linewidth 2 lc rgb "black" dt 4 nohead
#set label "DUNE" at 45.5,170 right
set label "$99\\%$ C.L." at 41,145 left

set label "$S_4$ Modular" at 41,360 left




splot [40.:53][125:392] "MM_th23_dCP(MOMENT).dat" u 1:2:3 w l lt 6 lw 2 notitle
