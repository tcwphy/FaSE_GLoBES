reset
set terminal context color standalone size 12cm,9cm font ",17"
set output 'x_eta.tex'
set colorsequence classic


set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot
set cntrparam levels discrete 11.83  # Draw contours at 3 sigma

set tmargin at screen 0.9; set bmargin at screen 0.2; set lmargin at screen 0.2; set rmargin at screen 0.8
set ytics offset -29.5
unset ztics

set xlabel "$x$" offset 0,0.6
set xtics offset 0,0.8
set label "$\\eta/\\pi$" at graph -0.2,0.4 rotate by 90
set label at -3.65029,1.13067 " " point pointtype 7 lc 0 pointsize 2.

set arrow from -4.,1.85 to -3.2,1.85 linewidth 2 lc rgb "black" nohead
set label "MOMENT" at -4.1,1.85 right
#set arrow from -4.,1.75 to -3.2,1.75 linewidth 2 dt 4 lc rgb "black" nohead
#set label "DUNE" at -4.1,1.75 right

set label "$3\\sigma$; TDLS" at -8.7,.95 left


splot[-9:-3][0.8:2]\
"constraint_x_eta_NT(MOMENT).dat" u 1:2:3 w l lt 6 lw 2 notitle


