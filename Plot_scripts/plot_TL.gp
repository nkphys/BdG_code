set term epslatex standalone color 
unset key
# unset xtics
# unset ytics
set size square

set xrange [-0.5:45.5]
set yrange [-0.5:45.5]
set cbrange [-1.0:1.0] 
set output 'plot_TL.tex'

set border 4095 lw 3 
set style arrow 1 head filled size screen 0.02,20 lw 10 lc pal 
set style arrow 3 head filled size screen 0.03,5,5  lw 3 lc pal 
scale = 0.7 
p 'MySkyrmion.txt' u (( $1 + (0.5*$2)  )-0.4*$4):(( sqrt(3.0)*0.5*($2) )-0.4*$5):(($4)*scale):(($5)*scale):(($3)) w vec arrowstyle 3 notitle
