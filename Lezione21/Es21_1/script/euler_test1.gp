set encoding utf8
set mxtics 1
set mytics 1
set grid xtics ytics mxtics mytics lw 0.4 lc rgb '#AAAAAA'
set border 3 lw 0.5 lc rgb "black"
set xtics nomirror scale 1
set ytics nomirror scale 1
set key box lw 0.5 lc rgb 'black' textcolor rgb 'black' spacing 1.5 opaque width 0.5
set for [i=1:8] linetype i lw 3 pt 7 ps 1.5
set linetype 1 lc rgb "#66C2A5"
set linetype 2 lc rgb "#FC8D62"
set linetype 3 lc rgb "#8DA0CB"
set linetype 4 lc rgb "#E78AC3"
set linetype 5 lc rgb "#A6D854"
set linetype 6 lc rgb "#FFD92F"
set linetype 7 lc rgb "#E5C494"
set linetype 8 lc rgb "#B3B3B3"
set title "Numerical solution of u'(t) = -2tu, u(0)=2" font ',20'
set xlabel "t"
set ylabel "u(t)"
set terminal qt size 800,500 font 'CMU Serif, 12' enhanced
plot "data/euler_test1_0_0.dat" with points title "Numerical solution", \
"data/euler_test1_0_1.dat" with lines title "Exact function"

pause mouse close
set terminal pdfcairo font "CMU Serif,12" enhanced
set output 'plot/euler_test1.pdf'
plot "data/euler_test1_0_0.dat" with points title "Numerical solution", \
"data/euler_test1_0_1.dat" with lines title "Exact function"

unset output
