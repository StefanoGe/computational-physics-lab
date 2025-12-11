set encoding utf8
set mxtics 1
set mytics 1
set grid xtics ytics mxtics mytics lw 0.4 lc rgb '#AAAAAA'
set border 3 lw 0.5 lc rgb "black"
set xtics nomirror scale 1
set ytics nomirror scale 1
set style line 1 lc rgb '#66C2A5' lw 3 pt 7 ps 1.5
set title 'Numerical solution of , u(0)=2' font ',20'
set xlabel 't'
set ylabel 'u(t)'
set terminal qt font 'CMU Serif, 12' enhanced
plot 'data/euler_test1.dat' with lp ls 1 title ''

pause mouse close
set terminal pdfcairo font "CMU Serif,12" enhanced
set output 'plot/euler_test1.pdf'
plot 'data/euler_test1.dat' with lp ls 1 title ''

unset output
