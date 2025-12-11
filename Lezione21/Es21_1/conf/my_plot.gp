set encoding utf8
set mxtics 1
set mytics 1
set grid xtics ytics mxtics mytics lw 0.4 lc rgb '#AAAAAA'
set border 3 lw 0.5 lc rgb "black"
set xtics nomirror scale 1
set ytics nomirror scale 1
set style line 1 lc rgb '#66C2A5' lw 3 pt 7 ps 1.5
set xlabel 'x'
set ylabel 'y'
set terminal qt font 'CMU Serif, 12' enhanced
plot 'data/my_plot.dat' with lp ls 1 title ''

pause mouse close
