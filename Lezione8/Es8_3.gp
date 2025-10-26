# Example plot with hackaday.io-like color theme
reset
set terminal qt font "CMU Serif, 12" enhanced #size 650,412 background rgb "#151515"

exp_sin(x) = exp(sin(x-1))

plot "data/expsin_fitted_plot.dat", exp_sin(x)

pause mouse close

#set tics textcolor rgb "#ffffff"
#set border lc rgb "#ffffff"
# cyan, magenta, and yellow for the first three line styles
#set style line 1 lc rgb "#00ffff"
#set style line 2 lc rgb "#ff00ff"
#set style line 3 lc rgb "#ffff00"
#set mxtics
#set mytics
#set my2tics
#set grid
#set xtics
#set ytics nomirror
#set ytics
#set y2tics
#set key textcolor rgb "#ffffff" bottom right
#set title textcolor rgb "#ffffff"
#set xlabel "Current ({/Symbol m}A)" textcolor rgb "#ffffff" font ',15'
#set ylabel "Luminous Intensity (Relative)" textcolor rgb "#ffffff" font ',15'
#set y2label "Overall Efficiency (Relative)" textcolor rgb "#ffffff" font ',15'
#set title "Luxeon Z (cyan) Relative Light Output" font ',20'
#stats 'luxeon_z_cyan.dat' u 1 nooutput name 'voltage'
#stats 'luxeon_z_cyan.dat' u 2 nooutput name 'current'
#stats 'luxeon_z_cyan.dat' u 3 nooutput name 'lux'
#plot 'luxeon_z_cyan.dat' using ($2*1e6):($3/lux_max) with linespoints ls 1 title 'intensity', \
#      'luxeon_z_cyan.dat' using ($2*1e6):($3/($1*$2))/(lux_max/(voltage_max*current_max)) with linespoints ls 2 title 'efficiency' axis x1y2
