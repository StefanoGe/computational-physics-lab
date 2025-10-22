# gnuplot script

# --- Interactive preview ---

set terminal qt font "CMU Serif, 12" enhanced

set key box
set title "Truncation error of exponential series"

set xlabel "x"
set ylabel "Absolute error"

#set logscale x
#set logscale y

#set yrange [0:2]

set multiplot

plot [0: 1]for [N = 1:1] "e_approx.dat" using 1:(column(N + 1)) with linespoints title "N = ".N, \
		x**(1+1)/(1+1)! with linespoints

pause mouse close

# --- Save to file ---
set terminal pdfcairo enhanced
set output 'Confronto convergenze in valore assoluto.pdf'
replot
unset output

unset multiplot

unset logscale x
unset logscale y

# gnuplot script

# --- Interactive preview ---

set terminal qt font "CMU Serif, 12" enhanced

set key box
set title "Truncation error of exponential series - Error over $\frac{x^{N+1}}{N+1}!$"

set xlabel "x"
set ylabel "Relative error"

unset yrange

plot for [N = 1:4] "rel_errors.dat" using 1:(column(N)+1) with linespoints title "N = ".N	

pause mouse close

# --- Save to file ---
set terminal pdfcairo enhanced
set output 'Confronto convergenze relative.pdf'
replot
unset output
