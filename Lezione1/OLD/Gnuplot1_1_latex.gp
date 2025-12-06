# gnuplot script

# --- Interactive preview ---

set terminal qt font "CMU Serif, 12" enhanced

set key left box
set title "Truncation error of exponential series"

set xlabel "x"
set ylabel "Absolute error"

#set logscale x
#set logscale y

#set yrange [0:2]

plot for [N = 1:4] "e_approx.dat" using 1:(column(N)+1) with linespoints title "N = ".N	

pause mouse close
#pause -1 "Press Enter to continue and save plot"

# --- Save to file ---
set terminal pdfcairo font "CMU Serif, 12" enhanced
set output 'Confronto convergenze in valore assoluto.pdf'
replot
unset output

unset logscale x
unset logscale y

# gnuplot script

# --- Interactive preview ---

set terminal qt font "CMU Serif, 12" enhanced

set key box
set title "Truncation error of exponential series - Error over x^{N+1}/N+1!"

set xlabel "x"
set ylabel "Relative error"

unset yrange

plot for [N = 1:4] "rel_errors.dat" using 1:(column(N)+1) with linespoints title "N = ".N	

# metto due pause perché se no mi ruba l'enter 

pause mouse close
#pause -1 "Press Enter to continue and save plot"

# --- Save to file ---
set terminal pdfcairo font "CMU Serif, 12" enhanced
set output 'Confronto convergenze relative.pdf'
replot
unset output
