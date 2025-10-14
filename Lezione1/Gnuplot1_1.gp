# gnuplot script

# --- Interactive preview ---

set terminal qt

set key box
set title "Truncation error of exponential series"

set xlabel "x"
set ylabel "Absolute error"

# set logscale x

p "e_approx.dat" u 1:2 w lp title "N = 1"

pause -1 "Press Enter to continue and save plot"

# --- Save to file ---
set terminal pdfcairo enhanced
set output 'Confronto convergenze'
replot
unset output
