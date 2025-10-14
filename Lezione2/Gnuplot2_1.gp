# gnuplot script

# --- Interactive preview ---

set terminal qt

set key box
set title "Truncation error comparison"

set xlabel "N"
set ylabel "True value deviation"

# set logscale x

p [2000:6000] [] "data.dat" u 1:2 w lp title "Right order", \
	"data.dat" u 1:3 w lp title "Reversed order"

pause -1 "Press Enter to continue and save plot"

# --- Save to file ---
set terminal pdfcairo enhanced
set output 'Confronto convergenze'
replot
unset output
