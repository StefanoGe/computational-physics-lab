# style_common.gp: gp module containing shared styling settings
# Use: inside main.gp load "path/style_common.gp"

set encoding utf8

# Set labels and title
set xlabel (exists("MY_XLABEL") ? MY_XLABEL : "x")
set ylabel (exists("MY_YLABEL") ? MY_YLABEL : "y")
set title  (exists("MY_TITLE")  ? MY_TITLE  : "Default title") font ',20'

# Setting ticks and border 
set mxtics 1
set mytics 1
set grid xtics ytics mxtics mytics lw 0.4 lc rgb '#AAAAAA'
set border 3 lw 0.5 lc rgb "black"
set xtics nomirror scale 1
set ytics nomirror scale 1

# Style preset
set style line 1 lc rgb '#66C2A5' lw 3 pt 7 ps 1.5

# Set logscale
if(exists("MY_LOGSCALE")){
	logsc_command = sprintf("set logscale %s", MY_LOGSCALE )
	eval(logsc_command)
}
