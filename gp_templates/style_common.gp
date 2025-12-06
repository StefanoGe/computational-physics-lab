
set xlabel (exists("MY_XLABEL") ? MY_XLABEL : "x") textcolor rgb "black
set ylabel (exists("MY_YLABEL") ? MY_YLABEL : "y")
set title  (exists("MY_TITLE")  ? MY_TITLE  : "Default title") font ',30'

if(exists("MY_LOGSCALE")){
	logsc_command = sprintf("set logscale %s", MY_LOGSCALE )
	eval(logsc_command)
}
