
# Detects number of columns 
stats MY_DATAFILE nooutput
N = STATS_columns

plot_command = ""
do for [i=2:N]{
	entry = sprintf("'%s' using 1:%d with lines", MY_DATAFILE, i, i)
	
	if (i == 2) {
	plot_command = "plot " . entry
	} else {
	plot_command = plot_command . ", " . entry
	}
}
