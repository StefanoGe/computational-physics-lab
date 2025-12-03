# one.gnu
# Gnuplot template for plotting one dataset

# Variables injected via -e
# MY_XLABEL, MY_YLABEL, MY_TITLE, MY_XLOGSCALE, MY_YLOGSCALE, MY_DATAFILE, 
# OUTPUT_NAME, MY_STYLE, MY_TERMINAL, etc...

set grid
set xlabel (exists("MY_XLABEL") ? MY_XLABEL : "x")
set ylabel (exists("MY_YLABEL") ? MY_YLABEL : "y")
set title  (exists("MY_TITLE")  ? MY_TITLE  : "Default title")

if(exists("MY_LOGSCALE")){
	logsc_command = sprintf("set logscale %s", MY_LOGSCALE )
	eval(logsc_command)
}

# Set plot style
plot_style = (exists("MY_STYLE") ? MY_STYLE : "lp" )

# Set plot label
plot_label = (exists("PLOT_LABEL") ? PLOT_LABEL : "data")

# Plot one dataset

plot_command = sprintf("plot '%s' with %s title '%s'", MY_DATAFILE, plot_style, plot_label )

exec_mode = "default"
if( exists("MY_TERMINAL") ){ exec_mode = "custom" }

if (exec_mode eq "custom") {

    eval(sprintf("set terminal %s font 'CMU Serif, 12' enhanced", MY_TERMINAL))

    if (exists("OUTPUT_NAME")) {
        set output OUTPUT_NAME
    }

    eval(plot_command)
    unset output

    if (exists("MY_TERMINAL") && MY_TERMINAL eq "qt") {
        pause mouse close
    }
}


if(exec_mode eq "default")
{
	set terminal qt font 'CMU Serif, 12' enhanced
	eval(plot_command)
	pause mouse close
	
	if(exists("OUTPUT_NAME"))
	{
		set terminal pdfcairo font "CMU Serif,12" enhanced
		set output OUTPUT_NAME
		eval(plot_command)
		unset output
	}
}
