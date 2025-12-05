
eval(sprintf("set terminal %s font 'CMU Serif, 12' enhanced", MY_TERMINAL))

if (exists("OUTPUT_NAME")) {
	set output OUTPUT_NAME
}

eval(plot_command)
unset output

if (exists("MY_TERMINAL") && MY_TERMINAL eq "qt") {
	pause mouse close
}
