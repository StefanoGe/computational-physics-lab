
set terminal qt font 'CMU Serif, 12' enhanced
eval(plot_command)
pause mouse close

if (exists("OUTPUT_NAME")){
	set terminal pdfcairo font "CMU Serif,12" enhanced
	set output OUTPUT_NAME
	eval(plot_command)
	unset output
}
