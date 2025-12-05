# Set plot style
plot_style = (exists("MY_STYLE") ? MY_STYLE : "lp" )

# Set plot label
plot_label = (exists("PLOT_LABEL") ? PLOT_LABEL : "data")

# Plot one dataset

plot_command = sprintf("plot '%s' with %s title '%s'", MY_DATAFILE, plot_style, plot_label )
