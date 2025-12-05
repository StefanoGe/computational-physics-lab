load "style_common.gp"

if( MODE == "one_datafile" ){load "one.gp"}

exec_mode = "default"
if( exists("MY_TERMINAL") ){ exec_mode = "custom" }

if (exec_mode eq "custom"){ load "terminal_custom.gp"}

if(exec_mode eq "default"){ load "terminal_default.gp"}
