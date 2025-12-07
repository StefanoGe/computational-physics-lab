load sprintf("%s/style_common.gp", TMP_PATH)


if( MODE eq "one" ){load sprintf("%s/one.gp", TMP_PATH)}

if(MODE eq "multi"){load sprintf("%s/multi_columns.gp", TMP_PATH)}

exec_mode = "default"
if( exists("MY_TERMINAL") ){ exec_mode = "custom" }

if (exec_mode eq "custom"){ load sprintf("%s/terminal_custom.gp", TMP_PATH) }

if(exec_mode eq "default"){ load sprintf("%s/terminal_default.gp", TMP_PATH)}
