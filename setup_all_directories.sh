#!/usr/bin/bash

create_dirs(){
	path="./$1/$2"
	echo "$path"
	mkdir -p "$path"
	mkdir -p "$path/data"
	mkdir -p "$path/plots"
	mkdir -p "$path/conf"
	mv "$path.c" $path
}


for dir in */; do
	if [[ $dir == Lezione* ]]
	then
		echo $dir
		for cfile in $dir/*.c; do
			[ -f "$cfile" ] || break
			cfile="${cfile#*//}"
			cfile="${cfile%.c}"
#			echo "$cfile"
			create_dirs $dir $cfile
		done
	fi
done
