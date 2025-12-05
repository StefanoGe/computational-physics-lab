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


create_cmake(){
	cfile="$1"
	path="$2"
	cfile="${cfile#*//}"
	touch "$path.CMakeLists.txt"
	"set(EXERCISE_SOURCES $cfile)" >> "$path.CMakeLists.txt"
	"add_executable($cfile ${EXERCISE_SOURCES})" > "$path.CMakeLists.txt"
	"target_link_libraries($cfile PRIVATE cp_library m)" > "$path.CMakeLists.txt"
	"target_include_directories($cfile PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)" > "$path.CMakeLists.txt"
}


create_subdirectories(){
	dir=$1
	for cfile in $dir/*.c; do
	[ -f "$cfile" ] || break
	cfile="${cfile#*//}"
	cfile="${cfile%.c}"
#			echo "$cfile"
	create_dirs $dir $cfile
	done
}

for lezdir in */; do
	if [[ $lezdir == LezioneTest ]]
	then
#		echo $lezdir
		create_subdirectories $lezdir
		for esdir in "$lezdir"/*/; do
			for cfile in $esdir/*.c; do
				[ -f "$cfile" ] || break
				echo "$cfile"
			done
		done
	fi
done





