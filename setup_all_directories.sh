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
	path=
	touch "$path.CMakeLists.txt"
	set(EXERCISE_SOURCES gauss_legendre_weights.c)
	add_executable(gauss_legendre_weights ${EXERCISE_SOURCES})
	target_link_libraries(gauss_legendre_weights PRIVATE cp_library m)
	target_include_directories(gauss_legendre_weights PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)

	
	
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

for dir in */; do
	if [[ $dir == Lezione* ]]
	then
		echo $dir
		create_subdirectories $dir
	fi
done
