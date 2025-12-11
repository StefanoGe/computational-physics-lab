#!/usr/bin/bash

create_dirs(){
	path="./$1/$2"
	echo "$path"
	mkdir -p "$path"
	mkdir -p "$path/data"
	mkdir -p "$path/plot"
	mkdir -p "$path/script"
	mv "$path.c" $path
}

update_dirs(){
	path="$1"
	mkdir -p "$path/data"
	mkdir -p "$path/plot"
	mkdir -p "$path/script"
}


create_cmake(){
	cfile="$1"
	path="$2"
	betterpath="${path/\/\//\/}"
	betterpath="${betterpath%/}"
	cfile="${cfile#*//}"
	cfile="${cfile#*//}"
	cfile="${cfile%.c}"
	touch "$path/CMakeLists.txt"
	echo "set(EXERCISE_SOURCES $cfile.c)" > "$path/CMakeLists.txt"
	echo "add_executable($cfile \${EXERCISE_SOURCES})" >> "$path/CMakeLists.txt"
	echo "target_link_libraries($cfile PRIVATE cp_library m)" >> "$path/CMakeLists.txt"
	echo "target_include_directories($cfile PRIVATE \${CMAKE_SOURCE_DIR}/Libraries/header)" >> "$path/CMakeLists.txt"
	echo "set_target_properties($cfile PROPERTIES RUNTIME_OUTPUT_DIRECTORY \"\${CMAKE_SOURCE_DIR}/$betterpath\")" >> "$path/CMakeLists.txt"
	echo "done"
	echo "$path/CMakeLists.txt"
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
	if [[ $lezdir == Lezione* ]]
	then
#		echo $lezdir
		create_subdirectories $lezdir
		for esdir in "$lezdir"/*/; do
			for cfile in $esdir/*.c; do
				[ -f "$cfile" ] || break
#				echo "$cfile"
				create_cmake $cfile $esdir
				update_dirs $esdir
				s="$esdir"
				s="${s/\/\//\/}"
				echo "$s"
				
				if grep -q "$s" "CMakeLists.txt"; then
					echo "is_there"
				else
					echo "$s is not there"
					echo "add_subdirectory($s)" >> "CMakeLists.txt"
				fi
			
			done
			
		done
	
		
	fi
done




