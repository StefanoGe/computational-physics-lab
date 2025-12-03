#!/usr/bin/bash

read -p "Enter parent directory name: " parent_dir

go_on=y

while [[ $go_on == y || $go_on == Y ]]; do
	
	read -p "Enter child directory: " child_dir
	path="./$parent_dir/$child_dir"
	mkdir -p "./$parent_dir"
	mkdir -p "$path"
	mkdir -p "$path/data"
	mkdir -p "$path/plots"
	mkdir -p "$path/conf"
#	cp --update=none "./template.gp" "$path/plots/template.gp"
	read -p "Go on?[y/n]: " go_on
done
