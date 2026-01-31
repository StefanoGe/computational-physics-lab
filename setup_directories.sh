#!/usr/bin/env bash
set -euo pipefail

PROGRAMS_DIR="programs"

create_dirs() {
  local macro_dir="$1"
  local exercise_name="$2"
  local exercise_dir="${macro_dir%/}/${exercise_name}"

  mkdir -p "${exercise_dir}/data" "${exercise_dir}/plot" "${exercise_dir}/script"
  echo "$exercise_dir"
}

update_dirs() {
  local exercise_dir="$1"
  mkdir -p "${exercise_dir%/}/data" "${exercise_dir%/}/plot" "${exercise_dir%/}/script"
}

create_subdirectories() {
  local macro_dir="${1%/}"

  for exercise_dir in "${macro_dir}"/*/; do
    [[ -d "$exercise_dir" ]] || continue

    update_dirs "$exercise_dir"

    shopt -s nullglob
    c_files=("${exercise_dir}"/*.c)
    shopt -u nullglob

    (( ${#c_files[@]} == 0 )) && continue

    cmake_file="${exercise_dir%/}/CMakeLists.txt"
    : > "$cmake_file"

    # Make output directory absolute at configure time, based on source tree.
    # This ensures executables end up inside the exercise folder.
    rel_from_root="${exercise_dir#./}"
    rel_from_root="${rel_from_root#/}"

    for src in "${c_files[@]}"; do
      name="$(basename "$src" .c)"

      echo "add_executable(${name} ${name}.c)" >> "$cmake_file"
      echo "target_link_libraries(${name} PRIVATE cp)" >> "$cmake_file"
      echo "set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY \"\${CMAKE_SOURCE_DIR}/${rel_from_root%/}\")" >> "$cmake_file"
      echo >> "$cmake_file"
    done
  done
}

main() {
  if [[ ! -d "$PROGRAMS_DIR" ]]; then
    echo "Error: '${PROGRAMS_DIR}/' not found (run from repo root)."
    exit 1
  fi

  for macro_dir in "${PROGRAMS_DIR}"/*/; do
    [[ -d "$macro_dir" ]] || continue
    create_subdirectories "$macro_dir"
  done

  echo "Done."
}

main "$@"


