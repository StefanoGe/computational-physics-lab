# Computational Physics Laboratory (C)

<!-- Short description of the repository -->
A collection of exercises and small programs for a Computational Physics Laboratory course.
The code is written in **C23** and built using **CMake**.

> **Development status**
>
> This repository is **work in progress**.
> Most programs have been implemented and are functional,
> but the accompanying **report is still being written**
> and most of the codebase is planned to be **refactored**.

## Repository structure

- `libs/` — shared library code
- `programs/` — executables / lab exercises
- `CMakeLists.txt` — top-level CMake configuration file
- `report/` ­— final report (WIP)

## Requirements

- C compiler supporting C23
- CMake **version 3.10 or newer**
- Make or Ninja (on Linux/macOS)

## Installation of requirements

### Ubuntu / Debian
```bash
sudo apt update
sudo apt install -y build-essential cmake
````

## Downloading the repository

### Option 1: Clone using Git (recommended)

```bash
git clone <REPOSITORY_URL>
cd <REPOSITORY_FOLDER>
```

### Option 2: Download ZIP

1. Click **Code → Download ZIP** on GitHub
2. Unzip the archive
3. Open a terminal inside the extracted folder


## Building the project


From the repository root directory:

```bash
mkdir build
cmake -S . -B build
cmake --build build
```

## Running the programs


After building, executables are located inside the `build/` directory
(or in custom directories if specified via `RUNTIME_OUTPUT_DIRECTORY`).

Example:

```bash
./build/programs/<executable_name>
```

<!-- Help users discover executables -->

To list available build targets:

```bash
cmake --build build --target help
```

