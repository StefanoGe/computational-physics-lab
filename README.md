# Computational Physics Laboratory (C)

<!-- Short description of the repository -->
A collection of exercises and small programs, written by Stefano Ge for the Computational Physics Laboratory course at University of Milano-Bicocca, during 2025/2026 winter semester.
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

## Downloading the repository

### Option 1: Clone using Git (recommended)

```bash
git clone https://github.com/StefanoGe/computational-physics-lab.git
cd computational-physics-lab
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


After building, executables are located in the same foder as their corresponding .c source.

Example, if a file is in:

```bash
./build/programs/error_analysis/Es1_1/exponential_approximation
```
Then, to execute:

```bash
cd build/programs/error_analysis/Es1_1/
./exponential_approximation
```
