# Computational Physics Laboratory (C)

<!-- Short description of the repository -->
A collection of exercises and small programs for a Computational Physics Laboratory course.
The code is written in **C** and built using **CMake**.

---

## Repository structure

<!-- Explain the main folders so users can orient themselves -->
- `libs/` — shared library code
- `programs/` — executables / lab exercises
- `CMakeLists.txt` — top-level CMake configuration file
- `report/` ­— final report (WIP)

<!-- Mention build directory policy -->
> ⚠️ Build files are generated in a separate `build/` directory (out-of-source build).

---

## Requirements

<!-- List software dependencies -->
- C compiler (`gcc` or `clang`)
- CMake **version 3.10 or newer**
- Make or Ninja (on Linux/macOS)

---

## Installation of requirements

<!-- Example installation commands; users may adapt them -->
### Ubuntu / Debian
```bash
sudo apt update
sudo apt install -y build-essential cmake
````

### macOS (Homebrew)

```bash
brew install cmake
```

---

## Downloading the repository

<!-- Two standard ways of obtaining the code -->

### Option 1: Clone using Git (recommended)

```bash
git clone <REPOSITORY_URL>
cd <REPOSITORY_FOLDER>
```

### Option 2: Download ZIP

1. Click **Code → Download ZIP** on GitHub
2. Unzip the archive
3. Open a terminal inside the extracted folder

---

## Building the project

<!-- Out-of-source build instructions -->

From the repository root directory:

```bash
cmake -S . -B build
cmake --build build
```

<!-- Optional Ninja build -->

### Using Ninja (optional)

```bash
cmake -S . -B build -G Ninja
cmake --build build
```

---

## Running the programs

<!-- Explain where executables end up -->

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

To find all executables:

```bash
find . -type f -executable
```

---


