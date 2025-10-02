# `Compressible_Cavity`: Compressible Flow Solver for Lid-Driven Cavity Flow

This repository contains the source code for Mini-Project 2 for MAE557: Simulation and Modeling of Fluid Flows. The project consists of a C++ based numerical solver for the two-dimensional, unsteady, compressible Navier-Stokes equations in a lid-driven cavity configuration.

The solver uses a conservative finite difference method with a Forward Euler time integration scheme. Output data is saved in the HDF5 format.

## Requirements

To compile and run this project, you will need the following dependencies installed:

*   **C++ Compiler:** A modern C++ compiler that supports OpenMP (e.g., `g++`).
*   **HDF5 Library (with C++ bindings):** Required for reading and writing simulation data files.
*   **Python 3:** Required for running the analysis and plotting scripts.
*   **Python Libraries:** `numpy`, `matplotlib`, `h5py`, `tqdm`.
    *   These can be installed via pip: `pip install numpy matplotlib h5py tqdm`

## Code Structure

The project is organized into the following directories and files:

```
.
├── data/                     # Recommended directory for simulation output (.h5 files)
├── compareplot.py            # Script for Mach number comparison plots
├── convergence.py            # Script for spatial and temporal convergence plots
│
└── src/                      # Contains all C++ source code and the general plotter
    ├── plotter.py            # General-purpose plotting script for single simulations
    ├── main.cpp              # Main C++ driver, sets parameters and runs the simulation loop
    ├── grid.hpp              # Header for the SimpleMesh class, data structures, and I/O
    └── hydro.hpp             # Header for the physics implementation (BCs, timestep, etc.)
```

## How to Compile and Run

The workflow involves configuring simulation parameters directly in the source code, compiling, running the simulation, and then plotting the results.

### 1. Configure Simulation Parameters

Open `src/main.cpp` in a text editor. The primary simulation parameters are defined in the variables at the top of the `main` function. Adjust these values to set up your desired case.

| Variable | Description | Example Value |
| :--- | :--- | :--- |
| `nx1`, `nx2` | The number of grid points in the x and y directions. | `64` |
| `ng` | The number of ghost cells for boundary conditions. | `1` |
| `gamma` | The ratio of specific heats for the fluid. | `1.4` |
| `R` | The specific gas constant. | `287` |
| `CFL` | The Courant-Friedrichs-Lewy number for timestep stability. | `0.01` |
| `printfreq` | How often (in iterations) to print status to the console. | `100` |
| `Re` | The Reynolds number of the flow. | `100` |
| `T0` | The initial and reference temperature in Kelvin. | `300` |
| `Ma` | The Mach number of the flow. | `0.025` |
| `omegat_max`| The final simulation time in non-dimensional units ($\omega t$). | `1.0` |
| `savedt` | The time interval (in $\omega t$) between saving output files. | `0.1` |
| `wall_cond`| The thermal boundary condition for the walls. | `"adiabatic"` or `"isothermal"` |

### 2. Compile the Code

From the project's **root directory**, run the following command:

```bash
g++ -O3 -fopenmp -I/path/to/hdf5/include -L/path/to/hdf5/lib -lhdf5_cpp -lhdf5 src/main.cpp -o fluid
```

**Important:** You must replace `/path/to/hdf5/include` and `/path/to/hdf5/lib` with the actual paths to your HDF5 installation. For example, on macOS with Homebrew:

```bash
g++ -O3 -fopenmp -I/opt/homebrew/include -L/opt/homebrew/lib -lhdf5_cpp -lhdf5 src/main.cpp -o fluid
```

### 3. Run the Simulation

The executable will save `.h5` files in the directory from which it is run. It is highly recommended to organize your runs inside the `data/` directory.

```bash
# Create a directory for your run
mkdir -p data/my_first_run

# Navigate into that directory
cd data/my_first_run

# Run the executable (use ../ to point to the root directory)
../../fluid
```

## Analysis and Plotting

The repository includes three Python scripts for visualizing results.

### 1. General Plotting (`src/plotter.py`)

This script generates PNG images for each output file from a single simulation run.

**Usage:**
Run the script from the root directory, pointing it to your data folder.

```bash
python src/plotter.py -d <directory> -f <field> -n <ngrid> --dt <save_dt>
```

| Argument | Description | Example Value |
| :--- | :--- | :--- |
| `-d` | Directory containing the `.h5` files. | `data/my_first_run` |
| `-f` | The physical field to plot. | `VelMag` |
| `-n` | The number of grid cells (e.g., 64 for 64x64). | `64` |
| `--dt` | The non-dimensional time between saves (matches `savedt`). | `0.1` |

**Available Fields:** `Density`, `VelX1`, `VelX2`, `Press`, `Temp`, `VelMag`.

### 2. Convergence Plots (`convergence.py`)

This script generates the spatial and temporal convergence plots.

**Important:** This script contains **hardcoded paths** and expects a specific directory structure inside the `data/` folder. You must place your simulation output in folders with these exact names:
*   `data/convergence_16x16/`
*   `data/convergence_32x32/`
*   `data/convergence_64x64/`
*   `data/convergence_128x128/`
*   `data/time_conv_01/`
*   `data/time_conv_005/`
*   `data/time_conv_001/`
*   `data/time_conv_0005/`

**Usage:**
Once your data is organized correctly, run the script from the root directory:
```bash
python convergence.py
```
It will generate `spatial_convergence.png` and `temporal_convergence.png`.

### 3. Mach Number Comparison Plot (`compareplot.py`)

This script generates a single figure comparing the normalized velocity magnitude for four different Mach numbers at a specific time.

**Important:** Like the convergence script, this script contains **hardcoded paths**. It expects the following directories to exist in the root of the project:
*   `comparison_0.8/`
*   `comparison_0.5/`
*   `comparison_0.1/`
*   `comparison_0.025/`

You must place the corresponding simulation output in these folders. You can also edit the script to change the timestep (`Time = 9`) or the file paths.

**Usage:**
After organizing the data, run the script from the root directory:
```bash
python compareplot.py
```
It will generate `Ma_comparison.png`.
