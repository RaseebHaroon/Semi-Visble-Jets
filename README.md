# Semi-Visble-Jets

This analysis is a key part of two important academic endeavors: a master's major project named "Prospects of New Physics with Jet Observables" and a conference paper titled "A Simulation Study of Soft and Hard Radiation using Jets at the LHC." We presented our findings at the 6^th^ National Conference On Advanced Materials And Radiation Physics, hosted by the Department of Physics at the Sant Longowal Institute Of Engineering And Technology in Sangrur, Punjab, India.

- `MC_Darkmatter.cc`: This is the main analysis code written for the project.
- `input_file.hepmc`: Name of the generated HEPMC file
- `input_file.yoda`: The input data file that will be analyzed.
- `output_folder/`: The folder where the analysis results and visualizations will be stored.

## Step 1: Set Up Rivet Environment

Before running the analysis, you need to set up the environment, 
   ```bash
   source path/to/rivetenv.sh
   ```
## Step 2: Analysis Setup

Clone or download this repository for the analysis:
   ```bash
   git clone <repository_url>
   cd <repository_directory>
   ```
   
## Step 3: Build the Analysis Code

To build the analysis code:

   ```bash
   rivet-build RivetMC_Sub.so MC_Darkmatter.cc
   ```
   
## Step 4: Running the Analysis

Once you've set up the environment and built the analysis code, you can run the analysis using the following command:

   ```bash
      rivet -a MC_Darkmatter input_file.hepmc --pwd -o input_file.yoda --skip-weights
   ```

## Step 5: Generate HTML output for visualization

   ```bash
      rivet-mkhtml input_file.yoda -o output_folder
   ```
   
## Prerequisites

Before running the analysis, ensure you have the following dependencies installed:

1. **For HEPMC generation**
   - Download  [Madgraph](https://launchpad.net/mg5amcnlo) and [Pythia](https://pythia.org/), for HEPMC generation.

2. **Rivet**
   - Install Rivet by following the instructions in the [Rivet documentation](https://rivet.hepforge.org/).

3. **Semi-Visible Jets Dark Matter Model Repository**
   - Clone or download the Dark Matter Model repository from [here](https://github.com/smsharma/semi-visible-jets).

If you have any questions or encounter issues with the analysis code or repository, please feel free to open an issue or contact Raseeb Firdous Haroon at raseeb.haroon1@s.amity.edu .

## Happy analyzing! 🤗
