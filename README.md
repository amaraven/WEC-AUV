# WEC-AUV Power Modeling  
Version 2.0.0

## Overview  

This MATLAB codebase simulates power dynamics between Wave Energy Converters (WECs), fleets of Autonomous Underwater Vehicles (AUVs), and a central energy storage unit. It is designed to help users analyze interdependencies between WEC sizing, AUV mission energy requirements, and available wave energy resources. The model simulates power exchanges between system components to provide recommended configurations, including the optimal number of AUVs that can be sustainably supported and the minimum required central battery capacity, given user-defined hardware and wave resource specifications.

The codebase is object-oriented and GUI-driven. Users configure simulation parameters, hardware specifications, power resource inputs, and execute simulation batches through a graphical interface (powerModelGUI). Internally, the GUI organizes user-inputs into class objects, prompts pre- and post-simulation calculations through the function runConfig.m, which also calls the function simulateSystemOps.m to execute timestep calculations. At each timestep, power exchanges between the WEC, central battery, and AUV(s) are resolved and AUV operational states are updated according to a set of scheduling priorities and energy availability. 

Model results are saved in the outputData/ folder for post-processing.

## Features  

- **Flexible Simulation:** Supports various wave resource data types and AUV models.  
- **Object-Oriented Design:** Classes for WECs, AUVs, and energy storage encapsulate hardware properties and behaviors.  
- **Battery Tracking:** Monitors battery levels for the WEC, central battery, and each AUV over time.  
- **Mission Scheduling:** Generates AUV deployment schedules that prioritize battery health and maximize mission time.  
- **System Configuration Recommendations:** Calculates fleet size, central battery capacity, and provides example mission schedules.  
- **Performance Metrics:** Quantifies and compares performance of different systems.  

## Scientific Context  

This tool was originally developed as part of research into sustainable autonomous marine operations, as described in the paper "Development of a WEC - AUV power tracking model" (2026). The model facilitates the investigation of interdependencies between WEC power generation, AUV energy requirements, and hardware configurations, and is extensible to other microgrid optimization problems.  

## Folder Structure  

```
  wec-auv/
├── modelStartup.m          # Adds subfolders to path and launches GUI
├── Classes/                # MATLAB classes for WEC, AUV, EnergyStorage, etc.    
│   ├── AUV.m                   # AUV component class
│   ├── WEC.m                   # WEC component class
│   ├── EnergyStorage.m         # Central energy storage component class
│   ├── ModelInput.m            # Stores and validates all user-defined simulation parameters
│   ├── ModelOutput.m           # Stores simulation results and generates output plots
├── Functions/              # Utility and helper functions
│   ├── calcFleetSize.m         # Estimates the number of AUVs the system can support
│   ├── runConfig.m             # Performs pre- and post- simulation computations
│   ├── simulateSystemOps.m     # Executes timestep calculations to simulate system operation 
├── GUI/                    # Folder containing GUI script
│   ├── powerModelGUI.mlapp     # Graphical user interface for configuring and executing simulation batches
├── inputData/              # Folder containing input data files and AUV presets
├── outputData/             # Folder for storing model outputs
```

## Getting Started

1. **Requirements:**  
   - MATLAB R2021a or newer required.    

2. **Run modelStartup.m**
    - This adds all required subfolders to path and launches the model GUI

3. **Configure Inputs & Run Simulation**
    - Adjust simulation parameters and provide model inputs as needed.
    - Model outputs are saved as a single `*.mat` file in the `outputData` folder.

## Customization  

- **Simulation Length:**  
  Edit the simulation length in `powerModel.m` or through prompted user-inputs.  
- **AUV Mission Scheduling:**  
  Choose to incorporate a stagger between AUV deployments in `powerModel.m` or through prompted user-inputs.  
- **Wave Resource:**  
  Change the loaded data file or resource type in `powerModel.m` to use different wave or power generation scenarios.  
- **Output Comparisons**  
  Compare between systems experiencing different power generation profiles or employing different AUVs.  
- **AUV Models:**  
  Edit the AUV model list in `powerModel.m` and `auv.m` to add or remove AUV types. Edit `auv.m` to modify characteristics of existing AUV models.  
- **WEC Models:**  
  Edit `WEC.m` to add or modify characteristics of existing WEC models.  
- **Component Classes:**  
  Extend or modify the classes in `Component Classes` to reformat for a different microgrid application or add new features or behaviors.  


## Inputs & Outputs  

### Inputs  
- **Simulation time**  
    Single value for the desired simulation length in hours  
    
- **Simulation timestep**  
    Model computation timestep in seconds

- **AUV Fleet Size**  
    Choose to cap the AUV fleet to a set number 

- **AUV Fleet Mission Scheduling**
    Choose to incorporate a stagger between AUV deployments

- **Wave Resource or Power Generation Data**  

    1. Sea state: Single value 1-10.  
  
    - Uses internal power generation time series data, modeled for a two-body floating point absorber WEC in different sea states.  
    
    2. Power generation, time series  
  
    - `*.mat` file containing (nx1) vectors of time and power generation. 
    
    3. Power generation, mean value  
  
    - User-input value of WEC power generation.  
    
    4. Wave specifications, time series  
  
    - `*.csv` file containing a table with time series of significant wave height (SignificantWaveHeight), wave energy period (EnergyPeriod), peak period (PeakPeriod), and time (formatted into Year, Month, Day, and Hour columns) or `*.mat` file containing vectors of significant wave height, wave energy period, peak period, and time.
    
    5. Wave specifications, mean values  
  
    - User-input values of significant wave height, wave energy period, and peak period.  
    
    6. Power matrix with time series  
  
    - `*.mat` file containing power matrix of energy generation for set values of significant wave height and energy period. The first column of the power matrix must be a vector of increasing significant wave height, and the first row must be a vector of increasing wave energy period. File must also contain vectors of time, significant wave height, and wave energy period.   

- **Comparison Variable**  

    1. AUV Model: Runs simulation for each AUV model in list  
      
    2. WEC Power Generation: Runs simulation for each input dataset provided  
      
- **Simulation Properties**  
    - Maximum AUV fleet size: Enforces a maximum AUV fleet size, or set to `0` to disable  
    - Enforce AUV deployment stagger: If there are multiple AUVs in the fleet, enforcing a stagger (`1`) prevents grouped AUV deployments, and staggers deployments proportionally. Set to `0` to disable.  

- **Component Properties**  
    WEC, AUV, and Energy Storage component properties can all be configured in the GUI

- **Output Plots**
    Select plots to generate upon completion

### Outputs
- Fleet size (maximum number of AUVs that can be supported by the system)  
- Central battery size recommendation  
- Example AUV mission schedule  
- Battery level tracking for all components  
- Global comparison plots for system performance  

## Model Logic

1. **Pre-Simulation:**  
   - Calculates power generation from wave resource data.  
   - Determines optimal fleet size and central battery capacity.  

2. **Simulation:**  
   - Tracks battery levels and schedules AUV missions over time.  
   - Implements logic to prevent battery depletion and optimize mission time.  

3. **Post-Simulation:**  
   - Quality check and outputs performance metrics and plots.  

[//]: # (## References)

## Contact

For questions or contributions, contact Ama Hartman: AmaH@uw.edu.