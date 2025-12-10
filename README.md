# WEC-AUV Power Modeling  

## Overview  

This MATLAB codebase simulates and optimizes power dynamics between Wave Energy Converters (WECs), fleets of Autonomous Underwater Vehicles (AUVs), and a central energy storage unit. It is designed to help users analyze interdependencies between WEC sizing, AUV mission energy requirements, and available wave energy resources. The model simulates power exchanges between system components to provide recommended configurations, including the optimal number of AUVs that can be sustainably supported and the minimum required central battery capacity, given user-defined hardware and wave resource specifications.

The codebase is object-oriented, with a central script (`powerModel.m`) coordinating simulation logic and interacting with WEC, AUV, energy storage, and model output component classes. The simulation tracks battery levels, schedules AUV missions, and generates performance metrics to support system design and optimization.

## Features  

- **Flexible Simulation:** Supports various wave resource data types and AUV models.  
- **Object-Oriented Design:** Classes for WECs, AUVs, and energy storage encapsulate hardware properties and behaviors.  
- **Battery Tracking:** Monitors battery levels for the WEC, central battery, and each AUV over time.  
- **Mission Scheduling:** Generates AUV deployment schedules that prioritize battery health and maximize mission time.  
- **Optimization Outputs:** Calculates fleet size, central battery capacity, and provides example mission schedules.  
- **Performance Metrics:** Quantifies and compares performance of different systems.  

## Scientific Context  

This tool was developed as part of research into sustainable autonomous marine operations, as described in the paper "Development of a WEC - AUV power tracking model" (2026). The model facillitates the investigation of interdependencies between WEC power generation, AUV energy requirements, and hardware configurations, and is extensible to other microgrid optimization problems.  

## Folder Structure  

```
  wec-auv/
├── powerModel.m              # Main simulation script (central workflow)
├── Component Classes/        # MATLAB classes for WEC, AUV, EnergyStorage, etc.
├── Functions/                # Utility and helper functions
├── inputData/                # Wave resource data files
├── outputData/               # Folder for storing model outputs
```

## Getting Started

1. **Requirements:**  
   - MATLAB R2021a or newer recommended.  
   - Folders must be in the MATLAB path.  

2. **Input Data:**  
   - Save wave resource data files (e.g., `data.mat`) in the `inputData` folder.  

3. **Edit Component Specifications (optional)**  
   - Open component class file(s) (`WEC.m`, `AUV.m`, or `EnergyStorage.m`) in MATLAB.  
   - Adjust component parameters as desired.  

4. **Run the Simulation:**  
   - Open `powerModel.m` in MATLAB.  
   - Adjust simulation parameters (e.g., simulation length, timestep) as needed.  
   - Run the script to start the simulation.  

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
    
- **Wave Resource or Power Generation Data**  

    1. Sea state: Single value 1-10.  
  
    - Uses internal power generation time series data, modeled for a two-body floating point absorber WEC in different sea states.  
    
    2. Power generation, time series  
  
    - `*.mat` file containing (nx1) vectors of time and power generation  
    
    3. Power generation, mean value  
  
    - User-input mean value of WEC power generation  
    
    4. Wave specifications, time series  
  
    - `*.csv` file containing a table with time series of significant wave height (SignificantWaveHeight), wave energy period (EnergyPeriod), peak period (PeakPeriod), and time (formatted into Year, Month, Day, and Hour columns)  
    
    5. Wave specifications, mean values  
  
    - User-input mean values of significant wave height, wave energy period, and peak period  
    
    6. Power matrix with time series  
  
    - `*.mat` file containing power matrix of energy generation for set values of significant wave height and energy period. The first column of the power matrix must be a vector of increasing significant wave height, and the first row must be a vector of increasing wave energy period. File must also contain vectors of time, significant wave height, and wave energy period.   

- **Comparison Variable**  

    1. AUV Model: Runs simulation for each AUV model in list  
      
    2. WEC Power Generation: Runs simulation for each input dataset provided  
      
- **Simulation Properties**  
    - Maximum AUV fleet size: Enforces a maximum AUV fleet size, or set to `0` to disable  
    - Enforce AUV deployment stagger: If there are multiple AUVs in the fleet, enforcing a stagger (`1`) prevents grouped AUV deployments, and staggers deployments proportionally. Set to `0` to disable.  

- **Component Properties**  
  Component property edits can be made within the class files or to the generated class objects  
    - WEC Characteristic Dimension  
    - WEC, AUV, AUV dock, and central battery hotel loads  
    - WEC, AUV, and central battery capacities  
    - Efficiencies  

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
   - Re-evaluates fleet size and battery capacity.  
   - Outputs performance metrics and plots.  

[//]: # (## References)

## Contact

For questions or contributions, contact Ama Hartman: AmaH@uw.edu.
