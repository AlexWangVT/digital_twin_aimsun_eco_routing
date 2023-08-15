# Digital Twin with Aimsun

## Dependencies
- Aimsun 22.0.2
- C++ 11
- Python (for plotting)
  - pandas
  - seaborn
  - matplotlib


## Build the API target library for aimsun
- Open `AAPI.vcxproj` in Visual Studio 
- In file `AAPI.cxx`, update `PROJECT_DIR` with local directory of this repo
- Choose build type (Debug/Release) & (x64/x86). Release x64 is recommended.
    A Release build will generate a library named `AAPI_R.dll`, while a Debug build will have `AAPI_D.dll`
- Click Menu -> Build -> Build Solution

## Launch a simulation
- Open `aimsun_model/I66SubNetwork - DT/I66SubNetwork - DT.ang` with Aimsun
- In the project window, find the settings of the target scenarios under Project -> Scenarios. In the `Aimsun Next APIs` tab, disable all unknown APIs (excpet those), then add the locally built .dll library
- Save and run the simulation

## Directory architecture
- /aimsun_lib: official library from aimsun for API to link against
- /aimsun_model: aimsun models
  - /I66SubNetwork - DT: for development
  - /I66SubNetwork - DT_batch_demand50: only for batch experiments where demand is 50
  - /I66SubNetwork - DT_batch_demand100: only for batch experiments where demand is 100
  - /I66SubNetwork - DT_batch_demand150: only for batch experiments where demand is 150
  - /I66SubNetwork - DT_high_speed_100%ICE: increased highway speed limit
  - /I66SubNetwork - DT_presentation: for presentation use only
  - /I66SubNetwork - DT_original: original package from McMaster
  - /I66SubNetwork - DT_updated: updated package from McMaster

  note: carefully check the active API for every scenario is the expected one

- /data_history: historical data from a simulation (not used)
- /include: official header files from aimsun
- /ml_models: machine learning models (currently not supportted in aimsun)
- /plot_scripts: scripts to plot result figures
  - /figures: result figures
- /pre_built_library: some pre-built API libraries for tests or batch experiments
- /scripts: back up of useful scripts of aimsun scripts model
- /sim_logs: all experiment results
  - summary.csv.txt: column names of the summary file
  - summary_demand50.csv/summary_demand100.csv/summary_demand150.csv: batch results for every demand level
  - /batch_summary
    - /20230809: the first batch results
    - /20230817: the second batch results (adding MPG and MPGe)
- /src_backup: backup of important API implementations

## Conduct a batch experiment
Ther are 5 variables that will decide an experiment: Demand_Percentage(%),	Prediction_Horizon(min),	CAV_Penetration(%),	Eco_Routing_with_Travel_Time(0/1), and	Vehicle_Fleet. We need to manually set these parameters. 

Each Aimsun Scenario has 11 replications, which are for 11 different CAV_Penetration.

Vehicle_Fleet need to be manualy changed in the script `script_Modify_CAV_percentage_all_scenario_based_on_scenario_name` around line 15.

Demand_Percentage(%),	Prediction_Horizon(min) and	Eco_Routing_with_Travel_Time(0/1) of a scenario can be quickly set by the scenario name. This will be done by script `script_Modify_CAV_percentage_all_scenario_based_on_scenario_name`, which is already in the simulation package. For example, `DT_Mixed_Fleet_D100_P5` represent a 100% demand, 5min prediction horizon and eco-routing scenario; `DT_Mixed_Fleet_Travel_Time_D50_P10` represent a 50% demand, 10min prediction horizon and travel time routing scenario. (`Travel_Time` needs to be explicitly shown in the scenario name to represent the travel time routing). After setting the correct scenario name, right-click the target scenario and choose `scripts->script_Modify_CAV_percentage_all_scenario_based_on_scenario_name`, this will automaticly extract parameters from the scenario name and populate all needed fields and attributes to the experiments under the target scenario.

After setting all experiment variables, you can launch all experiments (now there are 11 in each scenario for 11 different CAV penetration) in batch mode in one scenario by right-clicking the scenario and choose `scripts->script_Execute_replications_in_the_target_scenario`. You may also use `script_Execute_all_replications_in_all_scenario` to launch multiple scenarios by providing scenario ids.

Note: currently, there are some memory issues in the simulation package, need to restart the Aimsun after about 50 replications, otherwise the system may crash due to memory leak.


## Plot figures
- 