# Digital Twin with Aimsun


## How to build the target
- Open `AAPI.vcxproj` in Visual Studio 
- In file `AAPI.cxx`, update `PROJECT_DIR` with local directory of this repo
- Choose build type (Debug/Release) & (x64/x86). Release x64 is recommended.
    A Release build will generate a library named `AAPI_R.dll`, while a Debug build will have `AAPI_D.dll`
- Click Menu -> Build -> Build Solution

## How to launch a simulation
- Open `aimsun_model/I66SubNetwork_aramco.ang` with Aimsun
- In the project window, find the settings of the target scenarios under Project -> Scenarios. In the `Aimsun Next APIs` tab, 
  disable all unknown APIs, then add the locally built .dll library
- Save and run the simulation

