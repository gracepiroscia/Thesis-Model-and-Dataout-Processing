# Thesis Modelling and Data Processing
A description of the files used to achieve the relevant methods described in the thesis report is outlined below.

### Model Creation
- Script "ModelAnalysis.m" contains the model functions used, produces plots to analyse the relationship between model input parameters, and produces the struct "dict.m" (containing 56 trajectories, their input parameters and other various markers that can be measured using the model)

### Model Validation
The scripts used to perform model validation should be executed in this order:
1. "ProcessBag.m": processes the ros bags from the [USyd Campus Dataset](https://dx.doi.org/10.21227/sk74-7419) to extract the linear X, Y, Z acceleration, linear X, Y, Z velocity, and corresponding timestamp information of the electric vehicle, storing the information in a MATLAB structure.
2. "TimeIntervalsOfInterest.m": processes the MATLAB structures (one-by-one) produced by "ProcessBag.m" and crops the data (by time) to include only the information representing a pedestrian crossing interaction. These times were determined by manually looking through the footage available in the Usyd Campus Dataset. This data was stored to a MATLAB structure.
3. "CombineIntervalsOfInterest.m": collate the structures produced in "TimeIntervalsOfInterest.m" into a single MATLAB structure "Combined_Trajectories.mat".
4. "DecelModelFitting.m": fits the velocity profiles from "Combined_Trajectories.mat" with the model, using a MATLAB optimizer. Outputs the goodness-of-fit (R^2) for each of the trajectories, as well as general statistics for how the model performed as a whole.

### Chosen Trajectories
- Script "TrajectoriesToModel.m" plots the five trajectories (velocity with time, velocity with distance, and distance with frames for animation programming) chosen to simulate for the thesis study (justification found in paper).

### Measures
The scripts used to process the participant data are outlined below and should be executed in the given order. All of these scripts are within the "Process_DataOut" directory
1. "csv_import.m": process the raw .csv files collected during participant trials for both the objective (unreal engine measurements) and subjective (Likert-Scale questionnaire) data, collating them into a single struct "data_out.mat".  There is the option to visualise this data per participant and press any key to move to the next participant plots, as controlled by setting the variable "plot_data" to true or false.
2. "process.m": script processes "data_out.mat" for the objective data, producing the markers used in the thesis (average distance-to-collision-zone, etc.).
  
### Analaysis of results
All of the results, including plots and p-values, given in the thesis paper, can be reproduced by running "analysis.m".
