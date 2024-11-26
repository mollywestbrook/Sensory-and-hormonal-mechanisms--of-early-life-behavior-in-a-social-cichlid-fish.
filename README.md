All R Code for this project was written by Molly Westbrook.

All R code can be found in the github repository associated with this project. Additionally, hard copies of the repository can be found on the Juntti lab fish room computer in BPS2202 (working directory below), as well as on Molly Westbrook’s external HDD.

In order to run R scripts, within each folder of the R code folder are associated project folders. Within each of these folders is an R project. As long as the script of interest and the data sheet of interest are located within the same folder as the R project (typically by temporarily moving your script and data to the project root folder), the script will run without having to worry about working directories (automatically coded into the script). 

The file structure is as follows:

Aggression Emergence Boris Data: all analysis for each experiment’s boris aggression behavior counts. Each experiment’s code, summary csvs and figures are located within a separate folder matching the file structure of the boris video archive. 

Machine-learning: the analysis for Figure 3. 
Script contents
  Scripts labelled ‘sleap_#’ are the flow of scripts used to take a raw h5 output from the sleap project and convert it into a csv used in analysis. Each script is labelled with its purpose.
  Sleap_1_organization: the initial round of organization. Converts a SLEAP h5 to a csv, interpolates NAs linearly and smooths data. Converts to mm. Outputs with a _sleap1output label. 
  Sleap_2_erroranalysis: a pdf output of the error metrics SLEAP includes in the h5.
  Sleap_3_inputborisdata: takes a cleaned boris output containing only the Observation ID, Behavior Code and the ‘Start’ timestamp column, and the output of sleap_1 and meshes them together to include a frame timestamp for when behaviors of interest are occurring. Outputs with a _sleap3output label.
  Sleap_4_featurecalculator: takes the output of sleap_3 and calculates behavioral features of interest. Outputs with a _sleap4output label.
  Sleap_5_instances: Takes the output of sleap_4 and isolates moments of behaviors, include a +/- 2 second window of data surrounding the timestamp of the behavior. Outputs aren’t standardized as these scripts are mostly in place for individual behavioral analysis. 
  Sleap_6_behaviorvisualizer: takes the output of Sleap_5 and performs a number of visualizations of the features calculated in sleap_4. Outputs figure pngs. 
  Final-thesis-figures-and-analysis: the script used for the actual figures in MW’s thesis.
A few datasets:
  Bigdf_forcorrelations: this is the final dataset used in the final-thesis-figures-and-analysis script. Contains organized sleap data for the 17 datasets of interest. 
  Totalcountofbehaviors: the summarized counts of aggressive behaviors of each video used in the above dataset. 
  SLEAP-data-deck: early analysis making sure SLEAP’s tracker matched with video data and we weren’t crazy. 
Other folders of interest:
  2024_08_finalanalysis-for-manuscript: the folder containing the datasets (all outputs of the above scripts) used in MW’s thesis in labelled folders. 
  2024_03_batchtrack-for-grand-talk: folder containing figures and data used in MW’s second year grand talk. 
  Final Manuscript Figures: the folder containing the figures used in MW’s thesis. 
  MM-conversion-calc-Reference-info: the initial photos used to perform the pixel to mm conversion for sleap_1. 
  Other scripts: old scripts made that are archived but no longer needed. Useful code may be extracted from them, but none necessary for the analysis of this project. 
  WT031623_26dpf: a significant amount of analysis was performed on this dataset to troubleshoot this project. This folder has been left containing that information, as well as all of the datasets associated with it (h5, output of the sleap organization scripts, boris data).
  
Sex Steroid Dosing Sex Ratios: the data and analysis used for Figure 4b. 
