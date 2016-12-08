# DeepDOM

This project contains processing scripts for the Hoarfrost & Arnosti project conducted on the DeepDOM cruise,
aimed at understanding patterns in organic carbon cycling across latitude and depth in the South Atlantic Ocean.

It also contains several reference text files necessary for running scripts.

To generate the hydrolysis rates and figures from the associated publication, clone this github repository and run the shell script:

`sh do_DeepDOM.sh`

The above shell script will run the five scripts:
1. DownloadData_DeepDOM.sh

  Will download the data inputs for this project:
  * csvs-for-rates/: gel permeation chromatography (GPC) data of substrate fluorescence
  * MonResultsMaster_DeepDOM.csv - contains monomeric substrate hydrolysis rate results, calculated manually

2. FlaRates_DeepDOM.R

  Calculates rates from slant-corrected data input in folder csvs-for-rates.

  **Inputs**:
  * csvs-for-rates folder and its 3292 .csv files
  * stdbins/ folder containing 11 stdbins* files
  * HydrolysisCutsInfo_DeepDOM.csv
  * FLAElapsedTime_DeepDOM.csv
  * FLAMasterEmpty_DeepDOM.csv
  * StdsForRates_DeepDOM.csv

  **Outputs**: FLAMasterRates_DeepDOM.csv, with calculated rates.   

3. FlaRatesSpecialCases.R

  Calculates new rates for particular incubations tagged for manual adjustment, "change to zero" and "ignore FLA".

  **Inputs**: Takes same inputs as FlaRates_DeepDOM.R.

  **Outputs**: Updated FLAMasterRates_DeepDOM.csv, with tagged rows adjusted.

4. FurtherProcessing_DeepDOM.R

  Takes final calculated FLA rates, adds factor labels for metadata and finds max rate timepoint. Also finds timepoints of max activity from monomeric hydrolysis rates.

  **Inputs**:
  * FLAMasterRates_DeepDOM.csv
  * FLAElapsedTime_DeepDOM.csv
  * MonResultsMaster_DeepDOM.csv

  **Outputs**: FlaRatesWithFactors_DeepDOM.csv, FlaMaxRates_DeepDOM.csv, and MonMaxRates_DeepDOM.csv. FlaRatesWithFactors contains all of the timepoints measured, the FlaMaxes contains only the timepoint that corresponded with maximum activity, and MonMaxRates contains timepoints of max activity and associated metadata.       

5. figures_DeepDOM.R

  Does any statistical analyses required and generates figures in publication in "figures" folder.

  Inputs:
  * ChemData_DeepDOM.csv
  * FlaMaxRates_DeepDOM.csv
  * MonMaxRates_DeepDOM.csv

  Outputs: figures 2-6 from the publication. Code for B/W versions of the figures are adjacent to each plotting function, but are commented out.


## Reference text files:

This repository contains the following reference text files containing information necessary for calculating hydrolysis rates:

* **ChemData_DeepDOM.csv**; *contains CTD data from associated cast, and nutrient from closest depth measured, for each sampling site*
* **stdbins/** folder containing 11 files beginning stdbins*; *contains info about what time standards of known molecular weight come off the GPC column*
* **HydrolysisCutsInfo_DeepDOM.csv**; *contains info about how many hydrolysis cuts need to be made to produce byproducts of standard molecular weights*
* **FLAElapsedTime_DeepDOM.csv**; *contains incubation sampling times, for rate calculations*
* **FLAMasterEmpty_DeepDOM.csv**; *contains manual tags of particular samples which needed special rate considerations after manual inspection of chromatograms*
* **StdsForRates_DeepDOM.csv**; *contains info about which standard sets to use for each sample, depending on when they were run on the GPC*
