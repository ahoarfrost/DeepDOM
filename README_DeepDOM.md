# DeepDOM

This project contains processing scripts for the Hoarfrost & Arnosti project conducted on the DeepDOM cruise,
aimed at understanding patterns in organic carbon cycling across latitude and depth in the South Atlantic Ocean.

The following six scripts should be run in order to calculate rates
    1. DownloadData_DeepDOM.sh
        Will download the data inputs for this project:
            * csvs for rates: gel permeation chromatography (GPC) data of polysaccharide fluorescence
            * MonResultsMaster_DeepDOM.csv - contains monomeric substrate results, calculated manually
            * CTD data and nutrient data
            #* ChemData_DeepDOM.csv - Nutrient and CTD data specific to sampling locations in this study, from the corresponding CTD and nutrient data from nearest depth measured for each sampling site

    2. FlaRates_DeepDOM.R
        Calculates rates from slant-corrected data input in folder csvs-for-rates.
        requires:
            * csvs-for-rates folder and its 3292 files
            * HydrolysisCutsInfo.csv
            * FLAElapsedTime_DeepDOM.csv
            * FLAMasterList_DeepDOM.csv
            * StdsForRates_DeepDOM.csv and 11 files beginning stdbins*
        Will produce FLAMasterListFinal_DeepDOM.csv, with calculated rates, as output.   

    3. FlaRatesSpecialCases.R
        Calculates new rates for particular incubations tagged for manual adjustment,
        "change to zero" and "ignore FLA".

        Takes same inputs as FlaRates_DeepDOM.R. For those rows tagged in FLAMasterListFinal_DeepDOM.csv,
        calculates appropriate rates and produces updated FLAMasterListFinal_DeepDOM.csv as output.

    4. FurtherProcessing_DeepDOM.R
        Takes final calculated FLA rates, adds factor labels for metadata and finds max rate
        timepoint, saves two new tables FlaRatesWithFactors.csv and FlaMaxes.csv. The former
        contains all of the timepoints measured, the latter contains only the timepoint that
        corresponded with maximum activity.

        Also takes Monomeric substrate data master sheet MonResultsMaster_DeepDOM.csv, finds timepoints of max activity
        and produces MonMaxRates_DeepDOM.csv as output.

        Takes as input:
            * FLAMasterListFinal_DeepDOM.csv
            * FLAElapsedTime_DeepDOM.csv
            * MonResultsMaster_DeepDOM.csv

    5. ChemData_DeepDOM.R
        Wrangles CTD bottom casts used to collect water and record physicochem data,
        and nutrient data measured in discrete samples and downloaded from BCO-DMO.

        Creates new table recording each sampling site with its associated physicochem data
        from the corresponding CTD and nutrient data from nearest depth measured for each sampling site,
        saves output as ChemData_DeepDOM.csv.

        Takes as input:
            * FlaMaxRates_DeepDOM.csv
            * chem-data folder and its 19 files (18 ctd casts, 1 nutrient file as .cdf)            

    6. figuresDeepDOM.R
        Does any statistical analyses required and generates figures in publication in "figures" folder.
