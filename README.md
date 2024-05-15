---
author:
- Tom D. Holden
---

# README document for the replication code for the paper "Robust Real Rate Rules" by [Tom D. Holden](https://www.tholden.org/), to be published in Econometrica

*The views expressed in this document are those of the author and do not represent the views of the Deutsche Bundesbank, the Eurosystem or its staﬀ.*

[![DOI](https://zenodo.org/badge/780997704.svg)](https://zenodo.org/doi/10.5281/zenodo.11198015)

This document and the corresponding code and data originate from [the GitHub repository linked here](https://github.com/tholden/RobustRealRateRulesReplication).

The [working paper version of this paper](https://www.tholden.org/assets/files/RobustRealRateRulesBody.pdf), its [online appendices](https://www.tholden.org/assets/files/RobustRealRateRulesOnlineAppendices.pdf) and [supplemental appendices](https://www.tholden.org/assets/files/RobustRealRateRulesSupplementalAppendices.pdf) are available from [this link](https://www.tholden.org/#robust-real-rate-rules).

## Overview

The code in this replication package produces the numerical results, tables and figures from the paper "Robust Real Rate Rules" by [Tom D. Holden](https://www.tholden.org/) and from its [supplemental appendices](https://zenodo.org/records/10037239) (Holden 2024). The `Inputs` folder contains the source data, from [FRED](https://fred.stlouisfed.org/), [ALFRED](https://alfred.stlouisfed.org/), [the Survey of Professional Forecasters](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/median-forecasts), [Bauer & Swanson (2023a; 2023b)](https://www.michaeldbauer.com/files/FOMC_Bauer_Swanson.xlsx) and [Känzig (2021; 2024)](https://github.com/dkaenzig/oilsupplynews). All processing is performed by the MATLAB script `Main.m`, which saves all results in the `Outputs` folder, and which completes in around 100 seconds on a machine with 28 cores. Completion time may be proportionally longer on machines with fewer cores.

## Data Availability and Provenance Statements

All source data is contained in the `Inputs` folder.

### Statement about Rights

- [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript. 
- [x] I certify that the author(s) of the manuscript have permission to redistribute/publish the data contained within this replication package, as all data is in the public domain except the data in the `SmetsWouters2007` folder, which is licensed under the Creative Commons Attribution 4.0 International Public License.

### Summary of Availability

- [x] All data **are** publicly available.

### Details on each Data Source

The files in the `Inputs` folder are listed below, along with their sources. Each file name links to its source webpage. Note that since the FRED and ALFRED files were downloaded using the FRED/ALFRED API, these files may differ in format from those that can be obtained directly from the source website. However, their constituent data is identical.

All data in the `Inputs` folder is in the public domain. No explicit copyright statements or restrictive licenses are given for any of these sources.

* From [FRED](https://fred.stlouisfed.org/) (Federal Reserve Bank of St. Louis):
  * [`A794RX0Q048SBEA.xlsx`](https://fred.stlouisfed.org/series/A794RX0Q048SBEA) (U.S. Bureau of Economic Analysis 2024a).
  * [`DGS5.xlsx`](https://fred.stlouisfed.org/series/DGS5) (Board of Governors of the Federal Reserve System (US) 2024a).
  * [`CPIAUCSL.xlsx`](https://fred.stlouisfed.org/series/CPIAUCSL) (U.S. Bureau of Labor Statistics 2024).
  * [`GS5.xlsx`](https://fred.stlouisfed.org/series/GS5) (Board of Governors of the Federal Reserve System (US) 2024b).
  * [`T5YIEM.xlsx`](https://fred.stlouisfed.org/series/T5YIEM) (Federal Reserve Bank of St. Louis 2024a).
  * [`T5YIFRM.xlsx`](https://fred.stlouisfed.org/series/T5YIFRM) (Federal Reserve Bank of St. Louis 2024b).
  * [`PCECTPICTMLR.xlsx`](https://fred.stlouisfed.org/series/PCECTPICTMLR) (U.S. Federal Open Market Committee and Federal Reserve Bank of St. Louis 2024b).
* From [ALFRED](https://alfred.stlouisfed.org/) (Federal Reserve Bank of St. Louis): (Note, each file contains three sheets. The `Dates` sheet give the row labels for the `Data` sheet, while the `Vintages` sheet gives the column labels for the `Data` sheet.)
  * [`HistoricalCPIAUCSL.xlsx`](https://alfred.stlouisfed.org/series?seid=CPIAUCSL) (U.S. Bureau of Labor Statistics 2024; Federal Reserve Bank of St. Louis 2024c).
  * [`HistoricalPCEPI.xlsx`](https://alfred.stlouisfed.org/series?seid=PCEPI) (U.S. Bureau of Economic Analysis 2024b; Federal Reserve Bank of St. Louis 2024d).
  * [`HistoricalPCECTPICTM.xlsx`](https://alfred.stlouisfed.org/series?seid=PCECTPICTM) (U.S. Federal Open Market Committee and Federal Reserve Bank of St. Louis 2024a; Federal Reserve Bank of St. Louis 2024d).
* From [the median forecasts of the Survey of Professional Forecasters](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/median-forecasts):
  * [`SPF.xlsx`]('https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/historical-data/medianlevel.xlsx') (Federal Reserve Bank of Philadelphia 2024). (Note, the original file name was `medianlevel.xlsx`. Note also that only the data in the sheets `CPI5YR` and `CPIF5` is used by the replication code.)
* From [the website of Michael Bauer](https://www.michaeldbauer.com/research/), with the replication data for Bauer & Swanson (2023a):
  * [`FOMC_Bauer_Swanson.xlsx`](https://www.michaeldbauer.com/files/FOMC_Bauer_Swanson.xlsx) (Bauer & Swanson 2023b).
* From [the `oilsupplynews` GitHub repository of Diego Känzig](https://github.com/dkaenzig/oilsupplynews):
  * [`OilSupplyNewsShocksKaenzig.xlsx`](https://github.com/dkaenzig/oilsupplynews/raw/master/oilSupplyNewsShocks_2023M06.xlsx) (Känzig 2024). (Note, the original file name was different. It was `oilSupplyNewsShocks_2023M06.xlsx`.)
  
The `SmetsWouters2007` folder contains one data file: `usmodel_data.mat`, which is licensed under the Creative Commons Attribution 4.0 International Public License. This file contains the data on which the Smets & Wouters (2007) model was originally estimated, which is fully documented in the replication package, Smets & Wouters (2019). This file is unchanged from the one provided in Pfeifer (2024).

## Computational requirements

### Software Requirements

The scripts `Main.m` is intended to be run in MATLAB. It was tested with MATLAB version R2024a. Later versions should also work perfectly.

The following commercial MATLAB toolboxes from MathWorks are required:

* Optimization Toolbox
* Statistics and Machine Learning Toolbox
* Datafeed Toolbox
* Parallel Computing Toolbox
* Econometrics Toolbox

In order to replicate the Smets & Wouters (2007) results from Online Appendix C, additionally Dynare version 6.0 (Adjemian et al. 2024) is required. The `matlab` sub-folder of Dynare's install directory needs to be added to the MATLAB path, for example, by executing `addpath( 'C:\dynare\6.0\matlab' );` in MATLAB before running `Main.m`. (This example gives the correct command if you are using Windows, with the standard Dynare install location. For other platforms, see the Dynare documentation.) If Dynare is not detected, then replicating the Smets & Wouters (2007) results will be skipped.

Other than these toolboxes, there are no other requirements. The code should work on any platform on which MATLAB is available.

### Controlled Randomness

- [x] A random seed is set at line `39` of the function `private/DistributionOfCentralTendency.m`.
- [x] No other pseudo-random generation is used in the script `Main.m` or in any other function it calls.

### Memory, Runtime, Storage Requirements

#### Summary

Approximate time needed to reproduce the analyses on a standard (2024) desktop machine:

- [x] < 10 minutes

Approximate free memory needed to reproduce the analyses on a standard (2024) desktop machine:

- [x] < 8 GB

Approximate free storage space needed:

- [x] < 25 MBytes

#### Details

The code was last run on a **28-core Intel-based desktop with 256 GB of RAM, running Windows 11, with more than 1TB of free storage space**. Running `Main.m` took about 100 seconds. Completion time may be proportionally longer on machines with fewer cores.

## Description of programs/code

* The MATLAB script `Main.m` generates all needed outputs, in the `Outputs` directory.
* The folder `SmetsWouters2007` contains replication code for the Smets & Wouters (2007) model, used to generate some results in Online Appendix C. This is derived from replication materials provided by Pfeifer (2024), itself based on the authors’ original replication code (Smets & Wouters 2019). The folder contains the following files:
  * `Smets_Wouters_2007_45.mod` - The Dynare (Adjemian et al. 2024) model file for the Smets & Wouters (2007) model. A few minor changes have been made to the file from Pfeifer (2024) so that it produces the outputs needed for the paper's Online Appendix C.
  * `usmodel_data.mat` - The data on which this model was originally estimated, which is fully documented in the replication package, Smets & Wouters (2019). This file is unchanged from the one provided in Pfeifer (2024).
  * `usmodel_mode.mat` - The initial point used for searching for the mode. This file is unchanged from the one provided in Pfeifer (2024).
* The functions in the folder `private` are called by `Main.m` to perform various sub-tasks. These functions are listed below:
  * `corrNoMean.m` - Computes the correlation between two mean zero variables.
  * `DistributionOfCentralTendency.m` - Performs simulations to assess the performance of the central tendency midpoint as a measure of location, for Supplemental Appendix J.1.
  * `GCVRidgeRegression.m` - Performs ridge regression with weight chosen by generalized cross validation. Only used to provide a rough initial estimate of the time varying parameter regression in Supplemental Appendix J.2.
  * `GridFMinBound.m` - Minimises a function by performing an initial grid search then a finer search.
  * `hacIV.m` - Performs instrumental variables regression by iterated GMM, with HAC standard errors.
  * `ObtainHistoricalData.m` - Obtains multiple vintages of a FRED dataset, either from FRED or from the `Inputs` folder.
  * `ObtainRealTimeGrowthRate.m` - Obtains real time data on the growth rate of a FRED variable.
  * `PiStarParameterFunction.m` - Returns the transition matrices needed by the Kalman filter for the estimation exercise in Supplemental Appendix J.4.
  * `ReadFromFREDOrInputsFolder.m` - Obtains current data on a single series from FRED or the `Inputs` folder.
* The source GitHub repository also contains the MATLAB script `MakeRelease.m` which builds the release package based on the source files.

### License for Code

The code is provided here licensed under Version 3 of the GNU General Public License. See [LICENSE.md](LICENSE.md) for the full terms and conditions. Please contact the author should you wish to license the code under other terms.

## Instructions to Replicators

* Open MATLAB, and navigate to the folder containing this file.
* Optional: Ensure Dynare is in your MATLAB path by executing `addpath( 'C:\dynare\6.0\matlab' );` if your Dynare is installed in the folder `C:\dynare\6.0` or `addpath( 'PATH_TO_DYNARE_MATLAB_FOLDER' );` more generally, where `PATH_TO_DYNARE_MATLAB_FOLDER` is the path of the `matlab` sub-folder within your Dynare install.
* Type `Main` at the MATLAB command line, then press enter.
* This file will create all tables and figures used in the paper or its supplemental appendices, as well as all numbers directly included in the text of either document.

### Optional procedure to update the source data

To download the latest versions of the source data, replace the line `DownloadLatestData = false;` at the top of `Main.m` with the line `DownloadLatestData = true;`, then run the MATLAB script `Main.m`. The script will then redownload each of the source files from the internet, recreating the contents of the `Inputs` folder.

After performing this step, the data may be newer than the data used in the published paper, and so the results may differ.

## List of tables and programs

The provided script `Main.m` reproduces:

- [x] All numbers provided in text in the paper.
- [x] All tables and figures in the paper.

The following output files are generated by `Main.m` within the `Outputs` directory:

  * `Log.txt` - A complete log of all MATLAB command output during the script run. Check this file for numbers included in the text of the supplemental appendices or main paper.
  * `SupplementalTable1DistributionOfCentralTendencyT.xlsx` - Table 1 from the supplemental appendices.
  * `SupplementalTable2DistributionOfCentralTendencyN.xlsx` - Table 2 from the supplemental appendices.
  * `SupplementalFigure1BreakevenInflationVsSPF.svg` - Figure 1 from the supplemental appendices.
  * `SupplementalFigure2BreakevenInflationVsReality.svg` - Figure 2 from the supplemental appendices.
  * `SupplementalFigure3alpha.svg` - Figure 3 from the supplemental appendices.
  * `SupplementalFigure4beta.svg` - Figure 4 from the supplemental appendices.
  * `SupplementalFigure5rho_PCEPI_CPI.svg` - Figure 5 from the supplemental appendices.
  * `SupplementalFigure6AnnualExercise.svg` - Figure 6 from the supplemental appendices.
  * `SupplementalFigure7PiPiStar.svg` - Figure 7 from the supplemental appendices.
  * `SupplementalFigure8FiveYearExpectation.svg` - Figure 8 from the supplemental appendices.
  * `PaperFigure1QuarterlyMod.svg` - Figure 1 from the main paper.

## References

* Adjemian, Stéphane, Michel Juillard, Frédéric Karamé, Willi Mutschler, Johannes Pfeifer, Marco Ratto, Normann Rion and Sébastien Villemot (2024), “Dynare: Reference Manual, Version 6,” Dynare Working Papers, 80, CEPREMAP
* Bauer, Michael D. & Eric T. Swanson. 2023a. ‘A Reassessment of Monetary Policy Surprises and High-Frequency Identification’. *NBER Macroeconomics Annual* 37: 87–155.
* ———. 2023b. ‘High-Frequency FOMC Announcement and SVAR Data Used in the Paper: A Reassessment of Monetary Policy Surprises and High-Frequency Identification’. Retrieved from https://www.michaeldbauer.com/files/FOMC_Bauer_Swanson.xlsx, April 5 2024.
* Board of Governors of the Federal Reserve System (US). 2024a. ‘Market Yield on U.S. Treasury Securities at 5-Year Constant Maturity, Quoted on an Investment Basis \[DGS5\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/DGS5, April 5, 2024.
* ———. 2024b. ‘Market Yield on U.S. Treasury Securities at 5-Year Constant Maturity, Quoted on an Investment Basis \[GS5\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/GS5, April 5, 2024.
* Federal Reserve Bank of Philadelphia. 2024. ‘Median Forecast: Survey of Professional Forecasters’. Retrieved from https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/median-forecasts, April 5, 2024.
* Federal Reserve Bank of St. Louis. 2024a. ‘5-Year Breakeven Inflation Rate \[T5YIEM\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/T5YIE, April 5, 2024.
* ———. 2024b. ‘5-Year, 5-Year Forward Inflation Expectation Rate \[T5YIFRM\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/T5YIFRM, April 5, 2024.
* ———. 2024c. ‘Vintages of: Consumer Price Index for All Urban Consumers: All Items in U.S. City Average \[CPIAUCSL\]’. Retrieved from ALFRED, Federal Reserve Bank of St. Louis; https://alfred.stlouisfed.org/series?seid=CPIAUCSL, April 5, 2024.
* ———. 2024d. ‘Vintages of: Personal Consumption Expenditures: Chain-Type Price Index \[PCEPI\]’. Retrieved from ALFRED, Federal Reserve Bank of St. Louis; https://alfred.stlouisfed.org/series?seid=PCEPI, April 5, 2024.
* ———. 2024e. ‘Vintages of: FOMC Summary of Economic Projections for the Personal Consumption Expenditures Inflation Rate, Central Tendency, Midpoint \[PCECTPICTM\]’. Retrieved from ALFRED, Federal Reserve Bank of St. Louis; https://alfred.stlouisfed.org/series?seid=PCECTPICTM, April 5, 2024.
* Holden, Tom D. 2024. ‘Supplemental Appendices to: “Robust Real Rate Rules”’. Zenodo; https://zenodo.org/records/10037239.
* Känzig, Diego R. 2021. ‘The Macroeconomic Effects of Oil Supply News: Evidence from OPEC Announcements’. *American Economic Review* 111 (4): 1092–1125.
* ———. 2024. ‘Vintages for Oil Supply News Shocks Data’. Retrieved from https://github.com/dkaenzig/oilsupplynews, April 5, 2024.
* Pfeifer, Johannes. 2024. ‘DSGE_mod:  A Collection of Dynare Models’ (version v2.0.0). Zenodo; https://zenodo.org/records/10810290.
* Smets, Frank & Rafael Wouters. 2007. ‘Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach’. *American Economic Review* 97 (3): 586–606.
* ———. 2019. ‘Replication Data for: Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach’. Inter-university Consortium for Political and Social Research (ICPSR); https://www.openicpsr.org/openicpsr/project/116269/version/V1/view.
* U.S. Bureau of Economic Analysis. 2024a. ‘Real Personal Consumption Expenditures per Capita \[A794RX0Q048SBEA\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/A794RX0Q048SBEA, April 5, 2024.
* U.S. Bureau of Economic Analysis. 2024b. ‘Personal Consumption Expenditures: Chain-Type Price Index \[PCEPI\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/PCEPI, April 5, 2024.
* U.S. Bureau of Labor Statistics. 2024. ‘Consumer Price Index for All Urban Consumers: All Items in U.S. City Average \[CPIAUCSL\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/CPIAUCSL, April 5, 2024.
* U.S. Federal Open Market Committee and Federal Reserve Bank of St. Louis. 2024a. ‘FOMC Summary of Economic Projections for the Personal Consumption Expenditures Inflation Rate, Central Tendency, Midpoint \[PCECTPICTM\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/PCECTPICTM, April 5, 2024.
* ———. 2024b. ‘Longer Run FOMC Summary of Economic Projections for the Personal Consumption Expenditures Inflation Rate, Central Tendency, Midpoint \[PCECTPICTMLR\]’. Retrieved from FRED, Federal Reserve Bank of St. Louis; https://fred.stlouisfed.org/series/PCECTPICTMLR, April 5, 2024.
