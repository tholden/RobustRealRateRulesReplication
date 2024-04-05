---
author:
- Tom D. Holden
---

# README document for the replication code for the paper "Robust Real Rate Rules" by [Tom D. Holden](https://www.tholden.org/), to be published in Econometrica

**Tom D. Holden**

*The views expressed in this document are those of the author and do not represent the views of the Deutsche Bundesbank, the Eurosystem or its staﬀ.*

This document and the corresponding code and data originate from [the GitHub repository linked here](https://github.com/tholden/RobustRealRateRulesReplication).

The [working paper version of this paper](https://www.tholden.org/assets/files/RobustRealRateRulesBody.pdf), its [online appendices](https://www.tholden.org/assets/files/RobustRealRateRulesOnlineAppendices.pdf) and [supplemental appendices](https://www.tholden.org/assets/files/RobustRealRateRulesSupplementalAppendices.pdf) are available from [this link](https://www.tholden.org/#robust-real-rate-rules).

## Overview

The code in this replication package produces the numerical results, tables and figures from the paper "Robust Real Rate Rules" by [Tom D. Holden](https://www.tholden.org/) and from its [supplemental appendices](https://zenodo.org/records/10037239) (Holden 2024). The `Inputs` folder contains the source data, from [FRED](https://fred.stlouisfed.org/), [ALFRED](https://alfred.stlouisfed.org/), [the Survey of Professional Forecasters](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/median-forecasts), [Bauer & Swanson (2023)](https://www.michaeldbauer.com/files/FOMC_Bauer_Swanson.xlsx) and [Känzig (2021)](https://github.com/dkaenzig/oilsupplynews). All processing is performed by the MATLAB script `Main.m`, which saves all results in the `Outputs` folder, and which completes in around 90 seconds on a machine with 28 cores. Completion time may be proportionally longer on machines with fewer cores.

## Data Availability and Provenance Statements

All source data is contained in the `Inputs` folder.

Automatically converted CSV versions of the source data files are contained in the `InputsAsCSV` folder in the release. These are not contained in the source GitHub repository as they are generated as part of the release build script. Where the source data file had multiple sheets, these are saved in separate CSV files within one common folder.

### Statement about Rights

- [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript. 
- [x] I certify that the author(s) of the manuscript have permission to redistribute/publish the data contained within this replication package, as **all data is in the public domain**.

### Summary of Availability

- [x] All data **are** publicly available.

### Details on each Data Source

The files in the `Inputs` folder are listed below, along with their sources. Each file name links to its source webpage. Note that since the FRED and ALFRED files were downloaded using the FRED/ALFRED API, these files may differ in format from those that can be obtained directly from the source website. However, their constituent data is identical.

All data is in the public domain. No explicit copyright statements or restrictive licenses are given for any of these sources.

* From [FRED](https://fred.stlouisfed.org/) (Federal Reserve Bank of St. Louis):
  * [`DGS5.xlsx`](https://fred.stlouisfed.org/series/DGS5) (Board of Governors of the Federal Reserve System (US), Market Yield on U.S. Treasury Securities at 5-Year Constant Maturity, Quoted on an Investment Basis)
  * [`CPIAUCSL.xlsx`](https://fred.stlouisfed.org/series/CPIAUCSL) (U.S. Bureau of Labor Statistics, Consumer Price Index for All Urban Consumers: All Items in U.S. City Average.)
  * [`GS5.xlsx`](https://fred.stlouisfed.org/series/GS5) (Board of Governors of the Federal Reserve System (US), Market Yield on U.S. Treasury Securities at 5-Year Constant Maturity, Quoted on an Investment Basis.)
  * [`T5YIEM.xlsx`](https://fred.stlouisfed.org/series/T5YIEM) (Federal Reserve Bank of St. Louis, 5-Year Breakeven Inflation Rate.)
  * [`T5YIFRM.xlsx`](https://fred.stlouisfed.org/series/T5YIFRM) (Federal Reserve Bank of St. Louis, 5-Year, 5-Year Forward Inflation Expectation Rate.)
  * [`PCECTPICTMLR.xlsx`](https://fred.stlouisfed.org/series/PCECTPICTMLR) (U.S. Federal Open Market Committee and Federal Reserve Bank of St. Louis, Longer Run FOMC Summary of Economic Projections for the Personal Consumption Expenditures Inflation Rate, Central Tendency, Midpoint.)
* From [ALFRED](https://alfred.stlouisfed.org/) (Federal Reserve Bank of St. Louis): (Note, each file contains three sheets. The `Dates` sheet give the row labels for the `Data` sheet, while the `Vintages` sheet gives the column labels for the `Data` sheet.)
  * [`HistoricalCPIAUCSL.xlsx`](https://alfred.stlouisfed.org/series?seid=CPIAUCSL) (U.S. Bureau of Labor Statistics, Consumer Price Index for All Urban Consumers: All Items in U.S. City Average.)
  * [`HistoricalPCEPI.xlsx`](https://alfred.stlouisfed.org/series?seid=PCEPI) (U.S. Bureau of Economic Analysis, Personal Consumption Expenditures: Chain-type Price Index.)
  * [`HistoricalPCECTPICTM.xlsx`](https://alfred.stlouisfed.org/series?seid=PCECTPICTM) (U.S. Federal Open Market Committee and Federal Reserve Bank of St. Louis, FOMC Summary of Economic Projections for the Personal Consumption Expenditures Inflation Rate, Central Tendency, Midpoint.)
* From [the median forecasts of the Survey of Professional Forecasters](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/median-forecasts):
  * [`SPF.xlsx`]('https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/historical-data/medianlevel.xlsx') (Note, the original file name was `medianlevel.xlsx`. Note also that only the data in the sheets `CPI5YR` and `CPIF5` is used by the replication code.)
* From [the website of Michael Bauer](https://www.michaeldbauer.com/research/), with the replication data for Bauer & Swanson (2023):
  * [`FOMC_Bauer_Swanson.xlsx`](https://www.michaeldbauer.com/files/FOMC_Bauer_Swanson.xlsx)
* From [the `oilsupplynews` GitHub repository of Diego Känzig](https://github.com/dkaenzig/oilsupplynews):
  * [`OilSupplyNewsShocksKaenzig.xlsx`](https://github.com/dkaenzig/oilsupplynews/raw/master/oilSupplyNewsShocks_2023M06.xlsx) (Note, the original file name was `oilSupplyNewsShocks_2023M06.xlsx`.)

## Computational requirements

### Software Requirements

The scripts `Main.m` is intended to be run in MATLAB. It was tested with MATLAB version R2024a. Later versions should also work perfectly.

The following commercial MATLAB toolboxes from MathWorks are required:

* Optimization Toolbox
* Statistics and Machine Learning Toolbox
* Datafeed Toolbox
* Parallel Computing Toolbox
* Econometrics Toolbox

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

The code was last run on a **28-core Intel-based desktop with 256 GB of RAM, running Windows 11, with more than 1TB of free storage space**. Running `Main.m` took about 90 seconds. Completion time may be proportionally longer on machines with fewer cores.

## Description of programs/code

* The MATLAB script `Main.m` generates all needed outputs, in the `Outputs` directory.
* The functions in the folder `private` are called by `Main.m` to perform various sub-tasks.

The source GitHub repository also contains the MATLAB script `MakeRelease.m` which builds the release package based on the source files.

### License for Code

The code is provided here licensed under Version 3 of the GNU General Public License. See [LICENSE.md](LICENSE.md) for the full terms and conditions. Please contact the author should you wish to license the code under other terms.

## Instructions to Replicators

* Open MATLAB, and navigate to the folder containing this file.
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

Bauer, Michael D. & Eric T. Swanson. 2023. ‘A Reassessment of Monetary Policy Surprises and High-Frequency Identification’. *NBER Macroeconomics Annual* 37: 87–155.

Holden, Tom D. 2024. ‘Supplemental Appendices to: “Robust Real Rate Rules”’. https://zenodo.org/records/10037239

Känzig, Diego R. 2021. ‘The Macroeconomic Effects of Oil Supply News: Evidence from OPEC Announcements’. *American Economic Review* 111 (4): 1092–1125.
