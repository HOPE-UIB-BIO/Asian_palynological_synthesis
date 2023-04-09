# Exploring spatio-temporal patterns of palynological changes in Asia during the Holocene
## Authors
Kuber P. Bhatta, Ondřej Mottl, Vivian A. Felde, Suzette G.A. Flantua, Hilary H. Birks, Xianyong Cao, Fahu Chen, John-Arvid Grytnes, Alistair W. R. Seddon, H. John B. Birks

### Corresponding author
Kuber P. Bhatta (kuber.bhatta@uib.no)

## Abstract
Historical legacies influence present-day ecosystem composition and dynamics. It is therefore important to understand the long-term dynamics of ecosystems and their properties. Analysis of ecosystem properties during the Holocene using fossil pollen assemblages provides valuable insights into past ecosystem dynamics by summarising so-called pollen-assemblage properties (PAPs). Using 205 fossil pollen data-sets (records), we quantify eight PAPs (pollen-taxonomic richness, diversity, evenness, pollen-compositional turnover, pollen-compositional change, and rate of pollen-compositional change (RoC)) for the Asian continent at different spatial scales (in individual records, within and across climate-zones, and within the continent) and time (temporal patterns over the past 12000 years). Regression tree (RT) partitioning of the PAP-estimates using sample-age as a sole predictor revealed the ‘change-point(s)’ (time or sample-age of major change in a PAP). We estimated the density of RT and multivariate regression tree (MRT) change-points in 1000-year time bins during the Holocene. 
Pollen-compositional turnover (range of sample scores along the first DCCA axis) and change (number of MRT partitions) in each record reveal gradual spatial variation across latitude and a decline with longitude eastwards. Temporally, compositional turnover declines linearly throughout the Holocene at all spatial scales. Other PAPs are heterogeneous across and within spatial scales, being more detectable at coarser scales. RT and MRT change-point density is broadly consistent in climate-zones and the continent, increasing from the early- to mid-Holocene, and mostly decrease from the mid-Holocene to the present for all PAPs.
The heterogenous patterns in PAPs across the scales of study most likely reflect responses to variations in regional environmental conditions, anthropogenic land-use, and their interactions over space and time. Patterns at the climate-zone and continental scales indicate a gradual but congruent decline in major PAPs such as compositional turnover, rate of compositional change, and major temporal compositional changes (MRT) during the Holocene, especially during recent millennia, suggesting that vegetation in Asia has become progressively more homogenous. 
Data properties (e.g., spatial distribution of the records, distribution of samples within the records, and data-standardisation and analytical approaches) may also have partly influenced the results. It is critically important to evaluate the data properties and the approaches to data standardisation and summarisation.

## How to access the repo?
1. Clone the repo.
2. Run the script for initial project setup "___Init_project___.R"
3. Make sure that you can source the configuration file "00_Config_file.R" without any error message. You need to load this config file before running each R script.
4. Install required packges outside (independent of) the "00_Config_file.R", if required to run (source) the config file without any error. Please do not forget to install "REcopol" package for estimating the pollen-assemblage properties: 
### install.packages("devtools")
### devtools::install_github("HOPE-UIB-BIO/REcopol")
