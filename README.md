# Exploring spatio-temporal patterns of palynological changes in Asia during the Holocene
## Authors
Kuber P. Bhatta, Ondřej Mottl, Vivian A. Felde, Suzette G.A. Flantua, Hilary H. Birks, Xianyong Cao, Fahu Chen, John-Arvid Grytnes, Alistair W. R. Seddon, H. John B. Birks

### Corresponding author
Kuber P. Bhatta (kuber.bhatta@uib.no)

### Project
This publication is a part of European Research Council (ERC) under the European Union’s Horizon 2020 Research and Innovation Programme (grant agreement no. 741413) HOPE Humans on Planet Earth – longterm impacts on biosphere dynamics awarded to HJBB.

## Abstract
Historical legacies influence present-day ecosystem composition and dynamics. It is therefore important to understand the long-term dynamics of ecosystems and their properties. Analysis of ecosystem properties during the Holocene using fossil pollen assemblages provides valuable insights into past ecosystem dynamics by summarising so-called pollen-assemblage properties (PAPs). Using 205 fossil pollen data-sets (records), we quantify eight PAPs (pollen-taxonomic richness, diversity, evenness, pollen-compositional turnover, pollen-compositional change, and rate of pollen-compositional change (RoC)) for the Asian continent at different spatial scales (in individual records, within and across climate-zones, and within the continent) and time (temporal patterns over the past 12000 years). Regression tree (RT) partitioning of the PAP-estimates using sample-age as a sole predictor revealed the ‘change-point(s)’ (time or sample-age of major change in a PAP). We estimated the density of RT and multivariate regression tree (MRT) change-points in 1000-year time bins during the Holocene. 
Pollen-compositional turnover (range of sample scores along the first DCCA axis) and change (number of MRT partitions) in each record reveal gradual spatial variation across latitude and a decline with longitude eastwards. Temporally, compositional turnover declines linearly throughout the Holocene at all spatial scales. Other PAPs are heterogeneous across and within spatial scales, being more detectable at coarser scales. RT and MRT change-point density is broadly consistent in climate-zones and the continent, increasing from the early- to mid-Holocene, and mostly decrease from the mid-Holocene to the present for all PAPs.
The heterogenous patterns in PAPs across the scales of study most likely reflect responses to variations in regional environmental conditions, anthropogenic land-use, and their interactions over space and time. Patterns at the climate-zone and continental scales indicate a gradual but congruent decline in major PAPs such as compositional turnover, rate of compositional change, and major temporal compositional changes (MRT) during the Holocene, especially during recent millennia, suggesting that vegetation in Asia has become progressively more homogenous. 
Data properties (e.g., spatial distribution of the records, distribution of samples within the records, and data-standardisation and analytical approaches) may also have partly influenced the results. It is critically important to evaluate the data properties and the approaches to data standardisation and summarisation.

## How to access the repo?
"Asian_palynological_synthesis" ('R project' from here on) is accessible in two ways:
1. If a user has a [GitHub account](https://github.com/), the easiest way is to create your own GitHub repo using this Github link (https://github.com/HOPE-UIB-BIO/Asian_palynological_synthesis).
2. A user can download the latest *Release* of the R project as a zip file from the Asian_palynological_synthesis Release page (https://github.com/HOPE-UIB-BIO/Asian_palynological_synthesis/releases).

Different sections (folders) of the R project are as follows:
1. `Data/`: All the raw or pre-analysed and processed data are stored in this folder. 

1.1 The subfolder `Data/Input/` contains a geo-tif file for the biome or climate-zone classification (`Data/Input/Biomes_spatial/`), the pollen-taxa harmonisation tables (`Data/Input/Hamonisation_tables/`) used for standardising the pollen taxonomy, and the raw data used for estimating the pollem assemblage properties (PAPs). Please note that users will not find raw data here because it cannot be shared publicly due to an embargo of the data contributors. However, users can access all the PAPs estimated from the daw data in the subfolder `Data/Processed/Data_for_analyses/`.
1.2 The subfloder `Data/Processed/` contains all the processed data used for analyses of the spatio-temporal patterns. PAPs estimated from the raw data and the sample-ages (levels) are stored within the subfolder `Data/Processed/Data_for_analyses/`. Data of the estimates of individual PAPs are stored within the separate subfolders such as `Data/Processed/MVRT/`, `Data/Processed/Diversity/`..., where the name of the subfolder corresponds the abbreviated name of the PAP. Subfolder `Data/Processed/PAP_all/` contains all the PAP estimates combined together, and `Data/Processed/Metadata/` contains metadadata such as dataset ID, latitude, longitude, source of data (neotoma vs. other), number of samples in each dataset (fossil pollen record), if a dataset is pollen-counts ot percentages (TRUE/FALSE), number of chronology control points, ...

2. `Outputs/`: This folder contains all the outputs of the analyses in the form of data (`Outputs/Data/`) or figures (`Outputs/Data/`) or Tables (`Outputs/Tables/`).

3. `R/`: The R project consists of codes with individual scripts and functions. All scripts are stored in the `R/` folder and specific functions are stored in a subfolder `R/Functions/`.

3.1 `R/___Init_project___.R`: This script is useful for the initial project setup.
3.2 `R/00_Config.file.R`: This script (referred to as "*Config file*" from here on) is useful for setting up the preferences that are applied to the whole project. Here all settings (configurations) and criteria used throughout the project are predefined by the user before running all the R scripts in the project. In addition, it prepares the current session by loading the required packages and saving all settings throughout the project. Points in the Config file that require the user's attention are flagged by "**[USER]**", meaning that these are criteria that need to be checked by the user.
3.3 `R/Master.R`: If all the input data are available, this script runs the scripts of 01_Data_processing, 02_Main_analyses, and  03_Supplementary_analyses_analyses in a cascade manner. 
3.4 Subfolder `R/01_Data_processing` contains the R codes used for estimating the PAPs, preparing the metadata, and for making the estimated PAP data ready for the analyses. Please note that R scripts within this subfolder will not work because these require the raw data, which is not shared publicly due to the data contributors' embargo.
3.5 Subfolder `R/02_Main_analyses` contains R scripts used for the analyses of the spatio-temporal patterns of the PAPs.
3.6 Subfolder `R/03_Supplementary_analyses` contains R scripts used for the analyses that are used as supplementary material in the paper.

4. `renv/`: This folder stores all the installed packages with the record of their versions.

## How to use the repo?
Once a user obtains the R project, there are several steps to be done before using it:

- Update [R](https://en.wikipedia.org/wiki/R_(programming_language)) and [R-studio IDE](https://posit.co/products/open-source/rstudio/). There are many guides on how to do so (e.g. [here](https://jennhuck.github.io/workshops/install_update_R.html))

- Execute all individual steps with the `___Init_project___.R` script. This will result in the preparation of all R-packages using the [`{renv}` package](https://rstudio.github.io/renv/articles/renv.html), which is an R dependency management of your projects. Mainly it will install a crucial R-package [`{REcopol}`](https://github.com/HOPE-UIB-BIO/R-Ecopol-package) and and all its dependencies. The latest release of {REcopol} is automatically installed in the project set-up stage. Note that installing all packages can take a substantial amount of time.

- Run the `00_Config_file.R` at the beginning of running the R scripts so that the project configuration, required packages, and functions are loaded properly.

