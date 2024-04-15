# Supporting Information

This repository contains the source code to reproduce the calculations and plots of the following manuscript: *Maloney AE, Kopf SH, Zhang Z, McFarlin J, Nelson DB, Masterson AL, Zhang X. Large enrichments in fatty acid 2H/1H ratios distinguish respiration from aerobic fermentation in yeast Saccharomyces cerevisiae. Proceedings of the National Academy of Sciences.*

DOI: https://doi.org/10.1073/pnas.2310771121

# What can I do with this code?

In publishing this repository, our hope is that this code is useful to other members of the scientific community. This repository is released under a [Creative Commons BY (CC-BY) license](https://creativecommons.org/licenses/by/4.0/), which means that all code published here can be shared and adapted for any purposes so long as appropriate credit and citation of the original paper is given. See attribution section for details.

# How do I run this code?

1. Download and install [R for your operating system](https://cloud.r-project.org/).
2. Download and install [RStudio for your operating system](https://posit.co/download/rstudio-desktop/).
3. Download a [zip file of this repository](https://github.com/KopfLab/2024_maloney_et_al/archive/refs/heads/main.zip) and decompress it in a directory of your choosing on your computer.
4. Navigate to the `stats_model_SK` folder and open the `project.Rproj` file to start Rstudio and load this project's stats and model files.
5. Open the script `yeast_data_analysis.qmd`. 
6. Ensure that you have all of the required libraries installed by running `install.packages(c("tidyverse", "readxl", "cowplot"))` from the R console. If any libraries fail to install, note the name of the library and attempt to manually install its most recent version via CRAN or GitHub.
7. To generate an HTML report, select `File` --> `Render Document` from the menu. 
8. Alternatively, navigate to the `main_figs_AM/SCRIPTS` folder and run either `d2H d13C data (Fig 2).R` to reconstruct main text figure 2 or `E data (Fig3).R` to reconstruct main text figure 3. 



