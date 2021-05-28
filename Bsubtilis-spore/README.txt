Pipeline to analyze flowcytometry data of SYBR green stained culture of Bacillus subtilis.
Quatification relies on sample volume measurment by NovCyte machine.

use R version 3.6.3
https://support.rstudio.com/hc/en-us/articles/200486138-Changing-R-versions-for-RStudio-desktop

After switching the R version, open flow-cyto-methods/Bsubtilis-spore/Bsubtilis-spore.Rproj in Rsudio

This pipeline has many dependencies.
I have used the "renv" package to track all the dependencies.
learn more at: https://rstudio.github.io/renv/articles/renv.html

First time you wil need to insatll the packages by running the following commands (within project)
library(renv)
renv::restore

I also use the "here" package to call files using relative paths.
learn more at: https://here.r-lib.org/

To start analysis:
1. Export flow data to fcs file.
I assume there replicate data that are placed in folders with there names being the same except for well and rep IDs.
2. Copy the fcs files to: flow-cyto-methods/Bsubtilis-spore/data/FCM/fcs
3. To analyze specific groups of fcs files place them in a subfolder within fcs/ folder, and change the variavlble "day" assignment to name of folder.
4. Change the variavlble "dilution" assignment to culture dilution you made i forma "x100" for a 100-fold dilution.
5. change the variavlble "sample.var" assignment to list the samle info to be parsed form file name.
6. run the FOR loop

Results:
for each group of replicates the pipeline should make 4 plots in flow-cyto-methods/Bsubtilis-spore/figs/gate_plots folder:
> singlet: FSC height vs. area. used to exclude doublets
> noise: dustribution if FSC are values, with line marking the noise filtring threshold.
> scatterNoise: events scatter plot on SYBR fluoresence vs SSC area with events filtered out by noise marked in grey, and events passed to clustering in blue.
> cluster: Events scatter plot on SYBR fluoresence vs FSC area with assignment to spore and veg clusters given by color. Blue dots represnt the model based cluster center for veg and spore.

The actual counts. along with pipeline values collected along the way, are exported to flow-cyto-methods/Bsubtilis-spore/data/output as .csv file.

*Needs work still*
analysis_post_FCM.R can be used to bnd results from csv file to singlr data fram for further analysis.

 


