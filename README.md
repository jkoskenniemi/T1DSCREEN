# Genetic evidence for efficacy of targeting IL-2, IL-6, and TYK2 signaling in prevention of type 1 diabetes: A Mendelian randomization study  

Repository for analyses related to article Heikkilä et al. 2024 published in 
Diabetologia (add DOI, link to publisher website, and full bibliographical details once they are available).

1. To get data required for analyses, please run shell script code `/code/00-get-data.sh` by running (in unix systems)
   ```
   chmod u+x ./code/00-get-data.sh
   ./code/00-get-data.sh
   ````
   Note that this downloads data that are publicly available, but does not download data from DECODE, which can be accessed at https://www.decode.com/summarydata/
2. Run `Build-script.R` in R to replicate the code required for analyses
3. The final analyses can be found as Rmarkdown-files in folder `analyses`

  

Created with [workflowr](https://github.com/jdblischak/workflowr)
