# LDER-GE
We propose a statistical method to estimate the phenotypic variance explained by Gene-Envieonment interactions. The method is called  Linkage-Disequilibrium Eigenvalue Regression for Gene-Environment interactions (LDER-GE). LDER-GE extends existing LDSC-based methods to incorporate full LD information and enhances statistical efficiency of estimation.

:open_book: Citation:

Zihan Dong, Wei Jiang, Hongyu Li, Andrew T DeWan, Hongyu Zhao, LDER-GE estimates phenotypic variance component of gene–environment interactions in human complex traits accurately with GE interaction summary statistics and full LD information, Briefings in Bioinformatics, Volume 25, Issue 4, July 2024, bbae335, https://doi.org/10.1093/bib/bbae335

Acknowledgement: LDER-GE package is forked based on the original LDER package. If you are studying narrow-sense heritability, please use and refer to https://github.com/shuangsong0110/LDER). We modified the function calling procedure and adjusted the algorithm specifically for GE interaction analysis. For LD preparation, we largely maintain the original LDER framework.
We thank Shuang Song for sharing the original LDER code.

## Table of contents
* [Install and LD preparation](#hammer-install-and-ld-preparation)
* [Estimation of GE proportion](#rocket-estimation-of-ge-proportion)
* [Output](#bulb-output)
* [A Simplified Pipeline](#key-a-simplified-pipeline)

## :hammer: Install and LD preparation
LDER-GE R package requires R >= 3.5.0 and Python 3.6.
LDER-GE can be installed using the command:
```r
devtools::install_github('dongzhblake/LDER-GE')
```

We provide a function `plinkLD.py` for efficient LD information extraction and shrinkage based on Python. Users could either specify their own LD reference files with plink bfile format (.bim, .fam, .bed), or use the pre-computed LD information. When studying UKBB or any Biobank-scale dataset, we suggest use in-sample LD panel for better estimatuon results. Two examples are given below.

### Example 1: Use the pre-computed LD information

The pre-computed LD information for 396,330 hapmap3 variants from 276,050 UK Biobank European individuals can be manually downloaded from [https://drive.google.com/file/d/1mvDA79qPAoPXUjmUC1BQw-tInklZ4gPD/view?usp=drive_link](https://drive.google.com/file/d/1CCGil-ZnXourrk5JFJeqyEKIzLCvX2vh/view?usp=drive_link)

The pre-computed LD information for 966,766 hapmap3 and array variants from 307,259 UK Biobank European individuals can be manually downloaded from [https://drive.google.com/file/d/1UF1xP1Rg1JiFMozkFJ3bbJmFggDY8-5Y/view?usp=drive_link](https://drive.google.com/file/d/1FgikyxYd_jW05aLuHgj0yChW9i2ZgtT9/view?usp=sharing)

After downloading, decompress the files:

`unzip UKB396kvariant_hm3.zip`

`unzip UKB966kvariant_hm3.zip`

:exclamation: NOTE: Please keep the download path the SAME with that used in function `runLDER_GE`.


### Example 2: Use plink bfile as the input
The 1000 Genome Project reference panel (hg19) can be downloaded by:

`wget https://zenodo.org/record/7768714/files/1000G_Phase3_plinkfiles.tgz?download=1`

`tar -xvzf 1000G_Phase3_plinkfiles.tgz`


```r
library(LDERGE)
generateLD(assoc=GWAS_SUMMARY_STATISTICS (required), 
          path=OUTPUT_DIR (required),
          bfile_path=PATH_TO_LD_REFERENCE (required),
          cores=NUMBER_OF_CORES (optional),
          ethnic=ETHNIC (optional),
          plink_path=PATH_TO_PLINK_SOFTWARE (optional),
          python_path=PATH_TO_PYTHON_SOFTWARE (optional))                    
```
- GWAS_SUMMARY_STATISTICS (required): GWAS summary statistics, need to include `snp`, `chr`, `a0`, `a1`, `z` (header names are necessary)

- OUTPUT_DIR (required): The output path

- PATH_TO_LD_REFERENCE (required): The LD reference plink bfile. If the files are divided into chromosomes, the function should be run for each of the file.

- NUMBER_OF_CORES (optional): The number of cores for computation in parallel.

- ETHNIC (optional): Ethnic of the GWAS cohort; 'eur' for European ancestry.

- PATH_TO_PLINK_SOFTWARE (optional): The path to the plink software. If not specified, the function will use the default path (system("which plink"))

- PATH_TO_PYTHON_SOFTWARE (optional): The path to the python software. If not specified, the function will use the default path (system("which python"))

Note: The function will automatically download and install plinkLD functions. If it does not work, the packages can also be downloaded with `wget -O plinkLD.zip https://cloud.tsinghua.edu.cn/f/7c002e9b9539450182ef/?dl=1 --no-check-certificate` 


## :rocket: Estimation of GE proportion
The main funcion can be run with:

```r
runLDER_GE(assoc=GWAS_SUMMARY_STATISTICS (required), 
	n.gwas=SAMPLE_SIZE_OF_GWAS (required), 
	path=OUTPUT_DIR (required),
	LD.insample=IN_SAMPLE_LD (T/F, required),
	n.ld=SAMPLE_SIZE_OF_LD_REF (required), 
	method=METHOD (default='lder')
	cores=NUMBER_OF_CORES (optional))
```
- GWAS_SUMMARY_STATISTICS (required): GWAS summary statistics, need to include `snp`, `chr`, `a0`, `a1`, `z` (header is necessary)

- n.gwas (required): The sample size of the GWAS summary statistics

- OUTPUT_DIR (required): The output path (Note that the path should be SAME with that used in function `generateLD`)

- IN_SAMPLE_LD (required): T/F, whether the LD reference is estimated with the GWAS cohort (T) or external reference panel (e.g. 1000 Genome Project: F)

- SAMPLE_SIZE_OF_LD_REF (required): The sample size of the LD reference (e.g., 276,050 for 396K UKBB variants; 307,259 for 966K UKBB variants, 489 for 1000G)

- METHOD (optional): Default='lder'. We also provide a choice of 'both', which outputs the results for both LDER and LDSC.

- NUMBER_OF_CORES (optional): The number of cores for computation in parallel.




## :bulb: Output

If `method='lder'`, the `runLDER_GE` function returns a list with 5 elements:

`h2I`: Estimated GE proportion by LDER-GE

`h2I.se`: The standard error of estimated GE proportion with block-jackknife.

`h2I.p`: The P value for testing the estimated GE proportion.

`intecept`: Estimated intercept by LDER-GE

`intecept.se`: The standard error of estimated intercept with block-jackknife.

If `method='both'`, the `runLDER` function returns a list containing the results of both LDER-GE and LDSC-based methods.


## :key: A Simplified Pipeline
Download a sample GWIS summary statistics at [https://drive.google.com/file/d/1W1zaoOS3ob0dzSvJL9p2f4kdKUmgW5Ff/view?usp=drive_link](https://drive.google.com/file/d/1LSkzmfZo36WtiIpXBu0Cfw8iA4VwDss7/view?usp=sharing)


Download the pre-computed LD information for 396,330 hapmap3 variants from 276,050 UK Biobank European individuals from [https://drive.google.com/file/d/1mvDA79qPAoPXUjmUC1BQw-tInklZ4gPD/view?usp=drive_link](https://drive.google.com/file/d/1mvDA79qPAoPXUjmUC1BQw-tInklZ4gPD/view?usp=drive_link)


`unzip UKB396kvariant_hm3.zip`


Run with R:

```r
devtools::install_github('dongzhblake/LDER-GE')
library(LDERGE)
library(data.table)
path0 <- "UKB396kvariant_hm3" # the complete system path to this LD folder
assoc <- fread('LDER_GE_exampleGWIS.txt')
# The whole process of runLDER_GE is going to take a few minutes depending on the number of cores of the computer.
# If a higher number of cores are available, the parallel input of summary statistis will be faster.
res <- runLDER_GE(assoc, n.gwis=median(assoc$n), path=path0, n.ld=276050, cores=10, method='lder')

> unlist(res)
        lder.h2I    lder.intecept      lder.h2I.se       lder.h2I.p 
    5.193470e-02     1.143381e+00     7.510239e-03     4.672499e-12 
lder.intecept.se         ldsc.h2I    ldsc.intecept      ldsc.h2I.se 
    5.513937e-03     5.558019e-02     1.137570e+00     9.309304e-03 
      ldsc.h2I.p ldsc.intecept.se 
    2.366851e-09     9.410519e-03 

```


## :busts_in_silhouette: Maintainer

Please contact Zihan Dong (zihan.dong@yale.edu) if there are any problems or questions.


