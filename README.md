# LDER-GE
We propose a statistical method to estimate the phenotypic variance explained by Gene-Envieonment interactions. The method is called  Linkage-Disequilibrium Eigenvalue Regression for Gene-Environment interactions (LDER-GE). LDER-GE extends existing LDSC-based methods to incorporate full LD information and enhances statistical efficiency of estimation.

:open_book: Citation:

Dong, Z., Jiang, W., Li, H., Dewan, A. T., & Zhao, H. (2023). LDER-GE estimates phenotypic variance component of gene-environment interactions in human complex traits accurately with GE interaction summary statistics and full LD information. bioRxiv, 2023-11.

Acknowledgement: This LDER-GE package is modified based on the original LDER package. If you are doing narrow-sense heritability study, please use and refer to https://github.com/shuangsong0110/LDER). We modified function calling procedure and adjusted the algorithm specifically for GE interaction analysis. For LD preparation, we largely maintain the original LDER framework.
We thank Shuang Song for sharing the original LDER code.

## Table of contents
* [Install](#hammer-install)
* [LD prepared](#scroll-ld-prepared)
* [Estimation of GE proportion](#rocket-estimation-of-ge-proportion)
* [Output](#bulb-output)
* [A Simplified Pipeline](#key-a-simplified-pipeline)

## :hammer: Install
R >= 3.0.0

Python 3.6

LDER-GE is an R package which can be installed using the command:
```r
devtools::install_github('dongzhblake/LDER-GE')
```

## :scroll: LD prepared
We provide a function `plinkLD.py` for efficient LD information extraction and shrinkage based on Python. 
Users could either specify their own LD reference files with plink bfile format (.bim, .fam, .bed), or use the pre-computed LD information. We provide two examples here.


:exclamation: NOTE: We suggest users use plink bfile as the input, because the different numbers of SNPs in GWAS and in the reference panel may lead to a slight difference in the LD shrinkage.

### Example 1: Use plink bfile as the input (recommended)
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

### Example 2: Use the pre-computed LD information

The pre-computed LD information of 276,050 UK Biobank European individuals can be downloaded by

`wget -O LD.shrink.zip https://cloud.tsinghua.edu.cn/f/abf1020acb9c435eaa13/?dl=1 --no-check-certificate`

`wget -O LD.zip https://cloud.tsinghua.edu.cn/f/d93a4a7013fe461aa9fc/?dl=1 --no-check-certificate`

`unzip LD.shrink.zip`

`unzip LD.zip`

:exclamation: NOTE: Please keep the download path the SAME with that used in function `runLDER`.


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

- SAMPLE_SIZE_OF_LD_REF (required): The sample size of the LD reference (e.g., 489 for 1000G)

- METHOD (optional): Default='lder'. We also provide a choice of 'both', which outputs the results for both LDER and LDSC.

- NUMBER_OF_CORES (optional): The number of cores for computation in parallel.




## :bulb: Output

If `method='lder'`, the `runLDER_GE` function returns a list with 4 elements:

`h2`: Estimated GE proportion by LDER-GE

`inf`: Estimated intercept by LDER-GE

`h2.se`: The standard error of estimated GE proportion with block-jackknife.

`inf.se`: The standard error of estimated intercept with block-jackknife.

If `method='both'`, the `runLDER` function returns a list containing the results of both LDER-GE and LDSC-based methods.


## :key: A Simplified Pipeline
Download a sample GWAS summary statistics:

$ wget -O gwas_sample.txt https://cloud.tsinghua.edu.cn/f/828ab71c87d84dd28d47/?dl=1 --no-check-certificate

Download 1000G LD reference:

$ wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz

$ tar -xvzf 1000G_Phase3_plinkfiles.tgz


Run with R:

```r
devtools::install_github('dongzhblake/LDER-GE')
library(LDER)
library(data.table)
path0 <- getwd()
assoc <- fread('gwas_sample.txt')
for(chr in 1:22){
    generateLD(assoc, path = path0, bfile_path = paste0('./1000G_EUR_Phase3_plink/1000G.EUR.QC.', chr))
}
res <- runLDER_GE(assoc, n.gwas=2e4, path=path0, LD.insample=F, n.ld=489, cores=10, method='lder')

```


## :busts_in_silhouette: Maintainer

Please contact Zihan Dong (zihan.dong@yale.edu) if there are any problems or questions.


