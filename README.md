# HiFiMAP
High resolution fast identity-by-descent mapping (HiFiMAP)
<br />
Current version: 1.0.0

## Quick Installation 

`HiFiMAP` can be downloaded via `git clone https://https://github.com/baihongguo/HiFiMAP`

## Preparations
The software is developed using R and tested in Linux environments. The statistical computing software R (>=4.0.0) and the following R packages are required:
* [fastmatch](https://cran.r-project.org/web/packages/fastmatch/index.html) (>=1.1-6)
* [GMMAT](https://cran.r-project.org/web/packages/GMMAT/index.html) (>=1.3.2)
* [rhdf5](https://github.com/grimbough/rhdf5) (>=2.32.4)

## Usage
### Input Files
* #### Pheno file
The file containing the phenotype, individual ID and covariates:
```diff 
id covar pheno
0 -12.8758 0.137011197818403
1 -12.1141 -0.280733839831654
2 -12.3981 -0.103265188894815
3 -11.5917 0.224730939351052
```

* #### IBD segments file
The output IBD segment file from the IBD detection tools, containing at least the following columns:
```diff 
individual_1_id individual_2_id chromosome_id ibd_physical_position_start ibd_physical_position_end
1004 5 20 66532 354440
975 21 20 66532 354440
975 36 20 66532 354440
663 38 20 66532 354440
689 38 20 66532 354440
```


* #### Global IBD matrix
The global IBD(kinship) matrix saved in . RData format



### Step1 Preprocess the IBD segment file
The following script transform the IBD segment file into the high performance hdf5 data format to improve the computational efficiency for running HiFiMAP. This step will generate a .h5 IBD segment file and a .txt index file to record the local IBD changing positions (SNP locations) and indices of IBD segments overlapping each changing position(SNP).

```diff 
Rscript Step1_Preprocess_HiFiMAP.R toy_pheno.txt toy_global_kinship.RData toy_chr20_IBD.txt.gz toy 12345 1 2 3 4 5 0
```
where the inputs are:
| Value  | Description |
| ------------- | ------------- |
| toy_pheno.txt | The pheno file |
| toy_global_kinship.RData  | The global IBD(kinship) matrix  |
| toy_chr20_IBD.txt.gz  | The IBD segment file  |
| toy  | The prefix of the output files  |
| 12345  | Random seed  |
| 1 | Index of individual 1 column |
| 2 | Index of individual 2 column |
| 3 | Index of Chromosome column |
| 4 | Index of IBD start position column |
| 5 | Index of IBD end position column |
| 0 | 0 means unsorted IBD segment file (by start and end position column), 1 means sorted|

This script will generate a .h5 IBD segment file a changing position index file: toy_ibd_chr20.h5 and toy_ibd_chr20_CPos_index.txt

<br />

### Step2 Run HiFiMAP
The shell script Run_HiFiMAP_parallel.sh will perform parallelized HiFiMAP, and the inputs are:

```
#!/bin/bash
seed=12345
n_chunks=5
outfile_prefix="toy"
threads=4
min_jobs=5
script_name="Step2_Run_HiFiMAP.R"
chromosome=20 
#chromosome=seq(1 22) 
log_file="HiFiMAP_log.txt"
infile_prefix="toy"

```
| Argument  | Description |
| ------------- | ------------- |
| seed  | Random seed  |
| n_chunks  | Number of chunks to divide the IBD segment into, for the purpose of paralized processing. Recommend to set it to 5-10, set it to 1 means no paralized processing |
| outfile_prefix | The prefix of the outout file |
| min_jobs | This sets the minimum number of jobs currently running since each chunk finish analyis at different time. If less than min_jobs running, it will proceed to the next choromosme for the efficiency of genome-wide HiFiMAP |
| script_name | This will be Step2_Run_HiFiMAP.R|
| chromosome | chromosome |
| log_file | The path and name of the log file, which records the wall time and CPU time for each chromosome |
| infile_prefix | This will be the outfile_prefix in step 1 |

Please modify these inputs in Run_FiMAP_parallel.sh for your own analysis. This script will generate a log file and a result file containing the following columns:
```diff 
chr pos n.ibd.segs p.value
20 66532 8022 0.215780531425396
20 71723 8286 0.235451653255488
20 76361 8425 0.273023093580426
20 82720 8543 0.28431254447582
20 87903 8557 0.289574216210267
```

### Example
Here is an example running HiFiMAP using the toy data set in the example folder
* #### Step1 Preprocess the IBD segment file
```diff 
Rscript Step1_Preprocess_HiFiMAP.R toy_pheno.txt toy_global_kinship.RData toy_chr20_IBD.txt.gz toy 12345 1 2 3 4 5 0
```

* #### Step2 Run HiFiMAP
```diff 
bash Run_FiMAP_parallel.sh
```
