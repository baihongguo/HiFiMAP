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
The following script transform the IBD segment file into the high performance hdf5 data format to improve the computational efficiency for running HiFiMAP. This step will generate a index file to record the local IBD changing positions (SNP locations), and the index of IBD segments overlapping each changing position(SNP).

```diff 
Rscript Step1_Preprocess_HiFiMAP.R toy_pheno.txt toy_global_kinship.RData toy_chr20_IBD.txt.gz toy 12345 20 1 2 3 4 5 0
```
where the inputs are:
| Value  | Description |
| ------------- | ------------- |
| toy_pheno.txt | The pheno file |
| toy_global_kinship.RData  | The global IBD(kinship) matrix  |
| toy_chr20_IBD.txt.gz  | The IBD segment file  |
| toy  | The prefix of the outout files  |
| 12345  | Random seed  |
| 20  | Chromosome  |
| 1 | Index of individual 1 column |
| 2 | Index of individual 2 column |
| 3 | Index of Chromosome column |
| 4 | Index of IBD start position column |
| 5 | Index of IBD end position column |
| 0 | 0 means unsorted IBD segment file (by start and end position column), 1 means sorted|

This script will generate a .h5 IBD segment file a changing position index file: toy_ibd_chr20.h5 and toy_ibd_chr20_CPos_index.txt

<br />

### Step2 Run HiFiMAP
The following shell script (Run_HiFiMAP_parallel.sh) will perform parallelized HiFiMAP:
```
#!/bin/bash
seed=12345
n_chunks=5
outfile_prefix="/HGCNT95FS/ADDLIE/work/Test/FiMAP/simulation/toy/toy"
threads=4
min_jobs=5
script_name="Step2_Run_HiFiMAP.R"
chromosome=20 
#chromosome=seq(1 22) 
log_file="/HGCNT95FS/ADDLIE/work/Test/FiMAP/simulation/toy/FiMAP_log.txt"

# Ensure the log file exists
echo "Chr, Wall_Time(s), Total_CPU(s)" > "$log_file"

# Start time tracking
start_time=$(date +%s)
job_pids=()  # Store PIDs of running jobs
resource_logs=()  # Store temporary resource logs

for chr in $chromosome; do
#for chr in $(seq 1 22) 
    echo "Processing chromosome $chr "

    # Launch jobs for this chromosome
    for i in $(seq 1 $n_chunks); do
        # Create a unique log file for each Rscript execution
        tmp_log="tmp_log_chr${chr}_chunk${i}.txt"
        resource_logs+=("$tmp_log")
        # Launch job wrapped with /usr/bin/time -v and redirect stderr to the log file
        /usr/bin/time -v Rscript "$script_name" "$chr" "$seed" "$outfile_prefix" "$n_chunks" "$i" "$threads" 2> "$tmp_log" &
       
        job_pids+=($!)  # Store PID of the last background job
    done

    # Monitor running jobs for this chromosome
    while true; do
        running_jobs=$(pgrep -f "$script_name" | wc -l)
        if [ "$running_jobs" -lt "$min_jobs" ]; then
            echo "Less than $min_jobs jobs running. Proceeding to the next chromosome..."
            break
        fi
        sleep 10  # Check every 10 seconds
    done
done

echo "All chromosomes processed. Waiting for remaining jobs to finish..."

# Wait for all jobs for this pheno to finish before merging
while true; do
    running_jobs=$(pgrep -f "$script_name" | wc -l)
    if [ "$running_jobs" -eq 0 ]; then
        echo "All jobs finished. Proceeding with merging..."
        break
    fi
    sleep 10  # Check every 10 seconds
done

# Merge files for each chromosome

for chr in $chromosome; do
#for chr in $(seq 1 22); do
    chunk_files="${outfile_prefix}_chr${chr}_seed${seed}_chunk*.txt"
    output_file="${outfile_prefix}_chr${chr}_seed${seed}.txt"

    if ls $chunk_files 1> /dev/null 2>&1; then
        echo "Merging files for chromosome $chr"
        head -n 1 "$(ls $chunk_files | head -n 1)" > "$output_file"
        tail -n +2 -q $chunk_files | sort -k2,2n >> "$output_file"
        rm $chunk_files
        echo "Merging complete for chromosome $chr"
    else
        echo "No chunk files found for chromosome $chr, skipping merge."
    fi
done

# End time tracking
end_time=$(date +%s)
wall_time=$((end_time - start_time))

# Extract CPU time from all job logs for this pheno
total_cpu=0

for log in "${resource_logs[@]}"; do
        if [ -f "$log" ]; then
            cpu_time=$(awk '/User time/ {sum += $4} END {print sum}' "$log")

            total_cpu=$(echo "$total_cpu + $cpu_time" | bc)

            rm "$log"  # Clean up log file after extracting info
        fi
done

# Log the results
echo "$chromosome, $wall_time, $total_cpu" >> "$log_file"

echo "All processing and merging is complete"
```
where the inputs are:
| Argument  | Description |
| ------------- | ------------- |
| seed  | Random seed  |
| n_chunks  | Number of chunks to divide the IBD segment into, for the purpose of paralized processing. Recommend to set it to 5-10 and set it to 1 means no paralized processing |
| outfile_prefix | The prefix of the outout file |
| min_jobs | This sets the minimum number of jobs currently running since each chunk finish analyis at different time. If less than min_jobs running, it will proceed to the next choromosme for the efficiency of genome-wide HiFiMAP |
| script_name | This will be Step2_Run_HiFiMAP.R|
| chromosome | chromosome |
| log_file | The path and name of the log file, which records the wall time and CPU time for each chromosome |


