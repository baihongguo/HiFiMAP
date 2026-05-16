# HiFiMAP
High resolution fast identity-by-descent (IBD) mapping test
<br />
Current version: 2.0.0 (Optimized Streaming Pipeline)

## Overview
The new HiFiMAP pipeline has been rebuilt for scalability and efficiency. It uses a highly optimized Python preparser to convert IBD segments into chunked sparse matrices (`.mtx`) and differential updates (`.diff`). The main R package utilizes C++ (`Rcpp`) to stream these updates, allowing for fast, parallelized genome-wide association testing without frequently loading massive IBD matrices into RAM.

## Quick Installation 

`HiFiMAP` can be downloaded via:
```bash
git clone [https://github.com/hanchenlab/HiFiMAP](https://github.com/hanchenlab/HiFiMAP)
```

## Preparations
The software is developed and tested in Linux HPC environments. 

**Python Dependencies (Python >= 3.6):**
* `numpy`

**R Dependencies (R >= 4.0.0):**
* `data.table`
* `Matrix`
* `fastmatch`
* `GMMAT`
* `Rcpp`
* `CompQuadForm`

---

## Usage & Example Pipeline

Based on the repository structure, all test files are located in the `example/` directory, and all source codes are in the `src/` directory. Suppose you want to run a test analysis for chromosome 21 using the provided toy data:

* `example/chr21_toy.ibd.gz`: IBD segments generated from [hapIBD](https://github.com/browning-lab/hap-ibd).
* `example/chr21_toy.vcf.gz`: The corresponding VCF file used to generate the IBDs.
* `example/toy_phenotype.txt`: Contains individual IDs, phenotype, and covariates (e.g., age, sex).
* `example/toy_global_kinship.RData`: The global IBD (kinship) matrix.

### Step 0: Pre-parse the IBD Segments (Python)
This step transforms the raw IBD segments into high-performance sparse matrices and delta differences. It automatically chunks the chromosome to allow for parallel processing in Step 2.

Run the parser from the root of the repository:

```bash
python3 src/hapIBD_parsing/parsing_hapIBD.py \
    --ibd example/chr21_toy.ibd.gz \
    --vcf example/chr21_toy.vcf.gz \
    --output example/ibd_prep/chr21 \
    --n-checkpoints 20
```

| Argument | Description |
| ------------- | ------------- |
| `--ibd` | The IBD segment file output from hapIBD. |
| `--vcf` | The reference VCF file to extract exact SNP positions and sample IDs. |
| `--output` | The directory where the chunked matrices and diffs will be saved. |
| `--n-checkpoints` | Number of chunks to divide the chromosome into (e.g., `20`). This controls parallelization scaling in Step 2. Set this number to be less than the number of CPUs you have for better performance, recommended is 20-30. |

*(This will generate a folder `example/ibd_prep/chr21` containing `sites.txt` (bp of each testing site), `samples.txt` (individual IDs for matching with the null model), and the chunked `.mtx` / `.diff` files).*

<br />

### Step 1: Fit the Null GLMM Model (R)
Next, fit the null model using your phenotypes and global kinship matrix. **This step only needs to be run once per phenotype**, regardless of how many chromosomes you analyze. Please modify the glmmkin formula in Step1_Save_null_model_HiFiMAP.R based on your phenotype and covariate names. ID column should be named id.

```bash
Rscript src/HiFiMAP/Step1_Save_null_model_HiFiMAP.R \
    example/toy_phenotype.txt \
    example/toy_global_kinship.RData \
    toy_pheno \
    12345
```

| Argument | Description |
| ------------- | ------------- |
| `phenotype.txt` | The file containing ID, Phenotype, and covariates. |
| `global_kinship.RData` | The kinship matrix saved in `.RData` format. Row names and col names should be the corresponding subject IDs |
| `toy_pheno` | Prefix for the output file. |
| `12345` | Random seed for reproducibility. |

*(This generates the null model object: `toy_pheno_glmmkin2randomvec.rds` in your root directory. The subjects in the phenotype file can be a subset of the individuals included in your vcf files used for step0 and hapIBD).*

<br />

### Step 2: Run Parallelized HiFiMAP Scan (Bash wrapper)
To run the actual association scan efficiently across all chunks, you can use the `Run_HiFiMAP_parallel.sh` wrapper script. 

Make sure the configuration inside `Run_HiFiMAP_parallel.sh` points to the correct `src` paths and `example` outputs:

```bash
#!/bin/bash

# ==========================================
# CONFIGURATION
# ==========================================
seed=12345
# Set the number of threads for matrix operations in R via Sys.setenv(MKL_NUM_THREADS = threads).
# NOTE: This is only effective if R is compiled with the Intel Math Kernel Library (MKL).
threads=2
target_chromosomes="21"

OUT_BASE="./results"
IBD_PREP_BASE="./example/ibd_prep"
GLMM_RDS="toy_pheno_glmmkin2randomvec.rds" 
script_name="src/HiFiMAP/Step2_Run_HiFiMAP.R"

mkdir -p "$OUT_BASE"

# ==========================================
# EXECUTION LOOP
# ==========================================
for chr in $target_chromosomes; do
    ibd_out_dir="${IBD_PREP_BASE}/chr${chr}"
    
    echo "=================================================="
    echo "Processing Chromosome $chr..."
    echo "=================================================="

    if [ ! -d "$ibd_out_dir" ]; then
        echo "  [ERROR] Input dir not found for Chr $chr: $ibd_out_dir"
        continue
    fi

    # Read the ACTUAL site numbers from the freshly generated Python files
    all_checkpoints=($(ls "${ibd_out_dir}"/ibd_mat_*.mtx 2>/dev/null | xargs -n 1 basename | sed 's/ibd_mat_//;s/\.mtx//' | sort -n))
    num_avail_checkpoints=${#all_checkpoints[@]}

    if [ "$num_avail_checkpoints" -eq 0 ]; then
        echo "  [ERROR] No .mtx checkpoints found in $ibd_out_dir. Skipping."
        continue
    fi

    total_sites=$(tail -n +2 "${ibd_out_dir}/sites.txt" | wc -l)
    
    echo "  [LAUNCH] Chr $chr | Spawning exactly $num_avail_checkpoints parallel chunks..."

    # --- LAUNCH CHUNKS ---
    for (( j=0; j<num_avail_checkpoints; j++ )); do
        
        start_idx=${all_checkpoints[$j]}
        
        # Calculate End Index
        next_job_index=$((j + 1))
        if [ "$next_job_index" -lt "$num_avail_checkpoints" ]; then
            end_idx=$(( ${all_checkpoints[$next_job_index]} - 1 ))
        else
            end_idx=$((total_sites - 1))
        fi

        chunk=$((j+1))
        outfile_prefix="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${chunk}"

        Rscript "$script_name" \
            "$chr" \
            "$seed" \
            "$outfile_prefix" \
            "$threads" \
            "$GLMM_RDS" \
            "$ibd_out_dir" \
            "$start_idx" \
            "$end_idx" > "${outfile_prefix}.log" 2>&1 &
            
    done

    echo "  [MONITOR] All chunks dispatched. Waiting for completion..."
    wait 

    echo "  [MERGE] Merging results for Chromosome $chr..."
    final_out="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}.txt"
    chunk1="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk1.txt"
    
    if [ -f "$chunk1" ]; then
        head -n 1 "$chunk1" > "$final_out"
        count=0
        for (( j=1; j<=num_avail_checkpoints; j++ )); do
            f="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${j}.txt"
            if [ -f "$f" ]; then
                tail -n +2 -q "$f" >> "$final_out"
                rm "$f" 
                ((count++))
            fi
        done
        echo "  [SUCCESS] Chr $chr: Merged $count chunks -> $final_out"
    fi
done

echo "Pipeline Complete!"
```

### Final Output
The final merged result file (`results/HiFiMAP_toy_pheno_chr21.txt`) will contain the association statistics for every tested site:

```text 
chr     pos       n.ibd.segs    p.value
21	14240779	1052	0.8847846	
21	14242245	1054	0.8866947	
21	14250602	1056	0.8867795	
21	14256336	1057	0.8965234
```
