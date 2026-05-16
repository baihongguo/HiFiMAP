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
| `global_kinship.RData` | The global kinship matrix saved in `.RData` format. Row names and col names should be the corresponding subject IDs |
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
n_analysis_chunks=20   # MUST MATCH the --n-checkpoints used in Step 0
min_jobs=4             # Throttle: launch next chunks when active jobs < min_jobs
seed=12345
threads=4
target_chromosomes="21"

OUT_BASE="./results"
IBD_PREP_BASE="./example/ibd_prep"
GLMM_RDS="toy_pheno_glmmkin2randomvec.rds"
script_name="src/HiFiMAP/Step2_Run_HiFiMAP.R"

mkdir -p "$OUT_BASE"

# ==========================================
# EXECUTION
# ==========================================
for chr in $target_chromosomes; do
    echo "=========================================="
    echo "Processing Chromosome $chr..."
    echo "=========================================="
    
    # ---------------------------------------------------------
    # ALIGNMENT FIX: Rename Python outputs to sequential chunks
    # ---------------------------------------------------------
    echo "Aligning .mtx and .diff files to sequential chunks..."
    chunk_count=1
    
    # Use sort -V to ensure natural numeric sorting (e.g., 250 comes before 1000)
    for mtx_file in $(ls ${IBD_PREP_BASE}/chr${chr}/ibd_mat_*.mtx 2>/dev/null | sort -V); do
        # Extract the actual site number from the filename
        site_num=$(basename "$mtx_file" | grep -o -E '[0-9]+')
        
        # Only rename if it isn't already the correct sequential chunk number
        if [ "$site_num" != "$chunk_count" ]; then
            mv "${IBD_PREP_BASE}/chr${chr}/ibd_mat_${site_num}.mtx" "${IBD_PREP_BASE}/chr${chr}/ibd_mat_${chunk_count}.mtx"
            
            # Also rename the corresponding .diff file
            if [ -f "${IBD_PREP_BASE}/chr${chr}/delta_${site_num}.diff" ]; then
                mv "${IBD_PREP_BASE}/chr${chr}/delta_${site_num}.diff" "${IBD_PREP_BASE}/chr${chr}/delta_${chunk_count}.diff"
            fi
        fi
        ((chunk_count++))
    done
    echo "Alignment complete."
    # ---------------------------------------------------------

    echo "Launching jobs for Chromosome $chr..."
    
    for chunk in $(seq 1 $n_analysis_chunks); do
        
        # Queue Management
        running_jobs=$(pgrep -f "$script_name" | wc -l)
        if [ "$running_jobs" -ge "$min_jobs" ]; then
            while [ "$running_jobs" -ge "$min_jobs" ]; do
                sleep 5
                running_jobs=$(pgrep -f "$script_name" | wc -l)
            done
        fi
        
        # Launch Chunk in background
        out_file="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${chunk}.txt"
        
        Rscript "$script_name" \
            "$chr" "$seed" "toy_pheno" "$threads" \
            "$GLMM_RDS" "${IBD_PREP_BASE}/chr${chr}" \
            "$chunk" "$n_analysis_chunks" > "${out_file}.log" 2>&1 &
            
    done
    
    # Wait for all chunks on this chromosome to finish
    wait 
    
    # Merge chunks
    echo "Merging results for Chromosome $chr..."
    final_out="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}.txt"
    chunk1="${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk1.txt"
    
    if [ -f "$chunk1" ]; then
        head -n 1 "$chunk1" > "$final_out"
        for chunk in $(seq 1 $n_analysis_chunks); do
            tail -n +2 -q "${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${chunk}.txt" >> "$final_out"
            rm "${OUT_BASE}/HiFiMAP_toy_pheno_chr${chr}_chunk${chunk}.txt"* # Cleanup
        done
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
21      9411300   104811        0.00312
```
