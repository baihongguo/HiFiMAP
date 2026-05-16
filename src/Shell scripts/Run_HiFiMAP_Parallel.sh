#!/bin/bash

# ==========================================
# CONFIGURATION
# ==========================================
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
