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
IBD_PREP_BASE="./ibd_prep"
GLMM_RDS="toy_pheno_glmmkin2randomvec.rds"
script_name="Step2_Run_HiFiMAP_Direct_Parallel_flexible.R"

mkdir -p "$OUT_BASE"

# ==========================================
# EXECUTION
# ==========================================
for chr in $target_chromosomes; do
    echo "Launching chunks for Chromosome $chr..."
    
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
