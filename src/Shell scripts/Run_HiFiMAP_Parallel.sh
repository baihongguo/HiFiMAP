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
