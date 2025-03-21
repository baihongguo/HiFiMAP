#!/bin/bash
seed=12345
n_chunks=5
outfile_prefix="toy"
threads=4
min_jobs=5
script_name="Step2_Run_HiFiMAP.R"
chromosome=20 
#chromosome=seq(1 22) 
log_file="FiMAP_log.txt"
infile_prefix="toy"

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
        /usr/bin/time -v Rscript "$script_name" "$chr" "$seed" "$outfile_prefix" "$n_chunks" "$i" "$threads" "$infile_prefix" 2> "$tmp_log" &
       
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
