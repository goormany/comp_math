#!/bin/bash

GEN=./gen_matrix
SOLVER=./jacobi

RUNS=20
MIN_SIZE=3
MAX_SIZE=1000

types=("dd" "not_dd" "no_solution")

echo "Random size stress test"
echo "Runs per type: $RUNS"
echo "Size range: [$MIN_SIZE, $MAX_SIZE]"
echo "----------------------------------------"

for type in "${types[@]}"; do
    echo "Testing type: $type"

    total_time=0
    converged=0
    failed=0

    for ((i=1; i<=RUNS; i++)); do
        size=$(( RANDOM % (MAX_SIZE - MIN_SIZE + 1) + MIN_SIZE ))
        outfile="tmp_${type}.txt"

        echo "  Run $i | size=$size"

        # Generate matrix
        $GEN --outpath "$outfile" --size $size --type $type

        # Time using high precision
        start=$(date +%s.%N)

        output=$($SOLVER "$outfile" 2>&1)
        status=$?

        end=$(date +%s.%N)

        # Compute elapsed time (floating point)
        elapsed=$(echo "$end - $start" | bc)
        total_time=$(echo "$total_time + $elapsed" | bc)

        # Detect failure (crash / NaN / instability)
        if [[ $status -ne 0 ]] || echo "$output" | grep -qi "nan"; then
            ((failed++))
            continue
        fi

        # Detect convergence
        if echo "$output" | grep -q "Converged"; then
            ((converged++))
        fi
    done

    avg_time=$(echo "scale=4; $total_time / $RUNS" | bc)

    echo "----------------------------------------"
    echo "Results for $type:"
    echo "  Total time: $total_time sec"
    echo "  Avg time:   $avg_time sec"
    echo "  Converged:  $converged / $RUNS"
    echo "  Failed:     $failed / $RUNS"
    echo "----------------------------------------"
done

rm -f tmp_*.txt
