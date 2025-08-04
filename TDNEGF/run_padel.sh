set -uo pipefail

# Directory for timing outputs
OUTPUT_DIR="./timings"
mkdir -p "$OUTPUT_DIR"

# Logarithmically spaced t_end values between 10 and 1000
t_end_values=(10 31 100 316 1000 3162 10000)

# Path to your Julia script (adjust if located elsewhere)
script_path="./Padel.jl"

for t_end in "${t_end_values[@]}"; do
    echo "Running simulation with n=10, t_end=$t_end..."
    output_file="$OUTPUT_DIR/padel_time_${t_end}.txt"

    # Run Julia and capture both stdout and stderr (timing info)
    # Use || true to prevent early exit if Julia returns non-zero
    if ! julia "$script_path" 10 "$t_end" > "$output_file" 2>&1; then
        echo "⚠ Error encountered for t_end=$t_end; see $output_file for details"
    else
        echo "✓ Completed t_end=$t_end, results saved to $output_file"
    fi

done

echo "All runs attempted. Timing files are in the '$OUTPUT_DIR' folder."