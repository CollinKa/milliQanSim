#!/bin/bash

# Base directories
input_dir="/homes/dgimani/milliQanSim/inputData/Run3_G4_mcp_inputs/c_0p1"
#masses=("m_0p1" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" "m_1p7" )
masses=("m_0p7" "m_2p5")
#masses=("m_2p0" "m_2p3" "m_2p5" "m_2p7" "m_3p0" "m_3p5" "m_3p7" "m_4p0" "m_4p5" "m_4p7" "m_5p0" "m_5p5" "m_6p0" "m_6p5" "m_7p0" "m_8p0" "m_9p0" "m_10p0" "m_12p0" "m_15p0" "m_17p0" "m_20p0" "m_25p0" "m_27    p0" "m_33p0" "m_35p0" "m_40p0")

# Loop over each subdirectory 
for mass_dir in "$input_dir"/*; do
  echo "$mass_dir"
  #if [[ " ${masses[@]} " =~ " ${mass_dir} " ]]; then
  if [ -d "$mass_dir" ]; then
    # Extract the mass (subdirectory name)
    mass=$(basename "$mass_dir")

    echo "$mass"
    output_file="$input_dir/$mass/mcp_${mass}_G4_input.txt"
    > "$output_file"

    # Loop over each mode subdirectory inside the mass directory
    for mode_dir in "$mass_dir"/*; do
      if [ -d "$mode_dir" ]; then
        # Extract the mode (subdirectory name)
        mode=$(basename "$mode_dir")
	echo "$mode"
        # Ensure the corresponding directory exists in the output directory
        #mkdir -p "$output_dir/$mass/$mode"

	mode_file="$input_dir/$mass/$mode/output_mcp_${mass}_${mode}.txt"
	echo "$mode_file"

	# Check if the output file exists and concatenate it
        if [ -f "$mode_file" ]; then
          echo "Appending $mode_file to $output_file"
          cat "$mode_file" >> "$output_file"
        else
          echo "Filtered file not found: $mode_file"
        fi


      fi
    done
  fi
  #fi
done
