#!/bin/bash

# Paths
UNIPROT_DIR="/alphafold_data"
DAT_DIR="/dali_files/dat_files"
IMPORT="$(readlink -f ../../../pdb/DaliLite.v5/bin/import.pl)"
CHUNK_DIR="./chunks"
CHUNK_FILE="$CHUNK_DIR/target_chunk_$(printf "%03d" $SLURM_ARRAY_TASK_ID)"
MAPPING_FILE="./uniprot_to_pdb_id_map_2M.txt"
TMP_DAT_DIR="./tmp_dali_dat/tmp_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$TMP_DAT_DIR"

declare -A map
while read -r uniprot_name datname; do
    map["$uniprot_name"]="$datname"
done < "$MAPPING_FILE"

while read -r target; do
    ent_file="$UNIPROT_DIR/$target"
    out_file="${map[$target]}"

    if [[ -z "$out_file" ]]; then
        echo "[$SLURM_ARRAY_TASK_ID] No mapping for $target. Skipping."
        continue
    fi

    echo "[$SLURM_ARRAY_TASK_ID] Importing $target -> $out_file"

    rm -f "$TMP_DAT_DIR"/*.dat "$TMP_DAT_DIR/dali.lock"
    rm -f ./dali.lock  
    
    (
        cd "$TMP_DAT_DIR"
        "$IMPORT" --pdbfile "$ent_file" --pdbid "$out_file" --dat "." --clean
    )

    found_file=$(find "$TMP_DAT_DIR" -type f -name "${out_file}*.dat")

    if [[ -n "$found_file" ]]; then
        mv "$found_file" "$DAT_DIR/${out_file}A.dat"
        echo "[$SLURM_ARRAY_TASK_ID] Moved $found_file -> $DAT_DIR/$out_file"
    else
        echo "[$SLURM_ARRAY_TASK_ID] Warning: No .dat file found for $target"
    fi
done < "$CHUNK_FILE"

rm -rf "$TMP_DAT_DIR"

echo "[$SLURM_ARRAY_TASK_ID] Completed processing chunk $SLURM_ARRAY_TASK_ID"