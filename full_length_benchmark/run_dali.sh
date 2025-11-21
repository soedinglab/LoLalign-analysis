#!/bin/bash


ALIGNER_RESULTS="$SLURM_SUBMIT_DIR/3di_pairs.txt"
QUERY_LIST="$SLURM_SUBMIT_DIR/rep_ids_query.txt"
UNIPROT_MAP="$SLURM_SUBMIT_DIR/uniprot_to_pdb_id_map_2M.txt"
DAT_DIR="/dali_files/dat_files"
ALIGNMENTS_DIR="/dali_files/results"
DALI_BIN="$(realpath ../../../pdb/DaliLite.v5/bin/dali.pl)"

# Optional parameter
USE_FILTERED_TARGETS=true  # Set to false to use full target list instead of aligner hits

mkdir -p "$ALIGNMENTS_DIR"
mkdir -p "$SLURM_SUBMIT_DIR/logs_dali_run"

if [ ! -f "$QUERY_LIST" ]; then
  echo "ERROR: Query list file not found: $QUERY_LIST"
  exit 1
fi

if [ ! -f "$UNIPROT_MAP" ]; then
  echo "ERROR: mapping file not found: $UNIPROT_MAP"
  exit 1
fi

if [ ! -f "$DALI_BIN" ]; then
  echo "ERROR: DALI binary not found: $DALI_BIN"
  exit 1
fi

UNIPROT_QUERY=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$QUERY_LIST")

if [ -z "$UNIPROT_QUERY" ]; then
  echo "ERROR: No query found at index $SLURM_ARRAY_TASK_ID"
  exit 1
fi

echo "Query ID: $UNIPROT_QUERY"

PDB_QUERY=$(awk -v uniprot="$UNIPROT_QUERY" '$1 == uniprot {print $2}' "$UNIPROT_MAP")

if [ -z "$PDB_QUERY" ]; then
  echo "ERROR: No PDB mapping found for ID: $UNIPROT_QUERY"
  exit 1
fi

echo "Mapped PDB Query ID: $PDB_QUERY"

WORKDIR="${ALIGNMENTS_DIR}/${UNIPROT_QUERY}_${PDB_QUERY}_run"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

create_target_list() {
    local query=$1
    local target_file="$WORKDIR/targets_for_${query}.txt"
    
    echo "Extracting targets for query $query from aligner results..." >&2
    
    awk -v query="$query" '$1 == query {print $2}' "$ALIGNER_RESULTS" > "$target_file"
    
    local target_count=$(wc -l < "$target_file")
    echo "Found $target_count targets for query $query" >&2
    
    if [ $target_count -eq 0 ]; then
        echo "WARNING: No targets found for query $query in aligner results" >&2
        return 1
    fi
    
    echo "Using $target_count targets for DALI alignment" >&2
    echo "$target_file"
}

convert_targets_to_pdb() {
    local uniprot_targets_file=$1
    local pdb_targets_file="$WORKDIR/pdb_targets_for_${UNIPROT_QUERY}.txt"
    
    echo "Converting target IDs to PDB IDs..." >&2
    echo "Reading from: $uniprot_targets_file" >&2
    echo "Writing to: $pdb_targets_file" >&2
    
    if [ ! -f "$uniprot_targets_file" ]; then
        echo "ERROR: targets file not found: $uniprot_targets_file" >&2
        return 1
    fi
    
    rm -f "$pdb_targets_file"
    
    while IFS= read -r uniprot_target; do
        if [ -z "$uniprot_target" ]; then
            continue
        fi
        
        pdb_target=$(awk -v uniprot="$uniprot_target" '$1 == uniprot {print $2}' "$UNIPROT_MAP")
        if [ -n "$pdb_target" ]; then
            echo "${pdb_target}A" >> "$pdb_targets_file"
        else
            echo "WARNING: No PDB mapping found for target ID: $uniprot_target" >&2
        fi
    done < "$uniprot_targets_file"
    
    if [ ! -f "$pdb_targets_file" ]; then
        echo "ERROR: No valid PDB targets found after conversion" >&2
        return 1
    fi
    
    local pdb_count=$(wc -l < "$pdb_targets_file")
    echo "Successfully converted to $pdb_count PDB targets" >&2
    
    echo "$pdb_targets_file"
}


UNIPROT_TARGETS_FILE=$(create_target_list "$UNIPROT_QUERY")
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create target list for $UNIPROT_QUERY"
    exit 1
fi

PDB_TARGETS_FILE=$(convert_targets_to_pdb "$UNIPROT_TARGETS_FILE")
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to convert targets to PDB format"
    exit 1
fi

echo "Using target list: $PDB_TARGETS_FILE"
echo "Number of targets: $(wc -l < "$PDB_TARGETS_FILE")"

echo "Running DALI for query: $UNIPROT_QUERY -> $PDB_QUERY"

# Run DALI with PDB IDs
"$DALI_BIN" \
  --cd1 "${PDB_QUERY}A" \
  --db "$PDB_TARGETS_FILE" \
  --dat1 "$DAT_DIR" \
  --dat2 "$DAT_DIR" \
  > dali.stdout.log 2> dali.stderr.log

if [ $? -eq 0 ]; then
    echo "DALI job completed successfully for query: $UNIPROT_QUERY -> $PDB_QUERY"
else
    echo "ERROR: DALI job failed for query: $UNIPROT_QUERY -> $PDB_QUERY"
    echo "Check dali.stderr.log for details"
    exit 1
fi

wolf_file="$WORKDIR/wolf_output"
filter_file="$WORKDIR/filter_input"
if [[ -f "$filter_file" ]]; then
    echo "Deleting $filter_file"
    rm -f "$filter_file"
else
    echo "$filter_file not found"
fi

if [[ -f "$wolf_file" ]]; then
    echo "Deleting $wolf_file"
    rm -f "$wolf_file"
else
    echo "$wolf_file not found"
fi

echo "DALI job complete for query: $UNIPROT_QUERY"