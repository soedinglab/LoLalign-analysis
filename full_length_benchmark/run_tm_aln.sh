#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

query_index=$SLURM_ARRAY_TASK_ID

mkdir -p ./output_folder ./results_folder logs_tm_run

# Get unique queries
mapfile -t queryArray < <(awk '{print $1}' ./e_bench_3di_pairs.txt | sort -u)

if [[ $query_index -ge ${#queryArray[@]} ]]; then
    echo "Error: query_index $query_index out of bounds" >&2
    exit 1
fi

query=${queryArray[$query_index]}
outfile=./results_folder/${query}_results.csv
echo "query,target,qstart,tstart,q_tm_score,t_tm_score,q_len,t_len,aln_len,RMSD,seq_id,cigar" > "$outfile"

mapfile -t targets < <(awk -v q="$query" '$1 == q {print $2}' ./e_bench_3di_pairs.txt)

echo "Processing query: $query (${#targets[@]} targets)"

for target in "${targets[@]}"; do
    [[ -z "$target" ]] && continue

    tmpout=./output_folder/tmp_${query}_${target}.out

    if [[ ! -f "/alphafold_data/${query}" ]] || \
       [[ ! -f "/alphafold_data/${target}" ]]; then
        echo "Warning: missing file for $query vs $target" >&2
        continue
    fi

    TMalign /alphafold_data/${query} \
            /alphafold_data/${target} \
            -outfmt 1 > "$tmpout" 2>&1 || {
        echo "TMalign failed for $query vs $target" >&2
        rm -f "$tmpout"
        continue
    }

    query_header=$(grep "^>.*${query}" "$tmpout" | head -1 || true)
    target_header=$(grep "^>.*${target}" "$tmpout" | head -1 || true)

    query_alignment=$(awk -v q="$query" '
        /^>/ {if(index($0,q)>0){getline; gsub(/[ \t\r\n]/,""); print; exit}}' "$tmpout")
    target_alignment=$(awk -v t="$target" '
        /^>/ {if(index($0,t)>0){getline; gsub(/[ \t\r\n]/,""); print; exit}}' "$tmpout")

    query_alignment=$(echo "$query_alignment" | sed 's/[^A-Za-z-]//g')
    target_alignment=$(echo "$target_alignment" | sed 's/[^A-Za-z-]//g')

    if [[ -z "$query_alignment" || -z "$target_alignment" ]]; then
        echo "Warning: missing alignment for $query vs $target" >&2
        rm -f "$tmpout"
        continue
    fi

    if [[ ${#query_alignment} -ne ${#target_alignment} ]]; then
        echo "Warning: alignment lengths differ for $query vs $target" >&2
        rm -f "$tmpout"
        continue
    fi

    query_tm_score=$(echo "$query_header" | sed -n 's/.*TM-score=\([0-9.]*\).*/\1/p')
    target_tm_score=$(echo "$target_header" | sed -n 's/.*TM-score=\([0-9.]*\).*/\1/p')
    query_len=$(echo "$query_header" | sed -n 's/.*L=\([0-9]*\).*/\1/p')
    target_len=$(echo "$target_header" | sed -n 's/.*L=\([0-9]*\).*/\1/p')
    seq_id=$(echo "$query_header" | sed -n 's/.*seqID=\([0-9.]*\).*/\1/p')

    align_len=$(grep "^# Lali=" "$tmpout" | awk '{for(i=1;i<=NF;i++) if($i~/Lali=/){split($i,a,"="); print a[2]; exit}}')
    rmsd=$(grep "^# Lali=" "$tmpout" | awk '{for(i=1;i<=NF;i++) if($i~/RMSD=/){split($i,a,"="); print a[2]; exit}}')

    # Compute start positions & CIGAR
    read -r qstart tstart cigar <<< "$(python3 - <<PY
q = """${query_alignment}"""
t = """${target_alignment}"""
L = len(q)
left = None; right = None
for i in range(L):
    if q[i] != '-' and t[i] != '-':
        left = i; break
for i in range(L-1,-1,-1):
    if q[i] != '-' and t[i] != '-':
        right = i; break
if left is None or right is None:
    exit(1)
qstart = sum(1 for c in q[:left] if c != '-') + 1
tstart = sum(1 for c in t[:left] if c != '-') + 1
cigar = []
cur, count = None, 0
for a,b in zip(q[left:right+1], t[left:right+1]):
    if a == '-': op = 'D'
    elif b == '-': op = 'I'
    else: op = 'M'
    if cur == op: count += 1
    else:
        if cur: cigar.append(f"{count}{cur}")
        cur, count = op, 1
if cur: cigar.append(f"{count}{cur}")
print(f"{qstart}\t{tstart}\t{''.join(cigar)}")
PY
)"

    if [[ -z "$cigar" ]]; then
        echo "Warning: failed CIGAR for $query vs $target" >&2
        rm -f "$tmpout"
        continue
    fi

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$query" "$target" "$qstart" "$tstart" \
        "${query_tm_score:-}" "${target_tm_score:-}" \
        "${query_len:-}" "${target_len:-}" \
        "${align_len:-}" "${rmsd:-}" \
        "${seq_id:-}" "${cigar}" >> "$outfile"

    rm -f "$tmpout"
done

echo "Completed processing query: $query (${#targets[@]} targets checked)"
