./build/src/foldseek lolalign scope/target_db scope/target_db scope/allvsall_fakepref scope/allvsall_scope --threads 64 --lolalign-multidomain -0

./build/src/foldseek createtsv scope/target_db scope/target_db scope/allvsall_scope scope/lolalignaln.tsv

./bench.awk ../data/scop_lookup.fix.tsv <(cat ../alignResults/rawoutput/lolalignaln.tsv) > ../alignResults/rocx/lolalign.rocx

## calculate auc
 awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' ../alignResults/rocx/lolalign.rocx
