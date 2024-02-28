#!/bin/bash

gene_idlist=()
#gene_file="./all_blast_1e-10genelist.txt"
#gene_file="./01all_vinif_blast_1e-10genelist.txt"
gene_file="./01all_vinif_cold_1e-10_50genelist.txt"
gff_file="../pep_dir/GCF_000003745.3_12X_genomic.gff"

while IFS= read -r pep_id
do
  echo $pep_id
  gene_id=$(awk '$9 ~ /'"$pep_id"'/ {match($0, /;gene=([^;]+)/, str); print str[1]}' "../pep_dir/GCF_000003745.3_12X_genomic.gff" | head -n 1)
  echo "$gene_id"
  gene_idlist+=("$gene_id")
done < "$gene_file"

# 将数组内容写入txt文件
printf '%s\n' "${gene_idlist[@]}" > 01allcold_1e-10_50gene_ids.txt



