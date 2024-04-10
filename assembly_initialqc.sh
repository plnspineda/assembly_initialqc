#!/bin/bash
## This is version 2.1
### changes: made map and path file optional so could run with hifiasm assembly as well
### changes: trying to add parameters... still commented as not yet tested
## by PSPineda 2024.03.27 (polenpineda@gmail.com)

set -euo pipefail

helpFunction()
{
   echo ""
   echo "Usage: $0 -r <ref.fa> -q <qry.fa> -o <output_directory> [-m <assembly.scfmap>] [-p <assembly.paths.tsv>] [-c <telomere_count>] [-t <num_threads>]"
   echo "Description: To initially find the contigs equivalent to a chromosome, count telomeres and output a summary statistics."
   echo -e "\t-r reference genome"
   echo -e "\t-q query genome"
   echo -e "\t-o output directory"
   echo -e "\t-m assembly.scfmap output from verkko (optional)"
   echo -e "\t-p assembly.paths.tsv output from verkko (optional)"
   echo -e "\t-c telomere count cutoff (default 50)"
   echo -e "\t-t threads (default 2)"
   exit 1 # Exit script after printing help
}

while getopts "r:q:o:m:p:c:t:" opt
do
   case "$opt" in
      r ) ref="$OPTARG" ;;
      q ) qry="$OPTARG" ;;
      o ) dirname="$OPTARG" ;;
      m ) map="$OPTARG" ;;
      p ) path="$OPTARG" ;;
      c ) tel_cutoff="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$ref" ]
then
   echo "Please input reference genome"
   helpFunction
elif [ -z "$qry" ]
then
   echo "Please input query genome"
   helpFunction
elif [ -z "$dirname" ]
then
   echo "Please input output name"
   helpFunction
fi

echo -e "Reading input"
ref="$(realpath "$ref")"
qry="$(realpath "$qry")"
map="${map:-"0"}"
path="${path:-"0"}"
tel_cutoff="${tel_cutoff:-50}"
threads="${t:-2}"
outname="$(basename "$qry" .fa*)_tmp_asm.fasta"

QC_dir=$(dirname $0)

# if [ $# -eq 0 ]; then
#     echo "Usage: $0 <reference.fasta> <query.fasta> <output_dir> <assembly.scfmap> <assembly.paths.tsv>"
#     echo "Description: To initially find the contigs equivalent to a chromosome, count telomeres and output a summary statistics."
#     exit 1
# fi

# set -euo pipefail

# echo -e "Reading input"
# ref="$(realpath $1)"
# qry="$(realpath $2)"
# outname="$(basename "$2" .fa*)_tmp_asm.fasta"
# dirname="$3"
# map="$(realpath "${4:-0}")"
# path="$(realpath "${5:-0}")"
# t="${6:-2}"
# tel_cutoff="${7:-50}"

mkdir -p "$dirname"
dir=$(realpath $dirname)
if [ "$map" != "0" ] && [ "$path" != "0" ]; then
    map="$(realpath "$map")"
    path="$(realpath "$path")"
fi

echo -e "Reference genome: $(basename "$ref")"
echo -e "Query genome: $(basename "$qry")"
echo -e "Output fasta: $(basename "$outname")"
echo -e "Output directory: $dir"
echo -e "Input map file: $(basename "$map")"
echo -e "Input path file: $(basename "$path")"
echo -e "Telomere count cutoff: "$tel_cutoff""
echo -e "Number of threads: "$threads""

echo -e "Mapping query genome to reference genome"

if [ ! -f "$dir/minimap.paf" ]; then
    minimap2 "$ref" "$qry" -t "$threads" -x asm5 -o "$dir/minimap.paf"
fi
mkdir -p "$dir"/contig_list
mkdir -p "$dir"/archived

echo -e "Running Rscript to get homologous chromosomes and Computing assembly and reference propotion."
# Rscript findContigs.R "$ref" "$qry"

Rscript --version
Rscript "$QC_dir"/utils/stat_initialqc.R "$dir"

echo -e "Making combine folder and getting the contigs"
mkdir -p "$dir"/combine
while IFS=$'\t' read -r f1 f2 f3 f4; do
  echo -e "Getting $f1 contig and saving in folder"
  samtools faidx "$qry" "$f1" > "$dir"/combine/contig_"$f1".fa
done < <(cat "$dir"/contig_list/*list)

echo -e "Creating list for combinefasta"
for j in "$dir"/contig_list/*.list; do
	jbase=$(basename "$j" .list)
	awk '{print "contig_"$1".fa	" $2}' "$j" > "$dir"/combine/"$jbase"_order.list
done

echo -e "Running CombineFasta to create homologous chromosome"
cd "$dir"/combine || exit
for k in *_order.list; do
  kbase=$(basename "$k" _order.list)
  fasta_file="$kbase.fasta"
  if [ ! -f "$fasta_file" ]; then
    echo -e "Combining contigs for chromosome $kbase with $k as the order list."
    java -jar "$QC_dir"/tools/CombineFasta/store/CombineFasta.jar order \
      -i "$kbase"_order.list \
      -o "$fasta_file" \
      -p 100 \
      -n "$kbase"
  else
    echo "Skipping, $fasta_file already exists."
  fi
done

echo -e "Concatenating chromosomes into one fasta"
cd "$dir/combine"
if [ -e "../$outname" ]; then
    echo "$outname exists. Skipping."
else
    ls ./*fasta | sort -V
    # cat "$(ls -v combine/*[0-9].fasta)" > "$outname"
    cat $(ls -v ./*[0-9].fasta) > ../"$outname"
    X=$(stat -c %s "X.fasta")
    Y=$(stat -c %s "Y.fasta")

    if [ $X -gt 4096 ] && [ $Y -gt 4 ]; then
      cat "X.fasta" "Y.fasta" >> ../"$outname"
    elif [ $X -gt 4096 ]; then
      cat "X.fasta" >> ../"$outname"
    elif [ $Y -gt 4096 ]; then
      cat "Y.fasta" >> ../"$outname"
    fi
fi

cd ../

echo -e "Running seqtk to get gap coordinates..."
seqtk cutN -n 3 -g "$outname" > $(basename "$outname" .fasta).coor

echo -e "Running tidk to count telomeres."
# installed on conda activate centromere

tidkbase=$(basename "$outname" .fasta)

tidk search -s TTAGGG -o "$tidkbase"_bedgraph_tidk-search --dir "$dir" -e bedgraph "$dir"/"$outname" -w 10000
tidk search -s TTAGGG -o "$tidkbase"_tsv_tidk-search --dir "$dir" -e tsv "$dir"/"$outname" -w 10000
echo -e "Generating plot..."
tidk plot -o "$tidkbase"_tidk-plot -t "$dir"/"$tidkbase"_tsv_tidk-search_telomeric_repeat_windows.tsv

echo -e "Generating count for sex chromosomes."
tidk search -s TTAGGG -o X_tsv_tidk-search --dir "$dir" -e tsv "$dir"/combine/X.fasta -w 10000
tidk search -s TTAGGG -o Y_tsv_tidk-search --dir "$dir" -e tsv "$dir"/combine/Y.fasta -w 10000

echo -e "Generating plot for sex chromosomes."

x_size=$(stat -c %s "$dir"/X_tsv_tidk-search_telomeric_repeat_windows.tsv)
y_size=$(stat -c %s "$dir"/Y_tsv_tidk-search_telomeric_repeat_windows.tsv)

if [ "$x_size" -gt 4096 ]; then
	tidk plot -o X_tidk-plot -t "$dir"/X_tsv_tidk-search_telomeric_repeat_windows.tsv
	seqtk cutN -n 3 -g "$dir"/combine/X.fasta > "$dir"/X.coor
else
	echo -e "No telomeres found for X"
fi

if [ "$y_size" -gt 4096 ]; then
	tidk plot -o Y_tidk-plot -t "$dir"/Y_tsv_tidk-search_telomeric_repeat_windows.tsv
	seqtk cutN -n 3 -g "$dir"/combine/Y.fasta > "$dir"/Y.coor
else
	echo -e "No telomeres found for X"
fi

echo -e "Running stat_moreinfo.R..."

echo "path file: $path"
echo "map file: $map"

if [ "$path" = "0" ] && [ "$map" = "0" ]; then
  echo -e "No nodes path and map input, will not add nodes pathway to the all_STATS.tsv file."
  Rscript "$QC_dir"/utils/stat_moreinfo_nopathmap.R "$dir" $(basename "$outname" .fasta).coor "$tidkbase"_bedgraph_tidk-search_telomeric_repeat_windows.bedgraph "$tel_cutoff"
else
  echo -e "Nodes path and map provided, will add pathway to the all_STATS.tsv file."
  Rscript "$QC_dir"/utils/stat_moreinfo.R "$dir" "$path" "$map" $(basename "$outname" .fasta).coor "$tidkbase"_bedgraph_tidk-search_telomeric_repeat_windows.bedgraph "$tel_cutoff"
fi

echo "Done. :)"
