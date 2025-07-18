#!/bin/bash
## This is version 2.2
# - changed the telomere analysis to include telomere_analysis.sh from [vgp-pipeline](https://github.com/VGP/vgp-assembly). Basically, in 5k window, will look for occurence of "TTAGGG" and if the window contains 50% of this repeat, then it will count as telomere length. (on-going)
# - because of this, additional column in the *stat.txt for telomere length.
# - T2T chromosome standard will also change based on the telomere_analysis.sh result
## by PSPineda 2024.07.18 (polenpineda@gmail.com)

set -euo pipefail

helpFunction()
{
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
      r ) if [ -z "$OPTARG" ] ; then
      echo "Error: -r requires an argument"
      exit 1
      fi
      ref="$OPTARG" ;;
      q ) if [ -z "$OPTARG" ] ; then
      echo "Error: -q requires an argument"
      exit 1
      fi
      qry="$OPTARG" ;;
      o ) dirname="$OPTARG" ;;
      m ) map="$OPTARG" ;;
      p ) path="$OPTARG" ;;
      c ) tel_cutoff="$OPTARG" ;;
      t ) t="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

if [ ! $# -gt 3 ] ; then
    helpFunction
    exit 1
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
QC_dir=$(realpath $QC_dir)

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
  echo "Running minimap2 with the following command:"
  echo "minimap2 $ref $qry -t $threads -x asm5 -o $dir/minimap.paf"
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
samtools faidx "$outname"
seqtk cutN -n 3 -g "$outname" > $(basename "$outname" .fasta).coor

echo -e "Running telomere_analysis.sh" to estimate telomere length

export VGP_PIPELINE="$QC_dir/tools"
telomerebase=$(basename "$outname" .fasta)

echo -e "$VGP_PIPELINE/telomere/telomere_analysis.sh telomere_$telomerebase 0.5 5000 $dir/$outname"
"$VGP_PIPELINE"/telomere/telomere_analysis.sh telomere_"$telomerebase" 0.5 5000 "$dir"/"$outname" || true
# this is T2T Han telomere analysis parameter used

# echo -e "Running tidk to count telomeres."
# installed on conda activate centromere
# tidk search -s TTAGGG -o "$telomerebase"_bedgraph_tidk-search --dir "$dir" -e bedgraph "$dir"/"$outname" -w 1000
# tidk search -s TTAGGG -o "$telomerebase"_tsv_tidk-search --dir "$dir" -e tsv "$dir"/"$outname" -w 1000
# echo -e "Generating plot..."
# tidk plot -o "$telomerebase"_tidk-plot -t "$dir"/"$telomerebase"_tsv_tidk-search_telomeric_repeat_windows.tsv

# echo -e "Generating count for sex chromosomes."
# tidk search -s TTAGGG -o X_tsv_tidk-search --dir "$dir" -e tsv "$dir"/combine/X.fasta -w 1000
# tidk search -s TTAGGG -o Y_tsv_tidk-search --dir "$dir" -e tsv "$dir"/combine/Y.fasta -w 1000

# echo -e "Generating plot for sex chromosomes."

# x_size=$(stat -c %s "$dir"/X_tsv_tidk-search_telomeric_repeat_windows.tsv)
# y_size=$(stat -c %s "$dir"/Y_tsv_tidk-search_telomeric_repeat_windows.tsv)

# if [ "$x_size" -gt 4096 ]; then
# 	tidk plot -o X_tidk-plot -t "$dir"/X_tsv_tidk-search_telomeric_repeat_windows.tsv
# 	seqtk cutN -n 3 -g "$dir"/combine/X.fasta > "$dir"/X.coor
# else
# 	echo -e "No telomeres found for X"
# fi

# if [ "$y_size" -gt 4096 ]; then
# 	tidk plot -o Y_tidk-plot -t "$dir"/Y_tsv_tidk-search_telomeric_repeat_windows.tsv
# 	seqtk cutN -n 3 -g "$dir"/combine/Y.fasta > "$dir"/Y.coor
# else
# 	echo -e "No telomeres found for X"
# fi

echo -e "Extract the telomeric regions based on telomere_analysis.sh"

df1="$dir"/"$outname".fai
df2=$(find "telomere_$telomerebase" -name "*.windows.0.5.5kb.ends.bed")
output_bed="$dir/telomere_updated.bed"

# Read df1 into an associative array
declare -A df1_map
while read -r col1 col2 col3 col4 col5; do
    df1_map["$col1"]=$col2
done < "$df1"

while read -r col1 col2 col3; do
    if [ "$col2" != "0" ]; then
        if [ -n "${df1_map["$col1"]}" ]; then
            col3=${df1_map["$col1"]}
        fi
    fi
    echo -e "$col1\t$col2\t$col3" >> "$output_bed"
done < "$df2"


echo -e "Updated bed file of telomeric region. Saved to $output_bed."
echo -e "Extracting telomeres fasta..."

bedtools getfasta -fi $dir/$outname -bed $dir/telomere_updated.bed -fo $dir/telomere.fasta
"$QC_dir"/utils/count_telomere.sh $dir/telomere.fasta $dir/telomere.count.txt

echo -e "Running stat_moreinfo.R..."

echo "path file: $path"
echo "map file: $map"

if [ "$path" = "0" ] && [ "$map" = "0" ]; then
  echo -e "No nodes path and map input, will not add nodes pathway to the all_STATS.tsv file."
  Rscript "$QC_dir"/utils/stat_moreinfo_nopathmap.R "$dir" $(basename "$outname" .fasta).coor "$dir"/telomere.count.txt "$tel_cutoff" telomere_"$telomerebase"/*.windows.0.5.5kb.ends.bed
else
  echo -e "Nodes path and map provided, will add pathway to the all_STATS.tsv file."
  Rscript "$QC_dir"/utils/stat_moreinfo.R "$dir" "$path" "$map" $(basename "$outname" .fasta).coor "$dir"/telomere.count.txt "$tel_cutoff" telomere_"$telomerebase"/*.windows.0.5.5kb.ends.bed
fi

echo "Done. :)"
