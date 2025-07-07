#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <output_file>"
    exit 1
fi

fasta_file=$1
output_file=$2

if [ ! -f "$fasta_file" ]; then
    echo "File $fasta_file not found!"
    exit 1
fi

sequence=""
sequence_id=""

echo -e "chr\tp\tq" > "$output_file"

while IFS= read -r line; do
    if [[ $line == \>* ]]; then
        if [ -n "$sequence" ]; then
            count_TTAGGG=$(echo "$sequence" | grep -o "TTAGGG" | wc -l)
            count_CCCTAA=$(echo "$sequence" | grep -o "CCCTAA" | wc -l)
            echo -e "$sequence_id\t$count_CCCTAA\t$count_TTAGGG" >> "$output_file"
            sequence=""
        fi

        # Extract the sequence ID from the header
        sequence_id=$(echo "$line" | sed 's/^>//')
    else
        sequence="${sequence}${line}"
    fi
done < "$fasta_file"

if [ -n "$sequence" ]; then
    count_TTAGGG=$(echo "$sequence" | grep -o "TTAGGG" | wc -l)
    count_CCCTAA=$(echo "$sequence" | grep -o "CCCTAA" | wc -l)

    echo -e "$sequence_id\t$count_TTAGGG\t$count_CCCTAA" >> "$output_file"
fi

echo "Results saved to $output_file"
