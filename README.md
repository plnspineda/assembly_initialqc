# assembly_initialqc
Script for checking homologous chromosomes with a reference genome and checking gap and telomeric sequence presence after contig assembly. 

## REQUIREMENTS

- [R](https://www.r-project.org)
- [samtools](http://www.htslib.org)
- [minimap2](https://github.com/lh3/minimap2)
- [tidk](https://github.com/tolkit/telomeric-identifier)
- [seqtk](https://github.com/lh3/seqtk)

## Installation

1. Clone the repository:

        git clone https://github.com/plnspineda/assembly_initialqc.git
        cd assembly_initialqc

2. Create a conda environment:

        conda env create --file initialqc-env.yml

3. Run

        export QC="$PWD"
        $QC/assembly_initialqc.sh

## Usage

    $QC/assembly_initialqc.sh ref.fa asm.fa outdir assembly.scfmap assembly.paths.tsv 

wherein:

- `ref.fa` is the reference fasta genome
- `asm.fa` is the draft contig assembly
- `outdir` is the output directory (will make one if the folder do not exists)
- `assembly.scfmap` is an output from verkko along with `assembly.paths`. These are used to get the pathname and nodes for each chromosome/contig.

## Expected output files

file:

    - all_chr_assembly.tsv (incomplete information about each homologous contigs, can be deleted)
    - all_STATS.tsv (final output table with all information)
    - *_tmp_asm_bedgraph_tidk-search_telomeric_repeat_windows.bedgraph (bedgraph of the telomeres for autosomes)
    - *.coor (gaps coordinate)
    - *_tmp_asm.fasta (autosomal fasta file)
    - *.svg (telomeric signal plot)
    - *search_telomeric_repeat_windows.tsv (telomeric signal per 10kb window)
    - minimap.paf (alignment file)

folders:

    - archived (contains detailed information of the alignment per chromosomes)
    - contig_list (contains lists of contigs, orientation and proportion with the reference genome)
    - combine (contains fasta file of each chromosomes)

*still under development*
