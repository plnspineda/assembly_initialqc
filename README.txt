Installation:
Just copy the folder then

	cd assembly_initialqc
	export QC=$PWD

Dependencies (need this tools to run):
	minimap2 
	samtools
	R
	Java
	tidk
	CombineFasta (https://github.com/njdbickhart/CombineFasta)

Note: CombineFasta is hardcoded so change pathfile at line 50 of assembly_initialqc.sh.

Example run in slurm:

#!/bin/bash
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time=72:00:00
#SBATCH --mem=64GB

module purge
module use /apps/modules/all
module load R
module load SAMtools/1.17-GCC-11.2.0
module load minimap2/2.26-GCCcore-11.2.0
module load Java/17.0.6
#conda activate centromere

export QC="/hpcfs/groups/phoenix-hpc-avsci/Lloyd_Low/Tuli_x_Wagyu_data/scripts/assembly_initialqc"

ref="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/swamp_buffalo/reference_assembly/renamed_assemblies/ARS-UCD2.0_frozen.fna"
qry="/hpcfs/groups/phoenix-hpc-avsci/Lloyd_Low/Tuli_x_Wagyu_data/assembly_contig/verkko/TxWverkkoTrioT2T_gt7mat_asmonly/assembly.haplotype1.fasta"
dir="/hpcfs/groups/phoenix-hpc-avsci/Lloyd_Low/Tuli_x_Wagyu_data/assembly_contig/verkko/TxWverkkoTrioT2T_gt7mat_asmonly/assembly_qc/hap1"

$QC/assembly_initialqc.sh $ref $qry $dir
	
