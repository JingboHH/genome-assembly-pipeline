# De novo assembly pipeline for paired-end sequencing reads



## Preparation

Input files: paired-end reads sequenced by Illumina platform

Used software and tools:

Merging paired-end reads: [BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)

Assembling merged reads to generate contigs: [Unicycler](https://github.com/rrwick/Unicycler)

Annotating genes for contigs: [Prokka](https://github.com/tseemann/prokka)

Quality assessment of assemblies: [QUAST](https://github.com/ablab/quast)

## Merging paired-end reads

A pair of read may be split seperate into two folder: R1 and R2

The example of files name: DNA1_R1.fq.gz; DNA1_R2.fq.gz

```
READS_DIR_R1='reads/R1'
READS_DIR_R2='reads/R2'
```
Define the base folder and subfolders for all outputs

```
OUTPUT_DIR='/output'
MERGED_DIR="$OUTPUT_DIR/merged"
UNMERGED_DIR="$OUTPUT_DIR/unmerged"
IHIST_DIR="$OUTPUT_DIR/ihist"
```
To ensure output dirctories exist
```
mkdir -p $MERGED_DIR $UNMERGED_DIR $IHIST_DIR
```
Define a loop to extract the base name of all R1 files without the suffix of "_R1.fq.gz" for naming the output files and finding the corresponding R2

Define the name of output files

Run the script bbmerge.sh to merge the overlapping reads. Specify input files, output files

Log the completion of processing for each pair of reads, indicating the output files generated
```
for R1 in $READS_DIR/*_R1.fastq.gz; do
    BASENAME=$(basename $R1 _R1.fastq.gz)
    R2="$READS_DIR_R2/${BASENAME}_R2.fastq.gz"
    if [ ! -f "$R2" ]; then
        echo "Warning: Missing corresponding R2 file for $R1, skipping..."
        continue
    fi

OUT_MERGED="$MERGED_DIR/${BASENAME}_merged.fq.gz"
OUT_UNMERGED="$UNMERGED_DIR/${BASENAME}_unmerged.fq.gz"
OUT_IHIST="$IHIST_DIR/${BASENAME}_ihist.txt"

bbmerge.sh in1=$R1 in2=$R2 out=$OUT_MERGED outu=$OUT_UNMERGED ihist=$OUT_IHIST

echo "Processed $BASENAME: Merged -> $OUT_MERGED; Unmerged -> $OUT_UNMERGED; IHist -> $OUT_IHIST"
done
```
Additional options of BBMerge, attampting merge each pair of reads
If it is unsuccessful, then extend each read of 20bp and merge them again
```
bbmerge-auto.sh in1=$R1 in2=$R2 out=$OUT_MERGED outu=$OUT_UNMERGED ihist=$OUT_IHIST ecct extend2=20 iterations=5
```
## Assembling merged reads to generate contigs
Input: merged reads, for example, DNA1_merged.fq.gz

Output: assemblies, for example, assembly.fasta

Define the directory of output
```
UNICYCLER_OUTPUT_DIR='/unicycler_output'
```
Ensure the Unicycler output directory exists
```
mkdir -p "$UNICYCLER_OUTPUT_DIR"
```
Change directory to location of the merged files
```
cd /merged_files
```
Define a loop through all the merged reads, extract the base name without the suffix of ".fq.gz" to name the output

"-s" means run Unicycler in unpaired reads (merged reads) mode
```
for file in *.fq.gz; do
    base_name=$(basename "$file" .fq.gz)
    output_dir="$UNICYCLER_OUTPUT_DIR/${base_name}_unicycler_output"
    mkdir -p "$output_dir"
    unicycler -s "$file" -o "$output_dir"
    echo "Unicycler processing completed for $base_name"
done
```
Additional option: Unicycler can also run with paired-end reads directly
```
unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -o output_dir
```
## Annotating genes for contigs of assemblies
Input: assemblies, assembly.fasta
Output: annotated assemblies, DNA1_merged_

For each merged read, Unicycler outputs a folder that contains the assembly. For this reason, we define a base directory that contain all subfolders output from Unicycler
```
BASE_DIR="/assembled_folder"
```
Define Prokka annotation output
```
ANNOTATION_BASE_DIR="/annotation/output"
```
Ensure annotation directory exist
```
mkdir -p "$ANNOTATION_BASE_DIR"
```
Define a loop through each subfolder in the base directory, and construct the path to the assembly file within the base directory
```
for dir in "$BASE_DIR"/*/; do
    echo "Checking directory: $dir"
    assembly_file="${dir}assembly.fasta"
    echo "Looking for assembly file at: $assembly_file"
```
Check if the assembly.fasta exists, if exists, then extract the sample name from the directory path
```
    if [[ -f "$assembly_file" ]]; then
        sample_name=$(basename "$dir")
```
Define and create the output directory (folder name) for annotations based on the sample name
```
        outdir="$ANNOTATION_BASE_DIR/${sample_name}_annotation"
        echo "Creating output directory at: $outdir"
        mkdir -p "$outdir"
```
Run Prokka for annotation in reference genome priority mode with the reference genome ecoli_k12.gbff, specifying the input assembly file and output directory
```
        echo "Running Prokka for $sample_name"
        prokka --force --proteins "/reference_genome/ecoli_k12.gbff" --outdir "$outdir" --prefix "$sample_name" "$assembly_file"
        echo "Prokka annotation completed for $sample_name"
    else
        echo "assembly.fasta not found in $dir"
    fi
done
```
Additional option of Prokka, in default mode: prokka --outdir "$outdir" --prefix "$prefix" "${dir}/assembly.fasta"
## Assemblies quality assessment
## Comparing assemblies from different assemblers or pipelines
Input: assembly.fasta in each output subfolder from Unicycler
Output: QUAST analysis report
Define input and output directories
```
NEW_ASSEMBLY_DIR=/new_assemblies/
EXIST_ASSEMBLY_DIR_2=/exist_assemblies/
QUAST_OUTPUT_DIR=/quast_output/
```
Ensure QUAST output directory exists
```
mkdir -p "$QUAST_OUTPUT_DIR"
```
Define a loop to thorough all subfolders in NEW_ASSEMBLY_DIR
```
for dir in "$NEW_ASSEMBLY_DIR"/*/; do
    echo "Checking directory: $dir"
```
Extract the sample name from basename by removing the suffix "_merged_unicycler_output", and using sample name to define the name of new_assemly and exist_assembly, as well as output subdirectories
```
    sample_name=$(basename "$dir" | sed 's/_merged_unicycler_output//')
    new_assembly_file="${dir}assembly.fasta"
    exist_assembly_file="$EXIST_ASSEMBLY_DIR/${sample_name}.fasta"
    output_dir="$QUAST_OUTPUT_DIR/quast_results_${sample_name}"
```
After preparation of assemblies from two sources, then run the QUAST on them. "-r" means the reference genome for assembly quality, such as NGA50, or misassemblies metrics; while "-g" means the reference Genbank for gene assessment, such as the CDS number and gene number. "-t" defines the number of threads when running QUAST, "-1" means for single-thread: QUAST does not work with multi-threading when using the latest version of python, which is using python 3.7 and below.
```
    if [[ -f "$new_assembly_file" && -f "$exist_assembly_file" ]]; then
        quast.py -t 1 -G "/path/to/reference/genome/ecoli_k12.gbff" -o "$output_dir" "$new_assembly_file" "$exist_assembly_file"
        echo "QUAST analysis completed for $sample_name"
    else
        [[ ! -f "$new_assembly_file" ]] && echo "New assembly .gfa file not found in $dir"
        [[ ! -f "$exist_assembly_file" ]] && echo "Exist assembly .fasta file not found for $sample_name"
    fi
done
```
## Analyzing all assemblies from the same pipelines
MetaQUAST is another version of QUAST used for multiple assemblies, it generates several tables and graphs of different metrics for all assemblies
Locates to the directory of all assemblies
```
cd assemblies_folder
```
Run MetaQUAST for all assemblies
```
metaquast.py $(ls *.fasta)
-r /ecoliK12.fna
-g /ecoliK12.gff
-o output
-t 1

















