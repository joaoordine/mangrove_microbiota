# Evaluating the bioremediation potential of mangrove microbiota in polluted coastal ecosystems
We aimed to characterize, for the first time, the microbial communities of two highly impacted mangrove fragments in São Sebastião (São Paulo, Brazil), namely Araçá and Colhereiro mangroves. For this, we conducted microbiota analyses of sediments and whole-genome sequencing of isolated native bacteria. Additionally, we evaluated the biodegradation capacity of this bacteria against an urban landfill leachate, a common pollutant in mangrove fragments affected by human waste in urban centers. Results were validated by chemical analysis of the sediments, providing initial support for conserving these ecosystems. In this page, you can find all codes used during the development of the project. 
- Manuscript pre-print: (replace later with article doi)

# *1. Microbiota analysis (Nanopore Sequencing) pipeline*
- Code adapted from Zaramela and Zorz (https://github.com/Microbial-Ecosystems-Lab/asthma_obesity and https://github.com/jkzorz/SituSeq/). 

## R advice 
1. Open tidyverse environment (using conda), then open R in terminal 
2. Install required packages
3. When finished, open Rstudio, also in terminal. 
This workflow usually works better than simply oppening Rstudio 

# Preprocessing 
## Appending `Sample´ in all fastq files to get a pattern
```
for file in ID*fastq; do mv "$file" "Sample-$file"; done
```

## Quality control - before trimming/clipping
```
mkdir fastqc_before
fastqc -o fastqc_before *.fastq
multiqc fastqc_before
mv multiqc_report.html multiqc_report_before.html
```

Check script 01.Preprocessing.ipynb 

## Quality control - after trimming/clipping 
```
mkdir fastqc_after # I created this directory inside the `filtered`one 
fastqc -o fastqc_after *.fastq
multiqc fastqc_after
mv multiqc_report.html multiqc_report_after.html
```

# Custom database creation
Using Refseq NCBI 16S sequences 

## Download sequences' ID from NCBI - using edirect 
OBS.: Check IDs_16Ssequences.txt file to avoid running this part of the code

```
split -l 1000 IDs_16Ssequences.txt # spliting each line # if you run the following codes with the full IDs list, it`ll crash, so we're splitting this into separate files to make it easier
for i in ls x; do epost -db nuccore -input ${i} -format acc | efetch -format fasta > ${i}.fasta; done # used to post the IDs to the Entrez system (NCBI nuccore database). This step prepares the IDs for fetching; all files that were separated with the previous command start with 'x' so we use this to run the code in each file
cat .fasta > 16Ssequences.fasta # fetch the sequences associated with the posted IDs # concatenate all files into one
```

## Downloading taxonomical ID database - Refseq NCBI 16S sequences - it contains all hierarchical levels for plotting  
```
mkdir data
cd data
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -zxvf new_taxdump.tar.gz 
```

## Download the Refseq file available from Zaramela`s code and moved it to data directory
The reference file was obtained by parsing gff3 file provided by NCBI.
  
```
mv ~/Downloads/RefseqTaxID.txt .
```

Check script 02.Database_creation.ipynb

## Reference alignment 
KMA v1.4.14 - Clausen, Aarestrup & Lund. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics, 19, 307 (2018)

### Sending all my fastq files - filtered- to my data file 
```
conda install bioconda::kma
mkdir 16sequences_index_kma
mkdir output_16s_refseq
mkdir ../data/filtered_fastq | mv .fastq ../data/filtered_fastq 
```
### Create the databases needed to run KMA from a list of FASTA files
```
kma_index -i 16Ssequences.fasta -o 16sequences_index_kma/16Ssequences 
```

### Map and/or align raw reads to a template database created using kma_index
```
for file in filtered_fastq/.fastq; do
    filename=$(basename "$file")  # Extract the filename
    filename_no_ext="${filename%.}"  # Remove the file extension
    kma -i "$file" -o "output_16s_refseq/${filename_no_ext}_kma" -t_db 16sequences_index_kma/16Ssequences -bcNano -bc 0.7
done 
```

### Unzipping frag files to use as input for taxonomical assignation
```
gunzip Sample 
mv .frag ..
```

## Calculate annotation frequency 
I.e. the percentage of reads that were annotated using that database
Check script 03.Calc_annotation_freq.ipynb

# Taxonomy assignation to nanopore reads  
Check script 04.Taxonomy_assign.ipynb

# Diversity analysis
Check script 05.Diversity_analysis.ipynb

# Visualize taxonomy 
Check script 06.Visualize_taxonomy.ipynb

# Correlation analysis 
Check script 07.Correlation_heatmap.ipynb

# *2. Isolated bacteria - WGS (Illumina) pipeline*
Some mangrove isolates were selected for whole-genome sequencing (Illumina NovaSeq 6000) due to their features of interest.
- Pipeline adapted from Borelli, et al. "Combining functional genomics and whole-genome sequencing to detect antibiotic resistance genes in bacterial strains co-occurring simultaneously in a Brazilian hospital." Antibiotics 10.4 (2021): 419.

## Decompress raw files sent from sequencing company
```
tar xvzf NovogeneFeb24.tar.gz # decompressing files
mkdir mangrove_Genomes
```

```
for dir in */; do
    echo "Extracting files from $dir"
    cd "$dir"
    gunzip *.fq.gz
    cd ..
done # decompressing files, did this loop once for mangrove and once for valeria's directories 
```

```
mkdir all_fastq
for dir in */; do
    echo "Sending file copy from $dir"
    cd "$dir"
    scp *.fq ../all_fastq
    cd ..
done
```

## Installing required packages for analyses
```
conda create --name genome-analysis
 conda install bioconda::trimmomatic
 conda install bioconda::flash
 conda install bioconda::spades   
 conda install bioconda::fastqc
```

### For those packages with conflict, I installed them in individual environments
```
conda deactivate
conda create --name prokka-genome
conda activate prokka-genome
conda install -c conda-forge -c bioconda -c defaults prokka # conda install of Prokka can be sensitive to channel order. Installing like this tends to work better 
```
```
conda create --name bowtie-genome
conda activate bowtie-genome
conda install bioconda::bowtie2
```
```
conda create --name checkm-genome
conda activate checkm-genome 
conda install bioconda::checkm-genome
conda deactivate 
```
```
conda create --name multiqc
conda activate multiqc
conda install bioconda::multiqc
```

## Quality check of raw data before trimming/clipping
```
mkdir fastqc_before # Create a directory to store FastQC results before trimming
```
```
conda activate genome-analysis
for dir in */; do
    echo "Running FastQC in directory: $dir"
    fastqc -o fastqc_before "${dir}"*.fq
done # Run FastQC on all .fq files in each subdirectory
```
```
conda deactivate 
conda activate multiqc
multiqc fastqc_before -p -f # Generate a multiqc report for all FastQC results
mv multiqc_report.html fastqc_before/multiqc_report_before.html # move the file to the fastq directory and rename it 
mv multiqc_* fastqc_before/ 
```

## Processing raw sequencing files
### Removing adapters and low quality reads with Trimmomatic v0.39 
After checking the multiqc report, there's no need to remove adapters or to trim the reads (we already paid for for post-processing), but the general code for doing so is as follows:
```
java -jar /home/guilherme/miniconda3/envs/genome-analysis/share/trimmomatic-0.39-2/trimmomatic.jar PE -trimlog trimmomatic_mangrove.log  \
-basein SA1_CKDN240001995-1A_H2YCMDSXC_L2_ \
-baseout SA1_output_paired SA1_output_unpaired \
ILLUMINACLIP: \
LEADING: \
TRAILING: \
SLIDINGWINDOW: \
MINLEN: 
```

### Filter out possible human contamination from genomes 
#### Download the Human Reference Genome
```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz 
```
#### Create bowtie2 directory to store its files
```
mkdir bowtie_stuff 
mv GRCh38_latest_genomic.fna.gz bowtie_stuff
gunzip GRCh38_latest_genomic.fna.gz # extracting the fasta file 
```
```
mkdir all_fastq
find . -type f -name "*.fq" -exec scp {} all_fastq \; # I'm sending a copy of all fastq files to this subdirectory just to make things easier with bowtie2
```
#### Index the Human Reference Genome:
```
bowtie2-build GRCh38_latest_genomic.fna human_index
```
#### Align Reads Against the Human Reference Genome
```
mkdir bowtie_stuff/uncontaminated_reads # to store reads without human contamination
bowtie2 -x ~/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/bowtie_stuff/human_index -1 ~/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/all_fastq/*L2_1.fq -2 ~/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/all_fastq/*L2_2.fq -S human_aligned.sam --un bowtie_stuff/uncontaminated_reads -D 20 -R 3 -N 1 -L 20 --very-sensitive-local # 0.04% overall alignment rate, i.e. no contamination detected 
```

## Merging overlapped reads 
```
for file in *.fq; do
    filename=$(basename -- "$file") # Extract the filename without extensions
    filename_no_ext="${filename%-*}"
    flash "$filename_no_ext"-1A_H2YCMDSXC_L2_1.fq "$filename_no_ext"-1A_H2YCMDSXC_L2_2.fq -o "$filename_no_ext" -M 85
done
```
```
mkdir flash_output
mv *.fastq flash_output/
mv *.hist flash_output/
mv *.histogram flash_output/
```
Since I'm using the server, I can't visualize anything, so it's best to compact all histograms and send them to my computer so I can visualize them
```
zip histograms.zip *.histogram 
scp guilherme@143.107.194.159:/home/guilherme/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/all_fastq/flash_output/histograms.zip . 
```
On your own personal computer, run this  
```
mkdir hist
mv *.hist hist/
mv *.histogram hist/ # just cleaning up this directory a little bit more 
```

## Genome assembly 

### Facing AttributeError: module 'collections' has no attribute 'Hashable' when running Spades
When facing this error while using spades, run the following:
```
sed -i 's/collections\.Hashable/collections\.abc\.Hashable/g' /home/guilherme/miniconda3/envs/genome-analysis/share/spades-3.13.0-0/share/spades/pyyaml3/constructor.py
```

### Checking Spades installation before running it 
```
~/miniconda3/envs/genome-analysis/bin/spades.py --test # replace the pathway to the bin where you installed spades 
```
The following message should appear  
```
========= TEST PASSED CORRECTLY.
```

### Running Spades genome assembler 
```
mkdir assembled_genomes
```
```
input="/home/guilherme/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/all_fastq/flash_output"
for file in "$input"/*.fastq; do
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*.*}" 
    
    assembly_dir="/home/guilherme/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/assembled_genomes/${filename_no_ext}_assembly"
    mkdir -p "$assembly_dir" # Create a directory for the current genome assembly to store this output 
    
    ~/miniconda3/envs/genome-analysis/bin/spades.py -o "$assembly_dir" --merged "$filename_no_ext.extendedFrags.fastq" -s "$filename_no_ext.notCombined_1.fastq" -s "$filename_no_ext.notCombined_2.fastq"
done
```

### Running plasmidSPAdes 
To recover only plasmids from WGS 
```
mkdir plasmidSPAdes_output
```
```
input="/home/guilherme/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/all_fastq/flash_output"
for file in "$input"/*.fastq; do
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*.*}"

    assembly_dir="/home/guilherme/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/plasmidSPAdes_output/${filename_no_ext}_plasmid_assembly"
    mkdir -p "$assembly_dir" # Create a directory for the current genome assembly to store this output 
    
    ~/miniconda3/envs/genome-analysis/bin/spades.py -o "$assembly_dir" --merged "$filename_no_ext.extendedFrags.fastq" -s "$filename_no_ext.notCombined_1.fastq" -s "$filename_no_ext.notCombined_2.fastq" --plasmid
done
```

## Renaming scaffold files so they contain the isolate code in their name and are located un a different directory altogether 

### For genome assemblies
```
mkdir all_scaffolds
for dir in *_assembly/; do
    echo "Renaming/moving scaffold file from $dir"
    isolate=$(basename -- "$dir" | sed 's/_.*//') 
    cd "$dir"
    mv scaffolds.fasta ../all_scaffolds/"$isolate"_scaffold.fasta 
    cd ..
done
```

### For plasmid assemblies
```
mkdir all_scaffolds
for dir in *_plasmid_assembly/; do
    echo "Renaming/moving scaffold file from $dir"
    isolate=$(basename -- "$dir" | sed 's/_.*//') 
    cd "$dir"
    scp scaffolds.fasta ../all_scaffolds/"$isolate"_plasmids.fasta 
    cd ..
done
```

## Quality control of assembled genomes 
### with QUAST 
Produces better plots of several quality metrics compared to checkM (as well as using way less memory), but doesn't inform explicitly about % completeness/contamination
```
~/guilherme/miniconda3/envs/quast/bin/quast.py *_scaffold.fasta # for genomes
```
```
~/guilherme/miniconda3/envs/quast/bin/quast.py *_plasmids.fasta # for plasmids 
mv quast_results/ quast_results_plasmids
```

## Genome annotation 
```
for file in ./*_scaffold.fasta; do
    filename=$(basename -- "$file")
    filename_no_ext="${filename%_scaffold.fasta}"
    prokka "$file" --outdir "${filename_no_ext}_prokka_out" --force
done
```
### for plasmids
```
for file in ./*_plasmids.fasta; do
    prokka "$file" --outdir "${file}_prokkaOUT" --force --centre X --compliant
done
```

### Processing prokka output 
```
mkdir all_fna_genomes
mkdir all_gff_genomes
```
```
outdir="/home/guilherme/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/assembled_genomes/all_scaffolds/all_fna_genomes" # specifying the output directory
for dir in *_prokka_out/ ; do
    echo "Renaming/moving fna file from $dir"
    dirname=$(basename -- "$dir")
    scp "$dir"/PROKKA_03062024.fna "$outdir"/"$dirname".fna 
done # copying all .fna files to a separated directory 
```
```
outdir="/home/guilherme/Documents/novogene_genomes_feb24/NovogeneFeb2024/usftp21.novogene.com/01.RawData/mangrove_Genomes/assembled_genomes/all_scaffolds/all_gff_genomes" # specifying the output directory
for dir in *_prokka_out/ ; do
    echo "Renaming/moving fna file from $dir"
    dirname=$(basename -- "$dir")
    scp "$dir"/PROKKA_03062024.gff "$outdir"/"$dirname".gff 
done # copying all .gff files to a separated directory 
```

## Functional annotation 
### ABRICATE
#### Setting up environment
```
conda create --name abricate
conda install -c conda-forge -c bioconda -c defaults abricate
abricate --check
abricate --list
```

### Running abricate on different databases
NCBI database
```
abricate --db ncbi *_scaffold.fasta > ARGs-ncbi_genomes.faa # running abricate using NCBI database
abricate --summary ARGs-ncbi_genomes.faa > sum_ARGs_genomes.faa # create boolean table of ARG coverage on genomes
```
```
mkdir abricate_output
mv *.faa abricate_output/
```
Virulence genes 
```
abricate --db vfdb ../*_scaffold.fasta > ./ARGs-vfdb_genomes.faa # I didn't summarize the results because nothing really relevant was identified 
```
Plasmids 
```
abricate --db plasmidfinder ../*_scaffold.fasta > ./ARGs-plasmids_genomes.faa
```
CARD database
```
abricate --db card ../*_scaffold.fasta > ./ARGs-card_genomes.faa
abricate --summary ARGs-card_genomes.faa > sumCARD_ARGs_genomes.faa
```

### Annotation and analysis of secondary metabolite biosynthesis gene clusters 
AntiSMASH web server (default parameters): https://antismash.secondarymetabolites.org/#!/start 

## Taxonomical assignation of sequenced genomes
TYGS web server (default parameters): https://tygs.dsmz.de

## Submiting assembled genomes - NCBI Genome DB 
- BioProject PRJNA1042624 
When trying to submit the assembled genomes on Genome - NCBI database, I had to further trim the scaffolds to avoid issues. Along with that, I also used an NCBI tool to check for the presence of Illumina adaptors still present in the assembly.

### Filtering contigs with less than 200 nt
Check script 08.Filt_small_contig.ipynb 

### Removing contaminations identified by NCBI
FCS (Foreign Contamination Screening)-adaptor tool from NCBI - Astashyn, A., Tvedte, E.S., Sweeney, D. et al. Rapid and sensitive detection of genome contamination at scale with FCS-GX. Genome Biol 25, 60 (2024). https://doi.org/10.1186/s13059-024-03198-7

#### Creating directories to install tool
```
mkdir fcsadaptor
cd fcsadaptor
curl -LO https://github.com/ncbi/fcs/raw/main/dist/run_fcsadaptor.sh # Get a copy of the run script
chmod 755 run_fcsadaptor.sh # Change the permissions of run_fcsadaptor.sh
mkdir inputdir outputdir # Create input and output directories; move all fasta files inside inputdir
cd ~/Documents/IC3/Collaborations/Mangrove_Microbiota/WGS/NCBI_submission/fcsadaptor
```

#### Set singularity 
```
conda create --name singularity
conda activate singularity
conda install conda-forge::singularity
```

#### Download the Singularity image
```
curl https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-adaptor.sif -Lo fcs-adaptor.sif 
```

#### Run FCS-adaptor screen 
```
for file in inputdir/*.fasta; do
    filename=$(basename "$file")  # Extract the filename
    filename_no_ext="${filename%.fasta}"  # Remove the file extension
    mkdir -p "outputdir/$filename_no_ext" # Create the output directory for this file

    ./run_fcsadaptor.sh --fasta-input "$file" --output-dir "outputdir/$filename_no_ext" --prok --container-engine singularity --image fcs-adaptor.sif # Run the FCS adaptor tool on this FASTA file using Singularity
done # if the genomes are eukaryotic use --euk
```

### Alternative to remove contaminations identified by NCBI - manually (works better than the code above)

#### Download all contaminations reports from NCBI and move them to a specified folder
Trimming contamination reports 
```
for file in ./Contamination*; do
    tail -n 3 "$file" > "${file}.temp"
    mv "${file}.temp" "$file"
done
```
```
mkdir filtered_scaffolds
mv *_filtered.fasta filtered_scaffolds/
mv all_scaffolds/filtered_scaffolds NCBI_submission/
```

Removing contaminations 
#### For SA1
```
contamination_line=$(grep -n "NODE_49" filtered_scaffolds/SA1_scaffold_filtered.fasta | cut -d ":" -f 1)
next_line=$((contamination_line + 1)) # Identify the line numbers containing the contamination and the line below it
sed -e "${contamination_line},${next_line}d" filtered_scaffolds/SA1_scaffold_filtered.fasta > cleaned_SA1_scaffold_filtered.fasta # Exclude contamination 
grep NODE_49 cleaned_SA1_scaffold_filtered.fasta # checking if it worked
```

#### For the remaining isolates
```
contamination_line=$(grep -n "NODE_47" filtered_scaffolds/SB1_scaffold_filtered.fasta | cut -d ":" -f 1)
next_line=$((contamination_line + 1)) 
sed -e "${contamination_line},${next_line}d" filtered_scaffolds/SB1_scaffold_filtered.fasta > cleaned_SB1_scaffold_filtered.fasta 
grep NODE_47 cleaned_SB1_scaffold_filtered.fasta 
```
```
contamination_line=$(grep -n "NODE_26_length_5051_cov_287.398874" filtered_scaffolds/SB8_scaffold_filtered.fasta | cut -d ":" -f 1)
next_line=$((contamination_line + 1)) 
sed -e "${contamination_line},${next_line}d" filtered_scaffolds/SB8_scaffold_filtered.fasta > cleaned_SB8_scaffold_filtered.fasta 
grep NODE_26_length_5051_cov_287.398874 cleaned_SB8_scaffold_filtered.fasta 
```
```
contamination_line=$(grep -n "NODE_65_length_416_cov_233.578171" filtered_scaffolds/WB3_scaffold_filtered.fasta | cut -d ":" -f 1)
next_line=$((contamination_line + 1)) 
sed -e "${contamination_line},${next_line}d" filtered_scaffolds/WB3_scaffold_filtered.fasta > cleaned_WB3_scaffold_filtered.fasta 
grep NODE_65_length_416_cov_233.578171 cleaned_WB3_scaffold_filtered.fasta 
```
```
contamination_line=$(grep -n "NODE_22_length_2705_cov_188.989726" filtered_scaffolds/WB8II_scaffold_filtered.fasta | cut -d ":" -f 1)
next_line=$((contamination_line + 1)) 
sed -e "${contamination_line},${next_line}d" filtered_scaffolds/WB8II_scaffold_filtered.fasta > cleaned_WB8II_scaffold_filtered.fasta 
grep NODE_22_length_2705_cov_188.989726 cleaned_WB8II_scaffold_filtered.fasta 
```

#### Renaming them so I can submit them again on NCBI
```
rename 's/cleaned_//' cleaned_*
mv *.fasta filtered_scaffolds
```

# *3. Genome Mining - Environmental Bacteria*
Rapid annotation was performed using RAST web server: https://rast.nmpdr.org/

## Processing RAST output 
Check script 09.Process_RAST_out.ipynb
# *4. Bioremediation assessment* 
Performed by analyzing growth curves, more specifically their growth rates
Check script 10.Visualize_isolates_growth.ipynb
