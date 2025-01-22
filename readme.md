# Day 1 Protocol

Learned some Linux commands, most of them were familiar for me.

[Linux cheat sheet](https://files.fosswire.com/2007/08/fwunixref.pdf)

New ones:
- chmod (including numerical parameters 4 - read, 2 - write, 1 - execute). First group is the root user, then the group user, and lastly the users the file is shared with (world). \
  ``` chmod 777 filename.extension #all permissions for everyone ```

How to access the cluster of CAU:

```bash
ssh -X sunam___@caucluster.rz.uni-kiel.de
```

Then enter the password.

Learned about environments.

****************
# Day 2 Protocol
## Lecture highlights
- Workflow for metagenomic analysis
- De Brujn graph for contig generation from short reads
- Some history of genomics

## Quality control, trimming, assembly

1. Connect the caucluster.
2. Go to working directory.
3. Load the required modules.
4. Activate the modules.


This is the bash script template for running jobs in caucluster, which includes some of the aforementioned steps:

``` bash 
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --partition=base
#SBATCH --time=5:00:00
#SBATCH --reservation=biol217

#load necessary modules
module load gcc12-env/12.1.0
module load micromamba/1.4.2
eval "$(micromamba shell hook --shell=bash)"
export MAMBA_ROOT_PREFIX=$WORK/.micromamba

cd $WORK

micromamba activate .micromamba/envs/00_anvio/
# WRITE YOUR COMMANDS

# ##----------------- End -------------
module purge
jobinfo
```

**************

### Quality control:
``` bash
for file in /work_beegfs/sunam230/metagenomics/0_raw_reads/*.gz; do
    fastqc "${file}" -o /work_beegfs/sunam230/metagenomics/fastqc_output/
done
```
Result (from .html file):
![Image](./resourses/fastqcquality.png)

Output: fastqc.zip files with information about quality only, no raw reads there.

### Cleaning the data:

```
fastp -i /work_beegfs/sunam230/metagenomics/0_raw_reads/BGR_130305_mapped_R1.fastq.gz -I /work_beegfs/sunam230/metagenomics/0_raw_reads/BGR_130305_mapped_R2.fastq.gz -R fastp_report1 -o /work_beegfs/sunam230/metagenomics/fastp_output/130305_R1_clean.fastq.gz -O /work_beegfs/sunam230/metagenomics/fastp_output/130305_R2_clean.fastq.gz -t 6 -q 20
fastp -i /work_beegfs/sunam230/metagenomics/0_raw_reads/BGR_130527_mapped_R1.fastq.gz -I /work_beegfs/sunam230/metagenomics/0_raw_reads/BGR_130527_mapped_R2.fastq.gz -R fastp_report2 -o /work_beegfs/sunam230/metagenomics/fastp_output/130587_R1_clean.fastq.gz -O /work_beegfs/sunam230/metagenomics/fastp_output/130587_R2_clean.fastq.gz -t 6 -q 20
fastp -i /work_beegfs/sunam230/metagenomics/0_raw_reads/BGR_130708_mapped_R1.fastq.gz -I /work_beegfs/sunam230/metagenomics/0_raw_reads/BGR_130708_mapped_R2.fastq.gz -R fastp_report3 -o /work_beegfs/sunam230/metagenomics/fastp_output/130708_R1_clean.fastq.gz -O /work_beegfs/sunam230/metagenomics/fastp_output/130708_R2_clean.fastq.gz -t 6 -q 20
```
Output: fastq.gz files with cleaned (trimmed) reads.

### Assembly:
```
megahit -1 /work_beegfs/sunam230/metagenomics/fastp_output/130305_R1_clean.fastq.gz -1 /work_beegfs/sunam230/metagenomics/fastp_output/130587_R1_clean.fastq.gz -1 /work_beegfs/sunam230/metagenomics/fastp_output/130708_R1_clean.fastq.gz -2 /work_beegfs/sunam230/metagenomics/fastp_output/130305_R2_clean.fastq.gz -2 /work_beegfs/sunam230/metagenomics/fastp_output/130587_R2_clean.fastq.gz -2 /work_beegfs/sunam230/metagenomics/fastp_output/130708_R2_clean.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 -o /work_beegfs/sunam230/metagenomics/assembly_output/ -t 12  
```

Kept running overnight.
**************


### Errors&mistakes:
- Typos in path
  > no backslash in the beginning of the path
- Wrong directory of the .sh script

Used this [MD cheat sheet](https://enterprise.github.com/downloads/en/markdown-cheatsheet.pdf)


# Day 3 Protocol
## Assembly visualisation


```bash

```


```bash

```

![Image](./resourses/assembly_graph.png)

We used Bandage software to open .fastg files (we obtained them after converting the final assembly .fa files).

### Questions
We see a lot of small contigs and a small number of big ones (more than 1000 bp).

## Quality assessment of assembly

```bash

```

    What is your N50 value? Why is this value relevant?
    How many contigs are assembled?
    What is the total length of the contigs?
### N50
![Image]()

### Number of contigs
![Image]()

### Total length
![Image]()


## Genome Binning
### Preparations and mapping
format your fasta sequence IDs
```bash
anvi-script-reformat-fasta /work_beegfs/sunam230/metagenomics/assembly_output/final.contigs.fa -o /work_beegfs/sunam230/metagenomics/assembly_output/contigs.anvio.fa --min-len 1000 --simplify-names --report-file name_conversion.txt
```

Use the following command to index your mapping reference fasta file.
```bash
bowtie2-build /work_beegfs/sunam230/metagenomics/assembly_output/contigs.anvio.fa /work_beegfs/sunam230/metagenomics/assembly_output/contigs.anvio.fa.index
```
Now use bowtie2 for the actual mapping.
```bash
cd /work_beegfs/sunam230/metagenomics/fastp_output/

for i in *_R1_clean.fastq.gz; do
  base="${i%_R1_clean.fastq.gz}"; bowtie2 --very-fast -x /work_beegfs/sunam230/metagenomics/assembly_output/contigs.anvio.fa.index -1 $i -2 "$base"_R2_clean.fastq.gz -S /work_beegfs/sunam230/metagenomics/mapping/"${base}".sam 
done
```

sequence mapping file (SAM) with the .sam extension and which we convert to binary alignment and map (BAM). Then contigs.db generation and HMM search (adds hmm.hits into contigs.db).
```bash
for i in /work_beegfs/sunam230/metagenomics/mapping/*.sam; do samtools view -bS $i > "$i".bam; done

anvi-gen-contigs-database -f /work_beegfs/sunam230/metagenomics/assembly_output/contigs.anvio.fa -o contigs.db -n 'biol217'

anvi-run-hmms -c contigs.db --num-threads 12
```

### Visualisation of contigs.db in ANVI´O interactive
In terminal:
```bash
srun --pty --x11 --partition=interactive --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --mem=10G --time=01:00:00 /bin/bash
```
We started a server. 
`note down which node this command logged you on`

```
module load gcc12-env/12.1.0
module load micromamba
micromamba activate $WORK/.micromamba/envs/00_anvio/
```

`Run the command to display what you want`

Then open a new terminal
```
ssh -L 8080:localhost:8080 sunam230@caucluster.rz.uni-kiel.de
```
```
ssh -L 8080:localhost:8080 n# #replace the node number
```

http://127.0.0.1:8080/


Then close the connection with `ctrl c`
```
exit
```
# if the host is busy, try 8080 instead of 8060

http://127.0.0.1:8060/


****************

### Profiling with ANVI´O

>Sort and Index bam files
```bash
for i in *.bam; do anvi-init-bam $i -o "$i".sorted.bam; done
```

>Creating an Anvi’o Profile
```bash
cd /work_beegfs/sunam230/metagenomics/mapping/
for i in *.sorted.bam; do
    base=$(basename "$i" .sam.bam.sorted.bam)
    anvi-profile -i "$i" -c contigs.db -o ../profiling/${base}
done
```

>merge the profiles coming from your different samples into one profile:
```bash
anvi-merge /work_beegfs/sunam230/metagenomics/profiling/130305/PROFILE.db /work_beegfs/sunam230/metagenomics/profiling/130587/PROFILE.db /work_beegfs/sunam230/metagenomics/profiling/130708/PROFILE.db -o /work_beegfs/sunam230/metagenomics/profiling/merged_profiles/ -c /work_beegfs/sunam230/metagenomics/mapping/contigs.db --enforce-hierarchical-clustering
```

## Binning with Metabat2

- 3 archaea bins
