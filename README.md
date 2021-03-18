# paper_cipro_EC_AB_PA_2021
Supplementary code for "Morbidostat-Directed Experimental Evolution of Ciprofloxacin Resistance in Gram-negative Bacterial Pathogens" paper.

1.  **Software list**

> All computations were made on CentOS v.7.6.

  - R v. 3.6.0 (<https://www.r-project.org/>)

  - Python 3.6.7 (<https://www.python.org/>)

  - FastQC v. 0.11.8
    (<https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>)
    <sup>1</sup>

  - Trimmomatic v. 0.36
    (<http://www.usadellab.org/cms/?page=trimmomatic>) <sup>2</sup>

  - BWA v. 0.7.13 (<http://bio-bwa.sourceforge.net/>) <sup>3</sup>

  - Picard v. 2.2.1 (<https://broadinstitute.github.io/picard/>)

  - LoFreq v. 2.1.3.1 (<https://csb5.github.io/lofreq/>) <sup>4</sup>

  - Samtools + bcftools v. 1.3 (<http://www.htslib.org/>) <sup>3,
    5</sup>

  - Genome analysis Toolkit (GATK) v. 3.5
    (<https://software.broadinstitute.org/gatk/>) <sup>6, 7</sup>

  - Qualimap v. 2.2 (<http://qualimap.bioinfo.cipf.es/>) <sup>8</sup>

  - snpEff v. 4.3 (<http://snpeff.sourceforge.net/>) <sup>9</sup>

  - breseq v. 0.33.2 (<http://barricklab.org/breseq>) <sup>10</sup>

  - CNOGpro (R package, installation through CRAN) <sup>11</sup>

  - Albacore v. 2.3.4 (<https://community.nanoporetech.com/>; Customer
    login required)

  - Prechop v. 0.2.4 (<https://github.com/rrwick/Porechop>)

  - MinionQC v. 1.3.0 (<https://github.com/roblanf/minion_qc>)
    <sup>12</sup>

  - SPAdes v. 3.13.0 (<http://cab.spbu.ru/software/spades/>)
    <sup>13</sup>

  - MUMmer v. 3.1 (<http://mummer.sourceforge.net/>) <sup>14</sup>

  - iJump (<https://github.com/sleyn/ijump>)

**Custom scripts (available on GitHub):**

  - timeline.py

  - CNOGpro.R

  - breseq\_parser\_confident.py

  - spread\_gd\_tsv.R

  - make\_repeats.pl

<!-- end list -->

2.  **Reads quality control**

We accessed read quality with FastQC software <sup>1</sup>. The Q-scores
in all samples had a very good distribution across all positions in
reads and we did not trim reads based on Q-scores. The adapters were
trimmed with Trimmomatic <sup>2</sup>:

java -jar trimmomatic-0.36.jar PE \\

\-phred33 \\

\[Input reads 1\] \\

\[Input reads 2\] \\

\[Trimmed output paired reads 1\] \\

\[Trimmed output single reads 1\] \\

\[Trimmed output paired reads 1\] \\

\[Trimmed output single reads 2\] \\

ILLUMINACLIP:adapters.fasta:2:30:10

MINLEN:65

The adapter sequences used for trimming were copied from the Illumina
Adapter Sequences document for the corresponding library preparation
kit. Additionally, we added overrepresented sequences to the adapters
fasta file as we found that it improves read quality. These
over­­represented sequences in most cases were adapters with barcodes.

3.  **Population sequencing data analysis**
    
    1.  **Align reads on reference**

Processed reads were aligned to a reference genome with BWA-MEM
<sup>15</sup>. First, we should create read group header that is
required for Genome Analysis Toolkit (GATK) <sup>6</sup>:

read1=\[trimmed paired reads 1\]

rgid=$(echo $read1)

rgid=${rgid/\_1\_trimmed.fq.gz/} \# Leave only ID

rg="@RG\\\\tID:$rgid\\\\tSM:S1\\\\tPL:illumina\\\\tLB:L001\\\\tPU:R1"

Then run BWA-MEM itself:

bwa mem \\

\-t \[Number of CPUs\] \\

\-M \\ \# Parameter used for Picard compatibility

\-R "$rg" \\ \# Read group header

\-Y \\ \# Write all split hits as soft clipped

\[path to the reference fasta\] \\

\[Trimmed paired reads 1\] \\

\[Trimmed paired reads 2\] \> Sample.sam

Sort SAM file, save it as BAM and create index with Picard tools
(<https://broadinstitute.github.io/picard/>):

java -jar \~/ngsbin/picard-tools-2.2.1/picard.jar SortSam \\

I=Sample.sam \\

O=Sample.bam \\

SO=coordinate \\

CREATE\_INDEX=true

2.  **Alignment recalibration**

Realign reads with LoFreq Vitebri algorithm <sup>4</sup>. As reads
positions could be changed, we repeat sorting on the output bam:

lofreq viterbi \\

\--ref \[Path to the reference fasta\] \\

\--keepflags \\

\-o Sample\_v.bam \\

Sample.bam

samtools sort -o Sample\_i.bam Sample\_v.bam

samtools index -b Sample\_i.bam

We performed the recommended Base Quality Score Recalibration (BQSR)
with GATK. However, as there is no database of know SNP for *A.
baumannii* and *E. coli* are available. The official advise from GATK is
to make a first round of calling first and choose the most trusted
mutations. However, as we have a high coverage samples and most trusted
mutations would be anyway with high frequency we decided to skip it and
make BQSR with the empty SNP database:

java -jar \~/ngsbin/GATK-3.5/GenomeAnalysisTK.jar \\

\-T BaseRecalibrator \\

\-R Path to the reference fasta\] \\

\-I Sample\_i.bam \\

\-knownSites known.vcf \\

\-o Sample\_re.bqsr.grp

java -jar \~/ngsbin/GATK-3.5/GenomeAnalysisTK.jar \\

\-T PrintReads \\

\-R \[Path to the reference fasta\] \\

\-I Sample\_i.bam \\

\-BQSR Sample\_re.bqsr.grp \\

\-o Sample\_r\_bqsr.bam

To assess alignment quality we used Qualimap <sup>8</sup>:

qualimap bamqc \\

\-bam Sample\_r\_bqsr.bam \\

\-outfile BamQC.pdf

3.  **SNP and indel variant calling**

The variant calling of SNP and indels was performed with LoFreq:

lofreq call-parallel \\

\--pp-threads \[Number of CPUs\] \\

\--call-indels \\

\-f \[Path to the reference fasta\] \\

\-o variants.lofreq.vcf \\

Sample\_r\_bqsr.bam

4.  **IS elements rearrangement**

<!-- end list -->

  - The search of mobile elements rearrangement was done with iJump
    (<https://github.com/sleyn/ijump>) – the software that we developed
    to search for IS elements rearrangements in a population sequencing
    data. In general algorithm of iJump is following:
    
    1.  Find soft-clipped reads near boundaries of mobile elements.
    
    2.  Extract unaligned part of reads.
    
    3.  BLAST unaligned parts against the reference.
    
    4.  Find best hits with \>90% identity to the reference with unique
        highest bitscore. If two or more hits share the same bit score -
        the junction considered as ambiguous and skipped. Usually this
        is happening because aligner maps reads to several copies of
        mobile element.
    
    5.  Assess frequency of insertion by a formula:

\[Frequency = \frac{\frac{R_{l} + \ R_{r}}{2}*(1 + \ \frac{B_{\min}}{A_{\text{rlen}}})}{D_{t}*(1 - \ \frac{\text{mmatch}}{A_{\text{rlen}}})}\]

> where:
> 
> *R<sub>l</sub>* – a number of reads that support junction to the
> target on the "left" side of mobile element.
> 
> *R<sub>r</sub>* – a number of reads that support junction to the
> target on the "right" side of mobile element.
> 
> *Dt* – average depth of coverage of the target region.
> 
> *B<sub>min</sub>* – a threshold for the minimum length for unaligned
> parts of the reads used in the BLAST step. In algorithm we used 10nt.
> 
> *A<sub>rlen</sub>* – an average read length;
> 
> *mmatch* – a minimum length of the read part that could be aligned to
> reference. The *mmatch* parameter is accessed from the data as a
> minimum of longest clipped part of the read (e.g. for read in the SAM
> file with CIGAR string 10S120M30S mmatch is 30).

Before using iJump we made a BLAST search of reference genomes against
ISFinder database <sup>16</sup>. The output HTML files were processed
with the isfinder\_parse.py script that comes with the iJump.

python3 isfinder\_parse.py \\

\-i \[ISFinder BLAST HTML page\]

The script filters BLAST output and creates a table of IS elements
boundaries that is used by iJump. The iJump is running with the
following command:

python3 isjump.py \\

\-a \[Alignment BAM of SAM file\] \\

\-r \[Reference FNA file with nucleotide FASTA\] \\

\-g \[Reference GFF file with annotations\]

\-i \[IS coordinates\]

> Detailed description of iJump usage and input data formats could be on
> the GitHub page (<https://github.com/sleyn/ijump>).

1.  **Variant annotation and mutational dynamics analysis**

> To annotate variant calling results we used snpEff <sup>9</sup> and
> custom scripts. All VCF files from LoFreq output were copied to one
> directory and merged into one VCF file using bcftools:
> 
> for file in \*.vcf
> 
> do
> 
> gz=${file/.vcf/.gz}
> 
> bgzip -c $file \> $gz
> 
> tabix -p vcf $gz
> 
> done
> 
> bcftools merge \\
> 
> \-m none \\
> 
> \-o merged.vcf \\
> 
> \*.gz
> 
> snpEff was applied on the merged file:
> 
> java -jar \~/ngsbin/snpEff/snpEff.jar \\
> 
> \-ud 100 \\
> 
> \-no-downstream \\
> 
> \-no-upstream \\
> 
> \-v \[database for the selected reference genome\] \\
> 
> merged.vcf \> ann.vcf
> 
> Results of the LoFreq variant calling were combined into one
> comparative table using custom timeline.py script.
> 
> python3 timeline.py \\
> 
> \-g \[Reference GFF file with annotations\] \\
> 
> \-a ann.vcf \\
> 
> \-r \[File with repeats (see below)\] \\
> 
> \-v \[Path to the directory with VCF files\] \\
> 
> \-o summary\_table.txt
> 
> Additionally, to manual curation of the variant list we mark two cases
> of suspicious variants:

1.  Variants in the potential repetitive regions where mapped reads
    could be attributed to the wrong copy of the region by error.

> To do this we followed instructions from MUMmer <sup>14</sup> manual
> to find repeats in the reference genome:
> 
> prefix=\[Prefix of the output files\]
> 
> nucmer \\
> 
> \-p $prefix \\
> 
> \-nosimplify \\
> 
> \--maxmatch \\
> 
> \[Reference FNA file\] \\
> 
> \[Reference FNA file\]
> 
> show-coords \\
> 
> \-clrT \\
> 
> \-I 90 \\
> 
> \-L 65 \\
> 
> $prefix".delta" \> \\
> 
> $prefix".aln"
> 
> make\_repeats.pl $prefix
> 
> This procedure will create file with potential repeated regions in the
> genome with identity 90% or more and length at least 65nt.

2.  Variants that have stable less than 100% frequency in the
    population. It is another sign that alignment software has troubles
    with this region.

> To test uniformity of frequencies of a variant across samples for each
> reactor we preform chi-squared test for uniform distribution. This
> test id done automatically with the timeline.py script.
> 
> Results of the iJump IS element rearrangement calling were combined to
> a comparative table using combine\_results.py script from the iJump
> package.
> 
> python3 combine\_results.py \\
> 
> \-d \[Path to the directory with iJump report files\]

4.  **Clonal sequencing data analysis**
    
    1.  **Variant calling**

We used the breseq software <sup>10</sup> to call variants in clonal
data:

breseq \\

\-j \[number of CPUs\] \\

\-n \[Clone ID\] \\

\-r \[Reference GBK file\] \\

\--require-match-length 40 \\

\[Trimmed paired reads 1\] \\

\[Trimmed paired reads 2\]

The Copy Number Variation analysis was done using CNOGpro package for
the R programming language. The CNOGpro package uses GBK files but it is
limited for using files with one contig only. To use CNOGpro it is
necessary to split multiple contigs GBK file to the separate files for
each contig. For each contig the hitsfile containing coordiantes of
mapped reads should be generated (for details see CNOGpro manual):

\# in breseq data directory for each contig:

export CONTIG=\[Contig name\]

samtools view reference.bam | \\

perl -lane \\

'print "$F\[2\]\\t$F\[3\]" \\

if $F\[2\] =\~ m/$ENV{CONTIG}/' \\

\> out.hits

> CNOGpro was wrapped by us into an R script. Script filters out IS
> elements and outputs only CDS:
> 
> Rscript --vanilla CNOGpro.R \\
> 
> \[Reference GBK file with the selected contig\] \\
> 
> \[Contig name\] \\
> 
> \[iJump IS elements coordinates table to filter out them\] \\
> 
> \[Clone ID used for output file name\]

2.  **Mutational data comparison**

To collect breseq results in the table format we used custom script that
makes TSV table from the annotated GD files (\[breseq work
directory\]/output/evidence/annotated.gd):

python3 breseq\_parser\_confident.py \\

\-g \[directory with annotated gd files\] \\

\-o \[output file\]

To make a comparison table for variants found in the clonal sequencing
data we used gdtools with a custom script for reformatting the gdtools
output:

gdtools ANNOTATE \\

\-r \[Reference GBK file\] \\

\-o \[Output file\] \\

\-f TSV \\ \# TSV format

\[Directory with GD files\]/\*.gd

Rscript --vanilla spread\_gd\_tsv.R \\

\[gdtools output file\] \\

\[Output table\]

> To make a comparison table for copy number variations found with
> CNOGpro we used a custom script:
> 
> Rscript --vanilla Reshape\_cnv.R \\
> 
> \[Directory with the CNOGpro outputs\] \\
> 
> \[Output table\] \\
> 
> \[Minimum difference from 1 to show gene in the output. We used 0.3\]

5.  **Novel assembly of *Acenitobacter baumannii* ATCC 17978 genome**

> Initial variant calling against public reference revealed high number
> of variants in unevolved culture of *Acinetobacter baumannii* ATCC
> 19978 (number of variants with frequency \>85%: 16 variants against
> GCA\_001593425.2, 98.46% mapped reads; 87 variants against
> GCA\_000015425.1, 98.97% mapped reads). To improve variant calling we
> decided to make a hybrid assembly with Oxford Nanopore long reads.
> 
> One clone of unevolved *Acinetobacter baumannii* ATCC 19978 was
> sequences with ONT MinION sequencer using Nanopore Rapid Barcoding
> gDNA Sequencing kit (SQK-RBK004). Raw sequencing data QC was performed
> with MinionQC.
> 
> Reads were demultiplexed with Albacore basecaller (ONT community
> site):
> 
> read\_fast5\_basecaller.py \\
> 
> \-i \[Path to the FAST5 files directory\] \\
> 
> \-t \[Number of CPUs\] \\
> 
> \-s \[Path to the directory for output FASTQ file\] \\
> 
> \-f FLO-MIN106 \\
> 
> \-k SQK-RBK004 \\
> 
> \-r \\
> 
> \-o fastq \\
> 
> \-q 0
> 
> As in the library all barcoded samples were from the same organism to
> increase coverage, we decided to make second round of demultiplexing
> on unclassified reads with Porechop
> (<https://github.com/rrwick/Porechop>):
> 
> porechop \\
> 
> \-i \[Path to the unclassified/ folder\] \\
> 
> \-b Porechop\_barcodes/ \\
> 
> \--threads \[Number of CPUs\] \\
> 
> \--untrimmed \\
> 
> \--discard\_middle
> 
> Adapters from the reads demultiplexed by Albacore were trimmed by
> Porechop:
> 
> porechop \\
> 
> \-i \[FASTQ file\] \\
> 
> \-o \[Output trimmed FASTQ file\] \\
> 
> \--discard\_middle \\
> 
> \--threads \[Number of CPUs\]
> 
> Both reads demultiplexed by Albacore and Porechop were combined by
> simple cat command. The hybrid assembly of *Acinetobacter baumannii*
> ATCC 19978 genome was performed using SPAdes assembler <sup>13</sup>:
> 
> spades.py \\
> 
> \-k 21,33,55,77 \\
> 
> \--careful \\
> 
> \--nanopore \[Nanopore fasta files\] \\

\[Trimmed Illumina paired reads 1\] \\

\[Trimmed Illumina paired reads 2\]

> Assembled genome was annotated by RASTtk pipeline on the corresponding
> web server (<https://rast.nmpdr.org/>) <sup>17</sup>.

**References**

> 1\. Andrews S. FastQC. A quality control tool for high throughput
> sequence data.). Babraham Institute (2010).2. Bolger AM, Lohse M,
> Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data.
> *Bioinformatics* **30**, 2114-2120 (2014).3. Li H*, et al.* The
> Sequence Alignment/Map format and SAMtools. *Bioinformatics* **25**,
> 2078-2079 (2009).4. Wilm A*, et al.* LoFreq: a sequence-quality aware,
> ultra-sensitive variant caller for uncovering cell-population
> heterogeneity from high-throughput sequencing datasets. *Nucleic Acids
> Res* **40**, 11189-11201 (2012).5. Li H. A statistical framework for
> SNP calling, mutation discovery, association mapping and population
> genetical parameter estimation from sequencing data. *Bioinformatics*
> **27**, 2987-2993 (2011).6. DePristo MA*, et al.* A framework for
> variation discovery and genotyping using next-generation DNA
> sequencing data. *Nat Genet* **43**, 491-498 (2011).7. Van der Auwera
> GA*, et al.* From FastQ data to high confidence variant calls: the
> Genome Analysis Toolkit best practices pipeline. *Curr Protoc
> Bioinformatics* **43**, 11 10 11-11 10 33 (2013).8. Okonechnikov K,
> Conesa A, Garcia-Alcalde F. Qualimap 2: advanced multi-sample quality
> control for high-throughput sequencing data. *Bioinformatics* **32**,
> 292-294 (2016).9. Cingolani P*, et al.* A program for annotating and
> predicting the effects of single nucleotide polymorphisms, SnpEff:
> SNPs in the genome of Drosophila melanogaster strain w1118; iso-2;
> iso-3. *Fly (Austin)* **6**, 80-92 (2012).10. Deatherage DE, Barrick
> JE. Identification of mutations in laboratory-evolved microbes from
> next-generation sequencing data using breseq. *Methods Mol Biol*
> **1151**, 165-188 (2014).11. Brynildsrud O, Snipen LG, Bohlin J.
> CNOGpro: detection and quantification of CNVs in prokaryotic
> whole-genome sequencing data. *Bioinformatics* **31**, 1708-1715
> (2015).12. Lanfear R, Schalamun M, Kainer D, Wang W, Schwessinger B.
> MinIONQC: fast and simple quality control for MinION sequencing data.
> *Bioinformatics* **35**, 523-525 (2019).13. Bankevich A*, et al.*
> SPAdes: a new genome assembly algorithm and its applications to
> single-cell sequencing. *J Comput Biol* **19**, 455-477 (2012).14.
> Kurtz S*, et al.* Versatile and open software for comparing large
> genomes. *Genome Biol* **5**, R12 (2004).15. Li H, Durbin R. Fast and
> accurate short read alignment with Burrows-Wheeler transform.
> *Bioinformatics* **25**, 1754-1760 (2009).16. Siguier P, Perochon J,
> Lestrade L, Mahillon J, Chandler M. ISfinder: the reference centre for
> bacterial insertion sequences. *Nucleic Acids Res* **34**, D32-36
> (2006).17. Brettin T*, et al.* RASTtk: a modular and extensible
> implementation of the RAST algorithm for building custom annotation
> pipelines and annotating batches of genomes. *Sci Rep* **5**, 8365
> (2015).
