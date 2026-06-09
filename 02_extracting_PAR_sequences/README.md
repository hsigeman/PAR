# 2. Extracting sex-specific PAR gene sequences

## 2.1 Alignment and variant calling

Samples listed in Table S1 were aligned to the reference genome constructed from the male sample of each species, using a snakemake pipeline. 

```bash
from snakemake.utils import min_version
import os

min_version("4.4.0")

# Set directory paths
dir_path = os.getcwd()

###############################################
################## PATHS ######################
###############################################

ID = config["samples"]
SPECIES = config["species"]
FEMALE = config["female"]
MALE = config["male"]
FQ_DIR = config["fastq"]
REF_SPECIES = config["ref_species"]
REF_DIR = config["ref_dir"]
REF_NAME = config["ref_name"]
PREFIX = SPECIES + "_ref_" + REF_SPECIES

REF_PATH = REF_DIR + REF_NAME
REF_FASTA = REF_DIR + REF_NAME + ".fasta"
MAP_DIR = "intermediate/bwa/" + PREFIX + "/" 
GENCOV_DIR = "intermediate/bedtools" + PREFIX + "/"
VCF_DIR = "intermediate/freebayes" + PREFIX + "/"
RESULTDIR = "results/" + PREFIX + "/"

###############################################
################## RULES ######################
###############################################

rule all: 
    input: 
        REF_FASTA + ".bwt",
        REF_FASTA + ".fai",
        expand(MAP_DIR + "{S}" + ".sorted.status", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.status", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.bam.bai", S = ID),
        VCF_DIR + SPECIES + ".vcf.status",
        expand(MAP_DIR + "{S}" + ".sorted.flagstat", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.flagstat", S = ID),
        VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf",
        VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf.gz"

##########################################################  
##################### INDEX GENOME #######################      
##########################################################

rule index_fasta_bwa:
    input: 
        ref = REF_FASTA
    output:
        ref_bwt = REF_FASTA + ".bwt"
    priority: 80
    message: "Indexing {input} with BWA index."
    threads: 2
    shell:
        """
        bwa index {input}
        """

rule index_fasta_samtools: 
    input: 
        ref = REF_FASTA
    output: 
        ref_fai = REF_FASTA + ".fai"
    priority: 70
    threads: 2
    shell: 
        """
        samtools faidx {input}
        """

##########################################################  
######################## MAPPING #########################       
##########################################################   

rule map: 
    input: 
        R1 = FQ_DIR + "{S}_forward_paired.fq.gz",
        R2 = FQ_DIR + "{S}_reverse_paired.fq.gz",
        ref = REF_FASTA, 
        ref_bwt = REF_FASTA + ".bwt"
    output: 
        temp(MAP_DIR + "{S}" + ".bam")
    message: "Mapping reads to ref"
    threads: 20
    params:
        rg = "\"@RG\\tID:{S}\\tSM:{S}\""
    shell:
        """ 
        bwa mem -t {threads} -M -R {params.rg} {input.ref} {input.R1} {input.R2} | samtools view -Sb - > {output}
        """ 

rule sort_bam:
    input:
        MAP_DIR + "{S}" + ".bam"
    output:
        out = temp(MAP_DIR + "{S}" + ".sorted.bam"),
        log = MAP_DIR + "{S}" + ".sorted.status"
    threads: 15
    params:
        tmpdir = MAP_DIR + "{S}" + "_temp_sort/"
    shell:
        """
        mkdir {params.tmpdir}
        samtools sort -@ {threads} {input} -T {params.tmpdir} > {output.out}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
        """

rule remove_duplicates: 
    input: 
        MAP_DIR + "{S}" + ".sorted.bam"
    output: 
        out = MAP_DIR + "{S}" + ".sorted.nodup.bam",
        log = MAP_DIR + "{S}" + ".sorted.nodup.status"
    params:
        tmpdir = MAP_DIR + "{S}" + "_temp_dupl/"
    shell: 
        """
        mkdir {params.tmpdir}
        picard MarkDuplicates MAX_FILE_HANDLES=500 REMOVE_DUPLICATES=true I={input} O={output.out} M={input}_duplicatedata.txt TMP_DIR={params.tmpdir}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
        """

rule index_bam: 
    input: 
        MAP_DIR + "{S}" + ".sorted.nodup.bam"
    output: 
        MAP_DIR + "{S}" + ".sorted.nodup.bam.bai"
    threads: 1
    shell: 
        """
        samtools index {input}
        """

##########################################################  
#################### GENOME COVERAGE #####################       
########################################################## 

rule gencov_prepare_fasta:
    input: 
        ref_fai = REF_FASTA + ".fai"
    output: 
        GENCOV_DIR + "genome_5kb_windows.out"
    threads: 1
    shell: 
        """
        bedtools makewindows -g {input} -w 5000 -s 5000 > {output}
        """

rule gencov_bedtoolsall:
    input: 
        bam_f = MAP_DIR + FEMALE + ".sorted.nodup.bam",
        bai_f = MAP_DIR + FEMALE + ".sorted.nodup.bam.bai",
        bam_m = MAP_DIR + MALE + ".sorted.nodup.bam",
        bai_m = MAP_DIR + MALE + ".sorted.nodup.bam.bai",
        bed = GENCOV_DIR + "genome_5kb_windows.out"
    output: 
        GENCOV_DIR + "gencov.nodup.nm.all.out"
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam_f} {input.bam_m} -bed {input.bed} > {output}
        """

##########################################################  
#################### VARIANT CALLING #####################      
########################################################## 

rule freebayes_prep:
    input: 
        ref_fai = REF_FASTA + ".fai"
    output: 
        VCF_DIR + SPECIES + ".100kbp.regions"
    threads: 4
    shell: 
        """
        fasta_generate_regions.py {input} 100000 > {output}
        """

rule freebayes_parallel:
    input: 
        ref = REF_FASTA,
        regions = VCF_DIR + SPECIES + ".100kbp.regions",
        f = MAP_DIR + FEMALE + ".sorted.nodup.bam",
        m = MAP_DIR + MALE + ".sorted.nodup.bam"
    output: 
        vcf = VCF_DIR + SPECIES + ".vcf",
        log = VCF_DIR + SPECIES + ".vcf.status"
    threads: 18
    params:
        tmpdir = VCF_DIR + "temp/"
    shell: 
        """
        mkdir {params.tmpdir}
        export TMPDIR={params.tmpdir}
        freebayes-parallel {input.regions} {threads} -f {input.ref} {input.f} {input.m} > {output.vcf}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
        """

rule vcftools_singletons:
    input: 
        VCF_DIR + SPECIES + ".vcf"
    output: 
        VCF_DIR + SPECIES + ".singletons.bed"
    threads: 1
    shell: 
        """
        vcftools --vcf {input} --singletons --remove-filtered-geno-all --minQ 20 --minDP 3 --stdout | awk -v OFS="\\t" '{{print $1,$2,$2+1,$3,$4,$5}}' > {output}
        """
```

## 2.2 Extracting gene sequences

```bash
# Set directories and file paths
krakendir="kraken/bTaeGut1.pat.W.v2/"
genes="ZF_PAR.genes.new.list"
home=$PWD

# Step 1: Create BED file for PAR genes ordered per exon
ls $krakendir/ | while read sp; do
  cat $genes | while read gene; do
    grep $gene $krakendir/${sp}/mapped.gtf | \
    awk '$3=="exon" {print}' | \
    cut -f 1,4,5,9 | \
    sed 's/;/\t/g' | \
    awk '{print $1,$2,$3,$5,$7,$9}' | \
    tr -d "\"" | \
    sed 's/ /\t/g' | \
    sort -k5,5 -k6,6g
  done > $krakendir/${sp}/${sp}.PAR.genes.bed
done

# Step 2: Extract lines from the GTF file in the same order
ls $krakendir/ | while read sp; do
  cat $krakendir/${sp}/${sp}.PAR.genes.bed | while read scaff start end gene trans exonnr; do
    cat $krakendir/${sp}/mapped.gtf | \
    awk '{ if ($1=="'"$scaff"'" && $4=="'"$start"'" && $5=="'"$end"'" && $3=="exon") print $0}' | \
    grep $trans
  done > $krakendir/${sp}/${sp}.PAR.genes.gtf
done

# Step 3: Convert GTF to BED12 format and sort
ls $krakendir/ | while read sp; do
  gtf2bed --do-not-sort < $krakendir/${sp}/${sp}.PAR.genes.gtf > $krakendir/${sp}/${sp}.PAR.genes.sorted.bed
done

# Step 4: Generate FASTA sequences for PAR genes
cat samples_ref_genome.list | while read sp ref; do
  cd $krakendir/${sp}
  
  # Split BED file into separate list files by gene
  awk '{print >> $4 ".list"; close($4)}' ${sp}.PAR.genes.sorted.bed
  
  # Extract FASTA sequences using bedtools
  ls | grep list | while read gene; do
    echo ">$gene"
    bedtools getfasta -fi $home/../data/internal_raw/genome/${ref}.fasta -bed $gene -s -name | \
    grep -v ">" | tr -d "\n"
    echo
  done | sed 's/.list//' > ${sp}.PAR.fasta
  cd $home
done

# Step 5: Consensus sequence for male and female using VCF
cat samples_sex_sameline_ref.tsv | while read female male sp ref; do
  # Filter and process VCF files
  vcftools --vcf intermediate/freebayes_17nov2019_parallel/${sp}_ref_${sp}/${sp}.vcf \
    --non-ref-ac 1 --min-alleles 2 --max-alleles 2 --remove-filtered-all \
    --recode --stdout --minQ 20 --minDP 5 --bed $krakendir/${sp}/${sp}.PAR.genes.bed > \
    $krakendir/${sp}/${sp}.filt.vcf

  # Compress and index VCF
  bgzip -c $krakendir/${sp}/${sp}.filt.vcf > $krakendir/${sp}/${sp}.filt.vcf.gz
  tabix -p vcf $krakendir/${sp}/${sp}.filt.vcf.gz
done

# Step 6: Generate consensus sequences for male and female
cat samples_sex_sameline_ref.tsv | while read female male sp ref; do
  # Create female and male consensus FASTA sequences
  cat $home/../data/internal_raw/genome/${ref}.fasta | \
  bcftools consensus $krakendir/${sp}/${sp}.filt.vcf.gz -s $female > \
  $krakendir/${sp}/${sp}.female.fasta

  cat $home/../data/internal_raw/genome/${ref}.fasta | \
  bcftools consensus $krakendir/${sp}/${sp}.filt.vcf.gz -s $male > \
  $krakendir/${sp}/${sp}.male.fasta
done

# Step 7: Separate exons into separate sequences for male and female
cat samples_ref_genome.list | while read sp ref; do
  # Create exon-specific BED file for PAR genes
  cat $krakendir/${sp}/${sp}.PAR.genes.sorted.bed | \
  awk '{print $1,$2,$3,$4"_"$15,$5,$6}' | \
  sed 's/ /\t/g' | tr -d ";\"" > $krakendir/${sp}/${sp}.PAR.exonSeparate.genes.sorted.bed
 
  cd $krakendir/${sp}
  
  # Split the exons into list files
  awk '{print >> $4 ".list"; close($4)}' ${sp}.PAR.exonSeparate.genes.sorted.bed

  # Extract FASTA for female exons
  bedtools getfasta -fi ${sp}.female.fasta -bed ${sp}.PAR.exonSeparate.genes.sorted.bed -s -name | \
  sed 's/(+)//' | sed 's/(-)//' | sed 's/::/\t/' | cut -f 1 > ${sp}.PAR.exonSeparate.female.fasta

  # Extract FASTA for male exons
  bedtools getfasta -fi ${sp}.male.fasta -bed ${sp}.PAR.exonSeparate.genes.sorted.bed -s -name | \
  sed 's/(+)//' | sed 's/(-)//' | sed 's/::/\t/' | cut -f 1 > ${sp}.PAR.exonSeparate.male.fasta

  # Clean up list files
  cd $home
done

# Step 8: Reverse complement sequences in negative orientation
module load Fastx/0.0.14
oneline_fasta() {
  awk '/^>/ {printf("\n%s\n",$0);next;} { printf("%s",$0);} END {printf("\n");}' "$1"
}

ls $krakendir | while read sp ; do cat $krakendir/${sp}/${sp}.PAR.genes.sorted.bed | sed 's/ /\t/g' | awk '{print $1,$6,$4"_"$15}' | tr -d ";\"" | sed 's/ /\t/g' | awk '$2=="-" {print}' | while read scaff orientation exon ; do oneline_fasta $krakendir/${sp}/${sp}.exonSeparate.male.female.fasta | grep $exon$ -A 1 | grep -v "^--" | awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' | fastx_reverse_complement ; done ; done > allSp.PAR.exonSeparate.female.male.outgroup.fasta 

ls $krakendir | while read sp ; do cat $krakendir/${sp}/${sp}.PAR.genes.sorted.bed | sed 's/ /\t/g' | awk '{print $1,$6,$4"_"$15}' | tr -d ";\"" | sed 's/ /\t/g' | awk '$2=="+" {print}' | while read scaff orientation exon ; do oneline_fasta $krakendir/${sp}/${sp}.exonSeparate.male.female.fasta | grep $exon$ -A 1 | grep -v "^--" | awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' ; done ; done >> allSp.PAR.exonSeparate.female.male.outgroup.fasta 

# Step 9: Make multi-species fasta files for each gene

cat allSp.PAR.exonSeparate.female.male.outgroup.fasta | grep ">" | sed 's/male_/\t/' | cut -f 2 | sort | uniq > allPARexons.outgroup.list

mkdir exonSeparate_outgroup
cat allSp.PAR.exonSeparate.female.male.outgroup.fasta | grep ">" | sed 's/male_/\t/' | cut -f 2 | sort | uniq | while read gene ; do oneline_fasta allSp.PAR.exonSeparate.female.male.outgroup.fasta | grep ${gene}$ -A 1 | grep -v "^--" | sed 's/SylAtr_1EV02922/SylAtr/' > exonSeparate_outgroup/${gene}.fasta ; done
```

## 2.3 Calculating sequencing depth

```bash
krakendir="kraken/bTaeGut1.pat.W.v2/"

cat scratch/PAR/local_species_path.tsv | while read path female male sp ; do echo "java -jar ~/bin/jvarkit/dist/bamstats04.jar $path/$female.sorted.nodup.bam $path/$male.sorted.nodup.bam --bed $krakendir/$sp/${sp}.Z.genes.bed > $krakendir/$sp/${sp}.Z.genes.bamstat04.out" ; done 

ls $krakendir | while read sp ; do cat $krakendir/$sp/${sp}.Z.genes.bamstat04.out | awk '{print $0 "\t" "'"$sp"'"}' ; done > allSp.Z.genes.bamstat04.out

# And lastly summarize
bedtools intersect -a allSp.Z.genes.bamstat04.out -b allSp.Z.genes.bed -wa -wb | cut -f 5,6,7,8,9,12,16,17 > allSp.Z.genes.bamstat04.geneInfo.out
```

## 2.4 Private allele (singleton) analysis

```bash

ls $krakendir | while read sp ; do vcftools --vcf $krakendir/${sp}/${sp}.filt.vcf --singletons --out $krakendir/${sp}/${sp}.filt.singletons ; done

ls $krakendir | while read sp ; do cat $krakendir/${sp}/${sp}.filt.singletons.singletons | grep -v CHROM | awk '{print $1,$2-1,$2,$3,$4,$5}' | sed 's/ /\t/g' | bedtools intersect -a $krakendir/$sp/${sp}.Z.genes.bed -b stdin -wa -wb > $krakendir/$sp/${sp}.Z.genes.singletons.krakenSNPs.bed ; done

ls $krakendir | while read sp ; do cat $krakendir/$sp/${sp}.Z.genes.singletons.krakenSNPs.bed | awk '{print $0 "\t" "'"$sp"'"}' ; done > results/PAR/allSp.Z.genes.singletons.krakenSNPs.bed

cat data/meta/samples_sex.tsv | while read sample sp sex ; do cat results/PAR/allSp.Z.genes.singletons.krakenSNPs.bed | grep $sample | awk '{print $0 "\t" "'"$sex"'"}' ; done | awk '$10=="S" {print}' |  cut -f 4,12,13,14 | sort | uniq -c | grep -v GalMod | awk '{print $1,$2,$3,$4,$5}' | sed 's/ /\t/g' > results/PAR/allSp.Z.genes.singletons.S.sum.krakenSNPs.out

```
