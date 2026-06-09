# 1. Identification of PAR genes

## 1.1 PAR genes in zebra finch

```bash
###############################################################################
#                         DOWNLOAD EXTERNAL GENOME DATA                       #
###############################################################################

# Latest zebra finch genome and annotations from NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gtf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_cds_from_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/822/105/GCF_008822105.2_bTaeGut2.pat.W.v2/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gff.gz

###############################################################################
#                         EXTRACT CHROMOSOME Z SEQUENCE                       #
###############################################################################

gunzip *.gz

samtools faidx GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna
samtools faidx GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna NC_045027.1 > GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.chrZ.fasta

###############################################################################
#                        CLEAN AND PROCESS GTF FILE                           #
###############################################################################

zebra_gtf="GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.gtf"
gtf_clean="GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.small.clean.gtf"

cat "$zebra_gtf" |
  sed 's/exon_number/\texon_number/' |
  sed 's/protein_id/\tprotein_id/' |
  sed 's/;/\t/2' |
  cut -f1-9,11,12 |
  sed -e 's/\t/; /9' -e 's/\t//9' > "$gtf_clean"

###############################################################################
#                     FILTER LONGEST TRANSCRIPT PER GENE                      #
###############################################################################

module load bioinfo-tools CGAT/0.3.3
source /sw/apps/bioinfo/CGAT/0.3.3/rackham/conda-install/etc/profile.d/conda.sh
conda activate base; conda activate cgat-s

gtf_longest="GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.small.clean.longestTranscript.gtf"

grep -v "unknown" "$gtf_clean" |
  cgat gtf2gtf --method=filter --filter-method=longest-transcript > "$gtf_longest"

###############################################################################
#              GENERATE EXON FASTA FOR LONGEST TRANSCRIPTS                   #
###############################################################################

gffread -g GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna -w longestTranscripts.exons.fa -x longestTranscripts.exons.separate.fa "$gtf_longest"

###############################################################################
#              EXTRACT PAR AND CHROMOSOME Z GENES                             #
###############################################################################

chrZ="NC_045027.1"

# Genes in pseudoautosomal region (PAR) on Z
grep "$chrZ" "$gtf_longest" |
  awk '$3 == "exon" && $5 < 550000' |
  cut -f9 |
  awk '{print $2,$4}' |
  tr -d ';"' |
  sort -u |
  sed 's/ /\t/' |
  cut -f2 > ZF_PAR.genes.new.list

# Gene and transcript IDs in PAR
grep "$chrZ" "$gtf_longest" |
  awk '$3 == "exon" && $5 < 550000' |
  cut -f9 |
  awk '{print $2,$4}' |
  tr -d ';"' |
  sort -u |
  sed 's/ /\t/' > PAR_genes_geneID_transID.list

# All Z chromosome genes
grep "$chrZ" "$gtf_longest" |
  awk '$3 == "exon"' |
  cut -f9 |
  awk '{print $2,$4}' |
  tr -d ';"' |
  sort -u |
  sed 's/ /\t/' |
  cut -f2 > ZF_Z_chr.genes.list

# All Z genes with gene and transcript ID
grep "$chrZ" "$gtf_longest" |
  awk '$3 == "exon"' |
  cut -f9 |
  awk '{print $2,$4}' |
  tr -d ';"' |
  sort -u |
  sed 's/ /\t/' > ZF_Z_chr.genes_geneID_transID.list

# W-linked PAR genes
grep "NW_022611471.1" "$gtf_longest" |
  cut -f9 |
  awk '{print $2,$4}' |
  tr -d ';"' |
  sort -u |
  sed 's/ /\t/' > W_PAR_genes_geneID_transID.list

# Make table with gene names 
cat <<EOF > W_PAR_genes_geneID_transID_withGeneName.list
NEDD4L    LOC100217943    XM_030259050.2
ZNF532    LOC100218930    XM_030258323.2
ATP8B1    LOC100220790    XM_030259059.2
FECH      LOC100223728    XM_032744869.1
MALT1     LOC100224698    XM_012577803.3
ALPK2     LOC100227584    XM_030258318.2
ST8SIA3   LOC100229524    XM_030257523.2
WDR7      LOC100232465    XM_030258950.2
NARS1     LOC116806604    XM_032744866.1
TXNL1     LOC116807012    XM_032744865.1
ONECUT2   LOC116807017    XR_004366058.1
LMAN1     LOC116807022    XM_032744872.1
RAX       LOC116807023    XM_032744877.1
GRP       LOC116807024    XM_032744878.1
SEC11C    SEC11C          XM_032744862.1
EOF

cat <<EOF | awk -F',' '{OFS="\t"; print $1, $2, $3}' > W_PAR_genes_geneID_transID_withGeneName.updated.list
LMAN1,LOC116807022,XM_032744872.1
uncharacterized2,LOC115491185,XR_004366070.1
RAX,LOC116807023,XM_032744877.1
GRP,LOC116807024,XM_032744878.1
SEC11C,SEC11C,XM_032744862.1
ZNF532,LOC100218930,XM_030258323.2
MALT1,LOC100224698,XM_012577803.3
ALPK2,LOC100227584,XM_030258318.2
uncharacterized3,LOC116807018,XR_004366060.1
NEDD4L,LOC100217943,XM_030259050.2
uncharacterized4,LOC116807014,XR_004366054.1
ATP8B1,LOC100220790,XM_030259059.2
NARS1,LOC116806604,XM_032744866.1
FECH,LOC100223728,XM_032744869.1
ONECUT2,LOC116807017,XR_004366058.1
ST8SIA3,LOC100229524,XM_030257523.2
WDR7,LOC100232465,XM_030258950.2
TXNL1,LOC116807012,XM_032744865.1
EOF

# Make table with gene names including uncharacterized genes
cat <<EOF | awk -F',' '{OFS="\t"; print $1, $2, $3}' > PAR_genes_geneID_transID_geneName.list
uncharacterized1,LOC116806781,XR_004365834.1
LMAN1,LOC100220399,XM_030258642.2
uncharacterized2,LOC116806731,XR_004365759.1
RAX,RAX,NM_001243734.2
GRP,LOC105760884,XM_030258647.2
SEC11C,LOC116806596,XM_030257196.2
ZNF532,LOC116806602,XM_032744080.1
MALT1,LOC116806749,XM_032744259.1
ALPK2,LOC116806817,XM_032744489.1
uncharacterized3,LOC116806818,XR_004365896.1
MIR122,MIR122,NR_049051.1
NEDD4L,LOC116806603,XM_032744495.1
uncharacterized4,LOC115491286,XR_003957248.2
uncharacterized5,LOC116806819,XR_004365897.1
ATP8B1,LOC116806743,XM_032744197.1
NARS1,NARS1,XM_030259060.2
FECH,LOC116806605,XM_030259064.2
ONECUT2,LOC100226718,XR_003957244.2
ST8SIA3,LOC116806742,XM_032744196.1
uncharacterized6,LOC115491263,XR_003957222.2
WDR7,LOC116806851,XM_032744613.1
TXNL1,LOC100218989,XM_032744615.1
uncharacterized7,LOC116806852,XR_004365928.1
EOF

###############################################################################
#                     EXTRACT W PAR TRANSCRIPTS                               #
###############################################################################
oneline_fasta() {
  awk '/^>/ {printf("\n%s\n",$0);next;} { printf("%s",$0);} END {printf("\n");}' "$1"
}

gtf_path="GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.small.clean.longestTranscript.gtf"
exon_fasta="longestTranscripts.exons.fa"

grep "NW_022611471.1" "$gtf_path" |
  awk '{print $12}' |
  sort -u |
  tr -d '";' |
  while read gene; do
    oneline_fasta "$exon_fasta" | grep "$gene" -A1
  done > longestTranscripts.exons.PAR.W.fa

###############################################################################
#                               BLAST SEARCH                                  #
###############################################################################

module load bioinfo-tools blast

makeblastdb -in "$exon_fasta" -parse_seqids -dbtype nucl

blastn \
  -db "$exon_fasta" \
  -query longestTranscripts.exons.PAR.W.fa \
  -outfmt 6 \
  > zebrafinch_PAR_W_blast_results.out

```

## 1.2 Synteny analyses to the *Sylvioidea* reference genomes

```bash
###############################################################################
#                                 METADATA                                    #
###############################################################################

cat <<EOF | awk -F' ' '{OFS="\t"; print $1, $2}' > samples_ref_genome.list
AcrSch Acrocephalus_schoenobaenus
AegCau Aegithalos_caudatus
AlaArv skylark_min1kb
CetCet Cettia_cetti
CisJun Cisticola_juncidis
HirDau Hirundo_daurica
LocLus Locustella_luscinioides
PanBia panurus_min1kb
PhyCol Phylloscopus_collybita
PycBar Pycnonotus_barbatus
SylAtr_1EV02922	Sylvia_atricapilla_1EV02922
SylBra Sylvietta_brachyura
TurAlt Turdoides_altirostris
EreAlp horned_lark_min1kb
EOF

cat <<EOF | awk -F',' '{OFS="\t"; print $1, $2, $3}' > samples_sex_sameline.tsv
QF-1504-CP59475_S11_L004,QF-1504-BL37630_S12_L004,AegCau
QL-1681-19_S46_L006,QL-1681-21_S47_L006,AlaArv
QF-1504-P182137_S9_L003,QF-1504-2L18122_S4_L002,CetCet
QF-1504-CISJUN-2_S6_L002,QF-1504-RA5680_S5_L002,CisJun
QF-1504-P182141_S2_L001,QF-1504-P182142_S1_L001,HirDau
QF-1504-LOCLUS-43_S1_L001,QF-1504-LOCLUS-24_S3_L001,LocLus
QF-1504-2KR32024_S2_L001,QF-1504-1ET92164_S3_L001,PanBia
QF-1504-R86159_S5_L002,QF-1504-Z81303_S4_L002,PhyCol
1EL38952_S2_L001,1EV02922_S4_L002,SylAtr_1EV02922
QF-1504-CT90325_S17_L006,QF-1504-CT90312_S18_L006,AcrSch
QF-1504-H-19_S8_L003,QF-1504-H-88_S7_L003,EreAlp
SJ-2333-IB-2b_S32_L002,SJ-2333-IB-1a_S31_L002,TurAlt
SJ-2333-Pbar-197_S24_L002,SJ-2333-Pbar-421_S22_L002,PycBar
SJ-2333-Sbra-553_S28_L002,SJ-2333-Sbra-878_S26_L002,SylBra
EOF

while read -r sp ref; do
  awk -v species="$sp" -v refgen="$ref" '$3 == species { print $0, refgen }' \
    samples_sex_sameline.tsv
done < samples_ref_genome.list |
sed 's/ /\t/g' > samples_sex_sameline_ref.tsv
```

#### Kraken base configuration template (`kraken_base_config_bTaeGut2.pat.W.v2`):

```bash
[genomes]
sylvioidea    QUERY.fasta
zf    GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.chrZ.fasta

[pairwise-maps]
zf sylvioidea    satsuma/bTaeGut1.pat.W.v2/SPECIES/satsuma_summary.chained.out
```

#### Kraken base script (`kraken_Z_base.sh`)

```bash
#! /bin/bash -l 
#
#SBATCH -p core -n 2
#SBATCH -t 4:00:00  
#SBATCH -A snic2020-5-33 -J kraken_SPECIES
target="GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.chrZ.fasta"
query="QUERY"
output="satsuma/bTaeGut1.pat.W.v2/SPECIES"

mkdir -p kraken/bTaeGut1.pat.W.v2/SPECIES
rm kraken/bTaeGut1.pat.W.v2/SPECIES/mapped.gtf

~/bin/kraken/bin/RunKraken -c kraken_SPECIES_config_bTaeGut2.pat.W.v2 -s GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.small.clean.gtf -S zf -T sylvioidea -o kraken/bTaeGut1.pat.W.v2/SPECIES/mapped.gtf 
```

#### Prepare Kraken config for each species

```bash
cat samples_ref_genome.list | \
while read sp ref; do
  cat kraken_base_config_bTaeGut2.pat.W.v2 | sed "s/SPECIES/${sp}/" | sed "s|QUERY|${ref}|" > kraken_${sp}_config_bTaeGut2.pat.W.v2
done
```

#### Satsuma + kraken base script (`satsuma_Z_base.bTaeGut2.pat.W.v2.sbatch`)

```bash
#! /bin/bash -l
#
#SBATCH -p core -n 18
#SBATCH -t 1-12:00:00  
#SBATCH -A snic2020-5-33 -J satsuma_SPECIES
target="GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.chrZ.fasta"
query="QUERY"
output="satsuma/bTaeGut1.pat.W.v2/SPECIES"

~/bin/satsuma-code-0/SatsumaSynteny -t $target -q $query -o $output -n 17

mkdir -p kraken/bTaeGut1.pat.W.v2/SPECIES

~/bin/kraken/bin/RunKraken -c kraken_SPECIES_config_bTaeGut2.pat.W.v2 -s GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.edit.gtf -S zf -T sylvioidea -o kraken/bTaeGut1.pat.W.v2/SPECIES/mapped.gtf
```

#### Run analyses

```bash

cat samples_ref_genome.list | \
while read sp ref; do
  cat satsuma_Z_base.bTaeGut2.pat.W.v2.sbatch | sed "s/SPECIES/${sp}/g" | sed "s|QUERY|${ref}.fasta|g" > satsuma_Z_${sp}.bTaeGut2.pat.W.v2.sbatch
  sbatch satsuma_Z_${sp}.bTaeGut2.pat.W.v2.sbatch
done
```

## 1.3 Synteny analysis to the flycatcher genome

```bash
#! /bin/bash -l 
#
#SBATCH -p core -n 2
#SBATCH -t 4:00:00  
#SBATCH -A snic2020-5-33 -J satsuma_FicAlb
target="data/external_raw/genome/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.chrZ.fasta"
query="data/external_raw/genome/FicAlb1.5.editNames.fasta"
output="intermediate/satsuma/bTaeGut1.pat.W.v2/FicAlb"

~bin/SatsumaSynteny -t $target -q $query -o $output -n 17

mkdir -p intermediate/kraken/bTaeGut1.pat.W.v2/FicAlb

~/bin/kraken/bin/RunKraken -c code/kraken_FicAlb_config_bTaeGut2.pat.W.v2 -s data/external_raw/genome/GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.small.clean.gtf -S zf -T sylvioidea -o intermediate/kraken/bTaeGut1.pat.W.v2/FicAlb/mapped.gtf 
```