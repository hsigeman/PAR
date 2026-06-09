# 3. PAR gene tree analyses and W-to-Z pairwise distances

## 3.1 Make FASTA files for each exon

```bash
# Define directories and files
krakendir="kraken/bTaeGut1.pat.W.v2"
genes="ZF_PAR.genes.new.list"
home=$PWD

# Combine female and male PAR exon sequences from all species into a single file
cat samples_sex_sameline_ref.tsv | cut -f 3 | while read sp ; do
    cat $krakendir/${sp}/${sp}.PAR.exonSeparate.female.fasta | sed "s/>/>${sp}_female_/"
done > allSp.PAR.exonSeparate.female.male.fasta

cat samples_sex_sameline_ref.tsv | cut -f 3 | while read sp ; do
    cat $krakendir/${sp}/${sp}.PAR.exonSeparate.male.fasta | sed "s/>/>${sp}_male_/"
done >> allSp.PAR.exonSeparate.female.male.fasta

# Generate list of unique PAR exons
cat allSp.PAR.exonSeparate.female.male.fasta | grep ">" | sed 's/male_/\t/' | cut -f 2 | sort | uniq > allPARexons.list

# Function to format fasta files into single-line sequences
oneline_fasta() {
  awk '/^>/ {printf("\n%s\n",$0);next;} { printf("%s",$0);} END {printf("\n");}' "$1"
}

# Separate fasta files by gene and save them
mkdir exonSeparate
cat allSp.PAR.exonSeparate.female.male.fasta | grep ">" | sed 's/male_/\t/' | cut -f 2 | sort | uniq | while read gene ; do
    oneline_fasta allSp.PAR.exonSeparate.female.male.fasta | grep ${gene}$ -A 1 | grep -v "^--" | sed 's/SylAtr_1EV02922/SylAtr/' > exonSeparate/${gene}.fasta
done
```

## 3.2 Make alignments with PRANK

```bash
interactive -A naiss2024-5-340 --cores 7 --mem 50000 --tmp 100
salloc --cores=24 -t 2:00:00 -A naiss2024-5-340 -p main
mkdir exonSeparate/prank2
cd exonSeparate/prank2

module load bioinfo-tools
module load prank/170427
module load parallel

# Generate file list and align using prank
ls ../ | grep fasta$ > files
srun parallel -j 22 'prank -d=../{} -o={}.prank.aln -f=fasta -F ' :::: files
srun parallel -j 22 'prank -d=../{} -o={}.prank.F.aln -f=fasta +F ' :::: files

# Rename sequences in alignment files and make lists of files
mkdir IDfix
ls | grep fasta.prank.aln.best.fas | while read file ; do
    cat $file | sed 's/male_/male\t/' | cut -f 1 > IDfix/${file}
done

# Create file list for each gene in order
cd IDfix
ls | grep aln | sed 's/_/\t/' | cut -f 1 | sort | uniq | while read gene ; do
    ls | grep $gene | sed 's/_/\t/' | sed 's/.fasta.prank.aln.best.fas//' | cut -f 2 | sort -n | awk '{print "'"$gene"'""_"$1".fasta.prank.aln.best.fas"}' | tr '\r\n' ' '  
    echo
done > file_lists_wide.list

# Download catfasta2phyml.pl script here: https://github.com/nylander/catfasta2phyml

# Concatenate sequence files
mkdir concat
ls | grep aln | sed 's/_/\t/' | cut -f 1 | sort | uniq | while read gene ; do
    cat file_lists_wide.list | grep $gene | while read line ; do
        ../../../catfasta2phyml.pl -f $line -c > concat/${gene}.fasta.prank.aln.best.fas
    done
done

# Add zebra finch sequences to the alignment
mkdir mafft_ZF
mkdir TaeGut_sequences
mkdir renamed

cp concat/*fasta.prank.aln.best.fas mafft_ZF/

# Rename sequences and process zebra finch data
cd $home
cat PAR_genes_geneID_transID_geneName.list | while read name gene trans ; do
    cp exonSeparate/prank2/IDfix/concat/${gene}.fasta.prank.aln.best.fas exonSeparate/prank2/IDfix/renamed/${name}.fasta.prank.aln.best.fas
done

cat PAR_genes_geneID_transID_geneName.list | while read name gene trans ; do
    oneline_fasta ../data/external_raw/genome/longestTranscripts.exons.fa | grep $trans -A 1 | sed 's/>/>TaeGut_Z\t/' | cut -f 1 > exonSeparate/prank2/IDfix/TaeGut_sequences/${name}.TaeGut.fa
done

cat W_PAR_genes_geneID_transID_withGeneName.updated.list | while read name gene trans ; do
    oneline_fasta ../data/external_raw/genome/longestTranscripts.exons.fa | grep $trans -A 1 | sed 's/>/>TaeGut_W\t/' | cut -f 1 >> exonSeparate/prank2/IDfix/TaeGut_sequences/${name}.TaeGut.fa
done

# Re-align using MAFFT
cd exonSeparate/prank2/IDfix
ls renamed/ | sed 's/.fasta.prank.aln.best.fas//' | while read gene ; do
    mafft --reorder --add TaeGut_sequences/${gene}.TaeGut.fa --auto renamed/${gene}.fasta.prank.aln.best.fas > mafft_ZF/${gene}.realn.fasta
done
```

## 3.3 Trim alignments with trimAI

```bash

# Trim alignments with TrimAl
module load trimAl/1.4.1
cd mafft_ZF

ls | grep fasta$ | while read gene ; do trimal -in $gene -out ${gene}.auto.trim  -automated1 ; done

# Species-specific trimming

# RAX - removing PanBia (would remove half the sequence if kept)
trimal -in RAX.realn.fasta.auto.trim -selectseqs { 2,3 } > RAX.realn.fasta.auto.selectseqs.trim
iqtree -s RAX.realn.fasta.auto.selectseqs.trim -m TEST -bb 1000 -alrt 1000

# uncharacterized3 - remove PanBia
trimal -in uncharacterized3.realn.fasta.auto.trim -selectseqs { 28,29 } > uncharacterized3.realn.fasta.auto.selectseqs.trim

# uncharacterized4 - remove PycBar (short) and AcrSch (bad alignment)
trimal -in uncharacterized4.realn.fasta.auto.trim -selectseqs { 0,1,4,5 } > uncharacterized4.realn.fasta.auto.selectseqs.trim

# ST8SIA3 - remove SylAtr (short) 
trimal -in ST8SIA3.realn.fasta.auto.trim -selectseqs { 4,5 } > ST8SIA3.realn.fasta.auto.selectseqs.trim

ls | grep fasta$ | while read gene ; do trimal -in $gene.auto.selectseqs.trim -out ${gene}.auto.selectseqs.gt.0.8.trim -gt 0.8 ; done

module load bioinfo-tools
module load iqtree/1.5.3-omp

ls | grep fasta$ | while read gene ; do iqtree -s ${gene}.auto.selectseqs.gt.0.8.trim -m TEST -bb 1000 -alrt 1000 ; done

ls | grep fasta$ | while read gene ; do iqtree -s ${gene}.auto.gt.0.8.trim -m TEST -bb 1000 -alrt 1000 ; done

```

## 3.4 Phylogenetic trees with iqtree

```bash
#! /bin/bash -l 
#
#SBATCH -p core -n 10
#SBATCH -t 23:00:00  
#SBATCH -A snic2020-5-33 -J iqtree

source activate lark

home=$PWD

cd $home
cd results/PAR/genesSeparate_incl_ZF_W/prank

ls | grep resoverlap.seqoverlap.automated.trim$ > alignment_files.automated.list
parallel -j 9 'iqtree -s {} -m TEST -bb 1000 -alrt 1000' :::: alignment_files.automated.list

# = # = #

# Phylogenetic trees with iqtree
#! /bin/bash -l 
#
#SBATCH -p core -n 10
#SBATCH -t 12:00:00  
#SBATCH -A snic2020-5-33 -J iqtree

source activate lark

cd $home
cd results/PAR/genesSeparate_incl_ZF_W/clustalo
ls | grep trim$ > alignment_files.list

parallel -j 9 'iqtree -s {} -m TEST -bb 1000 -alrt 1000' :::: alignment_files.list

cd $home
sbatch code/iqtree_PAR_genes_clustalo_new.sbatch 

# = # = # = # = # = # = # = # = # =# =

```

## 3.5 Calculate Z-to-W sequence distance

```bash
ls | grep realn.fasta.auto.gt.0.8.trim$ | while read file ; do awk '
BEGIN {
    RS=">"; FS="\n"
}
NR > 1 {
    header = $1
    seq = ""
    for (i = 2; i <= NF; i++) seq = seq $i
    gsub(/[ \t\r]/, "", seq)

    prefix = header
    sub(/_.*/, "", prefix)

    names[prefix] = names[prefix] " " header
    seqs[prefix, header] = seq
}
END {
    print "Species\tIdentity_percent\tCompared_headers"
    for (p in names) {
        split(names[p], h, " ")
        if (length(h) == 2) {
            s1 = seqs[p, h[1]]
            s2 = seqs[p, h[2]]

            matches = 0
            total = 0

            for (i = 1; i <= length(s1); i++) {
                b1 = substr(s1, i, 1)
                b2 = substr(s2, i, 1)

                if (b1 != "-" && b2 != "-") {
                    total++
                    if (b1 == b2) matches++
                }
            }

            pct = (total > 0) ? 100 * matches / total : 0
            printf "%s\t%.2f\t%s vs %s\n", p, pct, h[1], h[2]
        }
    }
}
' $file | awk '{print "'"$file"'""\t"$0}'; done | sed 's/.realn.fasta.auto.gt.0.8.trim//' > seqIdentity.tsv

cut -f 1-3 seqIdentity.tsv | grep -v Identity_percent | awk '
{
    genes[$1]
    species[$2]
    data[$1, $2] = $3
}
END {
    # Print header
    printf "Gene"
    for (s in species) printf "\t%s", s
    print ""

    # Print rows
    for (g in genes) {
        printf "%s", g
        for (s in species) {
            key = g SUBSEP s
            printf "\t%s", (key in data ? data[key] : "NA")
        }
        print ""
    }
}
' > seqIdentity_wide.tsv

```
