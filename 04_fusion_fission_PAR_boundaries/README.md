# 4. Identification of fusion, fission and PAR boundary scaffolds

## 4.1 Synteny analyses between Sylvioidea species and zebra finch

```bash
# Synteny analysis with lastal between Sylvioidea and zebra finch
cat data/meta/samples_sex_sameline_ref_outgroup.tsv | while read female male sp ref ; do lastal data/external_raw/genome/bTaeGut2.pat.W.v2.db data/external_raw/genome/${ref}.fasta -P 15 | last-split > intermediate/lastal_ZF_bTaeGut1.pat.W.v2/${sp}_align
maf-convert psl intermediate/lastal_ZF_bTaeGut1.pat.W.v2/${sp}_align > intermediate/lastal_ZF_bTaeGut1.pat.W.v2/${sp}_align_converted
done

# Make genome windows files
cat data/meta/samples_sex_sameline_ref_outgroup.tsv | while read female male sp ref ; do cat intermediate/lastal_ZF_bTaeGut1.pat.W.v2/${sp}_align_converted | awk '{print $10,$12,$13,$14,$16,$17,$1}' | sed 's/ /\t/g' | bedtools intersect -a stdin -b intermediate/bedtools_17nov2019/${sp}_ref_${sp}/genome_5kb_windows.out -wa -wb | awk '{if($10-$9=="5000") print $8,$9,$10,$7,$1,$2,$3,$4,$5,$6}' | sed 's/ /\t/g' | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' > intermediate/synteny_bTaeGut1.pat.W.v2/${sp}.genome_windows.out ; done

# Find best match for each window
cat data/meta/samples_sex_sameline_ref_outgroup.tsv | while read female male sp ref ; do mkdir intermediate/synteny_bTaeGut1.pat.W.v2/${sp}_temp
cd intermediate/synteny_bTaeGut1.pat.W.v2/${sp}_temp
awk '{print >> $1; close($1)}' ../${sp}.genome_windows.out
ls | while read file; do cat $file | awk -v max=0 '{if($2>max){want=$0; max=$2}}END{print want}' ; done | sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' >> ../${sp}_bestMatch.list
cd /proj/sllstore2017102/nobackup/hanna/sylvioidea_sexchromosome
rm -r intermediate/synteny_bTaeGut1.pat.W.v2/${sp}_temp
done

# Here I create files with the synteny to the zebra finch chromosomes
ls | grep converted | sed 's/_align_converted//' | grep -v GRW | while read sp ; do cat ${sp}_align_converted | awk '$1 > 600 {print}' | cut -f 10- | grep NC_045027.1 | awk '$8<1000000 {print}' | cut -f 1 | sort | uniq | awk '{print NR, $1}' | sed 's/ /\t/' | while read nr scaff ; do cat ${sp}_align_converted | awk '$1 > 600 {print}' | grep $scaff | cut -f 10,14 | awk '{print "'"$nr"'" "\t" $1 "\t" $2 "\t" "'"$sp"'"}' ; done ; done > PAR_synteny.out

# Removing single matches
cat PAR_synteny.out | sort | uniq -c | awk '$1>1 {print $3,$4,$5}' | sed $'s/ /\t/g' | while read scaff chr sp ; do cat PAR_synteny.out | grep $scaff | grep $chr | grep $sp ; done > PAR_synteny_no_singles.out

ls | grep converted | sed 's/_align_converted//' | grep -v GRW | while read sp ; do cat ${sp}_align_converted | awk '$1 > 1000 {print}' | cut -f 10- | grep NC_045027.1 | awk '$8<1000000 {print}' | cut -f 1 | sort | uniq | awk '{print NR, $1}' | sed 's/ /\t/' | while read nr scaff ; do cat ${sp}_align_converted | awk '$1 > 1000 {print}' | grep $scaff | cut -f 10,14 | awk '{print "'"$nr"'" "\t" $1 "\t" $2 "\t" "'"$sp"'"}' ; done ; done > PAR_synteny_1kb.out

# Get all matches + match positions 
ls | grep converted | sed 's/_align_converted//' | grep -v GRW | while read sp ; do cat ${sp}_align_converted | cut -f 10- | grep NC_045027.1 | awk '$8<1000000 {print}' | cut -f 1 | sort | uniq | awk '{print NR, $1}' | sed 's/ /\t/' | while read nr scaff ; do cat ${sp}_align_converted | grep $scaff | cut -f 1,10,12,14,16 | awk '{print "'"$nr"'",$0,"'"$sp"'"}' | sed 's/ /\t/' ; done ; done > PAR_synteny_all_matches.out

ls | grep converted | sed 's/_align_converted//' | grep -v GRW | while read sp ; do cat ${sp}_align_converted | cut -f 10- | grep NC_045027.1 | awk '$8<1000000 {print}' | cut -f 1 | sort | uniq | awk '{print NR, $1}' | sed 's/ /\t/' | while read nr scaff ; do cat ${sp}_align_converted | grep $scaff | cut -f 1,10,12,13,14,16,17 | awk '{print "'"$nr"'",$0,"'"$sp"'"}' | sed 's/ /\t/' ; done ; done > PAR_synteny_all_matches_new.out

cat PAR_synteny_1kb.out | sort | uniq -c | awk '$1>1 {print $3,$4,$5}' | sed $'s/ /\t/g' | while read scaff chr sp ; do cat PAR_synteny_1kb.out | grep $scaff | grep $chr | grep $sp ; done > PAR_synteny_1kb_no_singles.out

# Bundle links
cd results
mkdir PAR_circos
cat PAR_synteny_all_matches_new.out | sed $'s/ /\t/' | cut -f 9 | sort | uniq | while read sp ; do cat PAR_synteny_all_matches_new.out | sed $'s/ /\t/' | awk '$2>200 && $9=="'"$sp"'" {print $3,$4,$5,$6,$7,$8}' | sed $'s/ /\t/g' > PAR_circos/links_${sp}.out ; done
cd PAR_circos/

# Will only consider last matches that fall within these regions
cat ../PAR_synteny_all_matches_new.out | sed $'s/ /\t/' | cut -f 9 | sort | uniq | while read sp ; do ~/bin/circos-tools-0.23/tools/bundlelinks/bin/bundlelinks -links links_${sp}.out -max_gap 10000 -min_bundle_size 2000 | awk '{print "'"$sp"'",$0}' ; done > allSp_bundlelinks_maxGap10kb_minBundleSize2kb.out

cat allSp_bundlelinks_maxGap10kb_minBundleSize2kb.out | sed $'s/ /\t/g' | cut -f 1-7 | while read sp scaff sStart sEnd chr cStart cEnd ; do cat ../PAR_synteny_all_matches_new.out | sed $'s/ /\t/' | grep $chr | grep $scaff | grep $sp  | awk -F "[\t]+"  -v MINSCAFF=$sStart -v MAXSCAFF=$sEnd -v MINCHR=$cStart -v MAXCHR=$cEnd '($4 >= MINSCAFF) && ($5 <= MAXSCAFF) && ($7 >= MINCHR) && ($8 <= MAXCHR) {print}' ; done > ../PAR_synteny_all_matches_new.bundleFiltered.out
```
