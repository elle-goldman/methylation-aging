<!--
RECOVERED FROM BACKUP — 2026-09-23

Source: /Volumes/EBU Time Machine snapshot 2026-06-01-071816,
        ~/Desktop/Projects/CR/AnnotationWorkflow2.md

This file no longer exists on the working disk; Desktop/Projects/CR, /Aging and
/Rhesus aging were all removed some time after 2026-06-01. Recovered verbatim.

Two caveats, left uncorrected so nothing is misattributed:
  * The shell snippets have spaces inserted around every flag ("bedtools
    intersect - a" rather than "-a", "sort - k1.4, 1.5h" rather than
    "-k1.4,1.5h"). They will not run as written — the mangling looks like a
    paste artefact, not intent.
  * The final "CpG Density Calculation" R block is an empty stub.
-->

# Workflow for Annotating Significant CR Diet × Precise Age at Start DMRs
Hypothesis:The metabolomic signature of CR is differentially affected by sex and age at diet start (total number of years on diet and number of post -
adolescent years on diet)
Hypothesis 2:A GLP - 1 signature of CR / weight loss should affect more upstream processes related to satiety as well as pancreatic function compared to a signature of adopting a CR diet, which should look more like preservation of more youthful cellular functionality.
## LiftOver Summary

DMRs from rheMac10 were liftOver - ed to hg38 using UCSC liftOver browser to create `liftedOver_rheMac10_hg38_CRdietXageDMRs.bed`.

To create `relaxed_parameters_liftOver_RheMac10_hg38_CRdietXpreciseAgeAtStart.bed`, the **
minimum ratio of bases that must remap ** was lowered to ** 0.8 ** and **
minimum hit size in query ** set to ** 10 ** .

## LiftOver Parameters Used

| Parameter | Value Used |
| ----------- | ------------ |
| Minimum ratio of bases that must remap | 0.90 |
| Minimum hit size in query | 20 |
| Minimum chain size in target | 50 |
|
Min ratio of alignment blocks / exons that must map (BED12 only) |
0.90 |
| Allow multiple output regions | Unchecked |

## Annotation Preprocessing

** Input ** :`relaxed_parameters_liftOver_RheMac10_hg38_CRdietX preciseAgeAtStart.bed`

** Command ** :#1. Get rid of deletions
```bash
awk '$4 != "deleted"' relaxed_parameters_liftOver_RheMac10_hg38_CRdietXpreciseAgeAtStart.bed |
\
sort - k1.4, 1.5h - k2, 2h > sorted_relaxed_parameters_liftOver_RheMac10_hg38_CRdietXpreciseAgeAtStart.bed
cp sorted_relaxed_parameters_liftOver_RheMac10_hg38_CRdietXpreciseAgeAtStart.bed human_dmrs_sorted.bed
```
## obtaining files using batch download fro ENCODE
The “files.txt” file can be copied to any server.
The following command using cURL can be used to download all the files in the list:xargs -
n 1 curl - O - L < files.txt

# 2. Sort & annotate
sortBed - i human_dmrs.bed > human_dmrs_sorted.bed
bedtools slop - i cpg_islands.bed - g hg38.chrom.sizes - b 2000 > cpg_shores.bed
bedtools intersect - a human_dmrs_sorted.bed - b cpg_islands.bed - c > dmrs_cpg_islands.txt
bedtools intersect - a human_dmrs_sorted.bed - b enhancers.bed - c > dmrs_enhancers.txt
bedtools closest - a human_dmrs_sorted.bed - b atac_peaks.bed - D ref > dmrs_atac_dist.txt

## CpG Density Calculation

```{r}




```