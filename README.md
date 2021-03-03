pafr demo
================

# Comparative genomics in fungi

``` bash
bash scripts/download_genomes.sh
echo 'export "MM_CPU=5"' > .shell_vars
```

    ## data/seqs/A_flav_SU-16.fna.gz: OK
    ## data/seqs/A_flav_CA14.fna.gz: OK
    ## data/seqs/kea.fna.gz: OK
    ## data/seqs/kakapo.fna.gz: OK

``` bash
source .shell_vars
if [ ! -f data/ali/aflav_ali.paf ];
then
    minimap2 -t ${MM_CPU} -cx asm5 data/seqs/A_flav_CA14.fna.gz data/seqs/A_flav_SU-16.fna.gz > data/ali/aflav_ali.paf
fi
```

``` r
library(pafr)

asper <- read_paf("data/ali/aflav_ali.paf")
asper
```

    ## pafr object with 702 alignments (36.5Mb)
    ##  9 query seqs
    ##  9 target seqs
    ##  12 tags: NM, ms, AS, nn, tp, cm, s1, s2, de, zd, rl ...

``` r
asper_clean <- filter_secondary_alignments(asper)
dotplot(asper_clean, order_by="qstart", label_seqs=TRUE)
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
plot_synteny(asper, t_chrom="CP061809.1", q_chrom="CP047251.1", rc=TRUE)
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
inv <- subset(asper, tname=="CP061809.1" & qname=="CP047251.1" & strand=="+")
inversion_intervals <-  query_bed(inv)
head(inversion_intervals)
```

    ## pafr object with 6 alignments (0Mb)
    ##  1 query seqs
    ##  0 target seqs
    ##  -6 tags:

``` r
write_bed(inversion_intervals, "data/results/inversion_query.bed")
```

``` sh
bedtools merge -d 5000 -i data/results/inversion_query.bed
```

    ## CP047251.1   4408713 4409218
    ## CP047251.1   4508061 4883353

# Comparing draft and finished genome assemblies

``` bash
source .shell_vars
if [ ! -f data/ali/parrots.paf.gz ]
then
    minimap2 -t ${MM_CPU} -cx asm20 data/seqs/kea.fna.gz data/seqs/kakapo.fna.gz | gzip > data/ali/parrots.paf.gz
fi
```

``` r
parrots <- read_paf("data/ali/parrots.paf.gz")
hq_parrots <- subset(filter_secondary_alignments(parrots), mapq > 50)
hq_parrots
```

    ## pafr object with 57638 alignments (1090.8Mb)
    ##  77 query seqs
    ##  34903 target seqs
    ##  11 tags: NM, ms, AS, nn, tp, cm, s1, s2, dv, zd, cg

``` r
dotplot(hq_parrots, order_by="qstart", dashes=FALSE ,label_seqs=TRUE) 
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
kea_wing <- "#D73202"
ggplot(hq_parrots, aes(alen/tlen)) + 
    geom_histogram(fill=kea_wing, colour='black', binwidth=0.05) + 
    scale_x_continuous()
```

![](README_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
subset(hq_parrots, alen /tlen > 0.99)
```

    ## pafr object with 20274 alignments (619.1Mb)
    ##  43 query seqs
    ##  20042 target seqs
    ##  11 tags: NM, ms, AS, nn, tp, cm, s1, s2, dv, zd, cg

``` r
sex_chrom_map <- c("CM013773.2" = "Z", "CM013763.2" = "W")
sex_chrom_ali <- subset(hq_parrots, qname %in% names(sex_chrom_map))
sex_chrom_ali$qname <- sex_chrom_map[sex_chrom_ali$qname]
dotplot(sex_chrom_ali, order_by="qstart" ,label_seqs=TRUE, dashes=FALSE)
```

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)
