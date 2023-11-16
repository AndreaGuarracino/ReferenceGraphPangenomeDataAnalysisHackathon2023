# Reference Graph Pangenome Data Analysis Hackathon 2023

## Understanding H3ABioNet pangenome graphs

You can find the `pggb` graphs for all chromosomes here: `/cbio/projects/031/andreaguarracino/6samples`.
Each graph was made with 6 diploid samples (so 12 haplotypes) plus 2 reference genomes (`GRCh38` and `CHM13`).

### Chromosome 6 pangenome graph

We extract the MHC locus and color bars by haplotype.
 
![chr6 MHC locus](images/chr6.pan.MHC.png)

The MHC locus is a bit broken:
- the worst sample, `CMI0Q`, covers the locus with multiple contigs;
- the best sample, `FPV0R`, covers the locus with 1 contig for each haplotype.

Let's take a look at the C4 locus.
In 1D, we color by copy number:
- white: 0 copies
- grey: 1 copy
- red: 2 copies
- orange: 3 copies

![chr6 C4 locus 1D](images/chr6.pan.C4.sorted.m.png)

The `EFB0A#1` haplotype has 3 copies of the C4 genes.

The graph layout is the [grafiocavallo](https://en.wikipedia.org/wiki/Caciocavallo) that we expect:

![chr6 C4 locus 2D](images/chr6.pan.C4.sorted.2D.png)

The big loop represents the copy number variation, while the nested loop differentiates the short and long forms of the C4 genes.

If we inject `GRCh38` annotations and untangle the graph, we can see a bit closer the C4 variation:

![chr6 C4 locus untangling](images/chr6.pan.C4.untangling.png)

The locus is inverted in the graph (not a problem).
