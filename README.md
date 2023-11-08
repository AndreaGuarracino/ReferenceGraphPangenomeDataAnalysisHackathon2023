# Reference Graph Pangenome Data Analysis Hackathon 2023

Material for the **Reference Graph Pangenome Data Analysis Hackathon 2023**, November 13-17 in Cape Town, South Africa.
The event was sponsored by [H3ABioNet](https://www.h3abionet.org/).

## Table of Contents

* [Reference-based pangenome graph building](#reference-based-pangenome-graph-building)
  * [Constructing and viewing your first graphs](#constructing-and-viewing-your-first-graphs)
* [Small pangenome graph building from sequence alignments](#small-pangenome-graph-building-from-sequence-alignments)
  * [Building HLA pangenome graphs](#building-hla-pangenome-graphs)
  * [Building LPA pangenome graphs](#building-lpa-pangenome-graphs)
* [Human pangenome graph building from sequence alignments](#human-pangenome-graph-building-from-sequence-alignments)
* [Understanding pangenome graphs](#understanding-pangenome-graphs)


## Reference-based pangenome graph building

### Learning objectives

- construct graphs using `vg construct`
- visualize graphs using `vg view` and `Bandage`
- convert graphs using `vg convert`

### Getting started

Make sure you have `vg` installed.
It is already available on the course workstations.
In this exercise, you will use small toy examples from the `test` directory.
So make sure you have checked out `vg` repository:

    cd ~
	git clone https://github.com/vgteam/vg.git

Now create a directory to work on for this tutorial:

    cd ~
	mkdir ref_based_pangenome_graph_building
	cd ref_based_pangenome_graph_building
	ln -s ~/vg/test/tiny

### Constructing and viewing your first graphs

First we will use `vg construct` to build our first graph.
Run it without parameters to get information on its usage:

	vg construct

We will construct a reference-based graph from the sequence in file `tiny/tiny.fa`, which looks like this:

	>x
	CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG

To construct a simple graph, run:

	vg construct -r tiny/tiny.fa -m 32 > tiny.ref.vg

This will create a graph that just consists of a linear chain of nodes, each with 32 characters.

The `-m` option tells `vg` to put at most 32 characters into each graph node.

To visualize a graph, you can use `vg view`.
By default, `vg view` will output a graph in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format.

	vg view tiny.ref.vg

Try to run `vf construct` with different values of `-m` and observe the different results.

Now let's build a new graph that has some variants built into it.
First, take a look at `tiny/tiny.vcf.gz`, which contains variants in (gzipped) [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format.

    zgrep '^##' -v tiny/tiny.vcf.gz | column -t
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO                           FORMAT  1
    x       9    .   G    A    99    .       AC=1;LEN=1;NA=1;NS=1;TYPE=snp  GT      1|0
    x       10   .   C    T    99    .       AC=2;LEN=1;NA=1;NS=1;TYPE=snp  GT      1|1
    x       14   .   G    A    99    .       AC=1;LEN=1;NA=1;NS=1;TYPE=snp  GT      1|0
    x       34   .   T    A    99    .       AC=2;LEN=1;NA=1;NS=1;TYPE=snp  GT      1|1
    x       39   .   T    A    99    .       AC=1;LEN=1;NA=1;NS=1;TYPE=snp  GT      1|0

then use `vg construct` to build the graph:

	vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz -m 32 > tiny.vg

We can write the graph in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format:

    vg view tiny.vg > tiny.gfa

In GFA format, each line is a separate record of some part of the graph.
The lines come in several types, which are indicated by the first character of the line.
What line types do you see? What do you think they indicate?

<details>
  <summary>Click me for the answers</summary>

You should see the following line types:
- `H`: a header.
- `S`: a "sequence" line, which is the sequence and ID of a node in the graph.
- `L`: a "link" line, which is an edge in the graph.
- `P`: a "path" line, which labels a path of interest in the graph. In this case, the path is the walk that the reference sequence takes through the graph.

Of note, the format does not specify that these lines come in a particular order.
</details>

A for visualizing (not too big) pangenome graphs is [BandageNG](https://github.com/asl/BandageNG).
It supports graphs in GFA format.
To use Bandage, just:

    wget -c https://github.com/asl/BandageNG/releases/download/v2022.09/BandageNG-9eb84c2-x86_64.AppImage
    chmod +x BandageNG-9eb84c2-x86_64.AppImage

Download the graph on your computer and try to visualize it locally with:

    ./BandageNG-9eb84c2-x86_64.AppImage load tiny.gfa 

Click the `Drawn graph` button to visualize the graph.
Click the 'Length Node' button under Node labels feature.
Select the longest node (it is in the middle, 19 bp) by clicking on it, then click the `Paths...` button to see the paths that pass through the selected node.
Check that the paths displayed are the same as you expect by checking the GFA. Any problems?

<details>
  <summary>Click me for the answers</summary>

The GFA format does not specify that its lines have to follow a particular order.
However, the downloaded version of Bandage seems to have a bug in reading the paths properly.
Edit the GFA file by moving the `P` at the end of the file.
Now `Bandage` should show correctly the paths traversing the selected node.
Never trust the tools!
</details>

Try to generate a PNG image with `Bandage` (take a look at the `BandageNG-9eb84c2-x86_64.AppImage -h` output).
This can be helpful to take at graphs that are too big to be directly loaded with `Bandage`.

Now, let's build a new complex graph.
First, take a look at `tiny/multi.vcf.gz`, which contains multiallelic variants. 
Try to generate with that input a graph in GFA format and visualize the output with `Bandage`. What are you able to see?


## Small pangenome graph building from sequence alignments

### Learning objectives

- build pangenome graphs using `pggb`
- explore `pggb`'s results
- understand how parameters affect pangenome graph building

### Getting started

Make sure you have `pggb` and its tools installed.
It is already available on the course workstations.
If you want to build everything on your laptop, follow the instructions at the [pggb homepage](https://github.com/pangenome/pggb#installation).
So make sure you have checked out `pggb` repository:

    cd ~
	git clone https://github.com/pangenome/pggb.git

Check out also `wfmash` repository (we need one of its scrips):

    cd ~
    git clone https://github.com/waveygang/wfmash.git

Now create a directory to work on for this tutorial:

    cd ~
	mkdir alignment_based_pangenome_graph_building_small
	cd alignment_based_pangenome_graph_building_small
	ln -s ~/pggb/data

### Building HLA pangenome graphs

The [human leukocyte antigen (HLA)](https://en.wikipedia.org/wiki/Human_leukocyte_antigen) system is a complex of genes on chromosome 6 in humans which encode cell-surface proteins responsible for the regulation of the immune system.

Let's build a pangenome graph from a collection of sequences of the DRB1-3123 gene:

    pggb -i data/HLA/DRB1-3123.fa.gz -o out_DRB1_3123.1 -n 12

Why did we specify `-n 12`?

<details>
  <summary>Click me for the answer</summary>

This parameter is important for the graph normalization with `smoothxg`.
It is used to determine the right partial order alignment (POA) problem size for the multiple sequence alignments.
</details>

How many pairwise alignments were used to build the graph (take a look at the `PAF` output)? Visualize the alignments:

    cd out_DRB1_3123.1
    ~/wfmash/scripts/paf2dotplot png large *paf

The last command will generate a `out.png` file with a visualization of the alignments.
Purple lines indicate that the 2 sequences are aligned in the same orientation.
Blue lines indicate that the 2 sequences are aligned in different orientation.

Take a look at the files in the `out_DRB1_3123.1` folder.

- `*.alignments.wfmash.paf`: sequence alignments;
- `*.log`: whole log;
- `*.params.yml`: `pggb`'s parameters in `YAML` format;
- `*.gfa`: final pangenome graph in GFA format;
- `*.og`: final pangenome graph in ODGI format (used in `odgi`);
- `*.lay`: graph layout in `LAY` format (used in `odgi`);
- `*.draw_multiqc.png`: static graph layout representation with sequences as colored lines;
- `*.draw.png`: static graph layout representation;
- `*.lay.tsv`: graph layout in `TSV` format;


`*.viz_*.png` images represent the graph in 1 dimension, with sequences represented as colored bars and graph links represented as black lines at the bottom.
Each image follow a different color scheme:
- `*.viz_multiqc.png`: each path has a different color, without any meaning;
- `*.viz_depth_multiqc.png`: paths are colored by depth. We define **node depth in a path** as the number of times the node is crossed by a path; 
- `*.viz_inv_multiqc.png`: paths are colored with respect to the strandness (black for forward, red for reverse);
- `*.viz_O_multiqc.png`: all paths are compressed into a single line, where we color by path coverage;
- `.viz_pos_multiqc.png`: paths are colored with respect to the node position in each path. Smooth color gradients highlight well-sorted graphs.
- `*.viz_uncalled_multiqc.png`: uncalled bases (`Ns`) are colored in green.

Try to visualize the graph also with `Bandage`.

Use `odgi stats` to obtain the graph length, and the number of nodes, edges, and paths:

    cd ~/alignment_based_pangenome_graph_building_small
    odgi stats -i out_DRB1_3123.1/DRB1-3123.fa.gz.bf3285f.eb0f3d3.9c6ea4f.smooth.final.og -S

Do you think the resulting pangenome graph represents the input sequences well?
Check the length and the number of the input sequences to answer this question.
To answer, check the length of the input sequences.

<details>
  <summary>Click me for the answer</summary>

The input sequences are ~13.6Kbp long, on average.
The graph is about 1.6X longer, so not much longer, then it is a good representation of the input sequences.
Pangenome graphs longer than the input sequences are expected because they contain the input sequences plus their variation.
</details>

`pggb`'s default parameters assume an average divergence of approximately 10% (`-p 90` by default).
Try building the same pangenome graph by specifying a higher percent identity

    cd ~/alignment_based_pangenome_graph_building_small
    pggb -i data/HLA/DRB1-3123.fa.gz -o out_DRB1_3123.2 -n 12 -p 95

Check the graph statistics.
Does this pangenome graph represent better or worse the input sequences than the previously produced graph?

<details>
  <summary>Click me for the answer</summary>

The graph is much longer than before, about ~4.1X longer than the input sequences.
This indicates a strong under-alignment of all the sequences.
This happens because the HLA locus is highly polymorphic in the population, with great genetic variability.
</details>

Try to increase and decrease the segment length (`-s 5000` by default):

    cd ~/alignment_based_pangenome_graph_building_small
    pggb -i data/HLA/DRB1-3123.fa.gz -o out_DRB1_3123.3 -n 12 -s 15000
    pggb -i data/HLA/DRB1-3123.fa.gz -o out_DRB1_3123.4 -n 12 -s 100

How is this affecting graph statistics?

<details>
  <summary>Click me for the answer</summary>
This parameter influences the sensitivity in detecting structural variants (SVs) and inversions.
Lower values lead to better resolution of SVs breakpoints and the possibility of detecting shorter inversions, but at the same time increase the complexity of the graph in terms of the number of nodes and edges.
This happens because short segment lengths lead to catching shorter homologies between the input sequences (that is, more mappings and then alignments).
Higher values reduce sensitivity, but lead to simpler graphs.
</details>

Choose another HLA gene from the `data` folder (`A-3105.fa.gz` for example) and explore how the statistics of the resulting graph change as` s` and `p` change.

### Building LPA pangenome graphs

[Lipoprotein(a) (LPA)](https://en.wikipedia.org/wiki/Lipoprotein(a)) is a low-density lipoprotein variant containing a protein called apolipoprotein(a).
Genetic and epidemiological studies have identified lipoprotein(a) as a risk factor for atherosclerosis and related diseases, such as coronary heart disease and stroke.

Try to make LPA pangenome graphs.
The input sequences are in `data/LPA/LPA.fa.gz`.
Sequences in this locus have a peculiarity: which one?
Hint: visualize the alignments and take a look at the graph layout with `Bandage` and/or in the `*.draw_multiqc.png` files.
The `*.draw_multiqc.png` files contain static representations of the graph layout.
They are similar to what `Bandage` shows (probably a little less attractive), but such visualizations can scale to larger pangenomic graphs.


## Human pangenome graph building from sequence alignments

### Learning objectives

- collect and preprocess *de novo* human assemblies
- partition the assembly contigs by chromosome
- built chromosome-specific pangenome graphs

### Getting started

Create a directory to work on for this tutorial:

    cd ~
	mkdir human_pangenome_graphs
	cd human_pangenome_graphs

Download 2 human references and 4 diploid human *de novo* assemblies from the Human Pangenome Reference Consortium (HPRC) data:

    mkdir -p ~/human_pangenome_graphs/assemblies
    cd ~/human_pangenome_graphs/assemblies

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.maternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.maternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.maternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.maternal.f1_assembly_v2_genbank.fa.gz

Decompress and index the assemblies:

    gunzip *genbank.fa.gz
    ls *genbank.fa | while read f; do echo $f; samtools faidx $f; done

When making pangenome graphs, you should always include at least one reference genome in order to use it as a coordinate system (for projecting variants and/or exploiting its annotations).
We usually use both GRCh38 and CHM13.

### Pangenome Sequence Naming

We follow the [PanSN-spec](https://github.com/pangenome/PanSN-spec) naming to simplify the identification of samples and haplotypes in pangenomes.
The HPRC samples already follow such a convention (`1` is the PATERNAL haplotype, `2` is the MATERNAL haplotype), but we need to apply PanSN naming to the reference genomes.
So, let's add a prefix to their sequence names.
We can do that by using [fastix](https://github.com/ekg/fastix):

    fastix -p 'grch38#1#' <(zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz) | bgzip -@ 16 > grch38_full.fa.gz
    samtools faidx grch38_full.fa.gz
    
    fastix -p 'chm13#1#' <(zcat chm13v2.0.fa.gz) | bgzip -@ 16 > chm13.fa.gz
    samtools faidx chm13.fa.gz

About GRCh38, we remove the unplaced contigs that are (hopefully) properly represented in CHM13:

    samtools faidx grch38_full.fa.gz $(cut -f 1 grch38_full.fa.gz.fai | grep -v _ ) | bgzip -@ 16 > grch38.fa.gz

Cleaning:

    rm chm13v2.0.fa.gz GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz grch38_full.fa.gz.*

Take a look at how sequence names are changed in the FASTA files.

### Sequence partitioning

To reduce analysis complexity, we partition assembly contigs by chromosome and generate chromosome-specific pangenome graphs.
For doing that, we first need to put the two reference genomes together

    zcat chm13.fa.gz grch38.fa.gz | bgzip -@ 16 > chm13+grch38.fa.gz && samtools faidx chm13+grch38.fa.gz

and then map each assembly against the two reference genomes:

    DIR_BASE=~/human_pangenome_graphs
    PATH_REFERENCES_FASTA=$DIR_BASE/assemblies/chm13+grch38.fa.gz

    mkdir -p $DIR_BASE/assemblies/partitioning

    ls $DIR_BASE/assemblies/*.f1_assembly_v2_genbank.fa | while read FASTA; do
      NAME=$(basename $FASTA .fa);
      echo $NAME

      PATH_PAF=$DIR_BASE/assemblies/partitioning/$NAME.vs.ref.paf
      wfmash $PATH_REFERENCES_FASTA $FASTA -m -N -t 32 > $PATH_PAF
    done

`wfmash` should take 4-5 minutes for each haplotype (each FASTA file).

Why are we using two reference genomes? 

<details>
  <summary>Click me for the answer</summary>
A single human genome cannot fully represent the genetic variability of the entire population.
Having multiple reference genomes allows greater genetic variability to be represented, allowing contigs to be better partitioned.
</details>

Recently, a new high-quality human diploid assembly was released at [HG002](https://github.com/marbl/HG002).
You could also try using these 2 new haplotypes (`HG002#MATERNAL` and `HG002#PATERNAL`) to partition assembly contigs and see if you can partition more.

Run `wfmash` without parameters to get information on the meaning of each parameter.

What does `-m` mean?

<details>
  <summary>Click me for the answer</summary>
We ask `wfmash` to compute only the mapping, not the alignment, to save computation time.
To partition chromosomes, we don't need the base-level alignments.
</details>

What does `-N` mean?

<details>
  <summary>Click me for the answer</summary>
With `-N` we ask `wfmash` to generate mappings that completely cover the assembly contigs.
Without `-N`, mappings can only cover parts of the contigs, and different parts of the same contig could map to difference reference chromosomes.
This is frequent for contigs belonging to acrocentric chromosomes or sex chromosomes because of their homologous regions (https://www.nature.com/articles/s41586-023-05976-y, https://doi.org/10.1093/hmg/7.13.1991).
</details>

For each haplotype (for each FASTA file), count how many contigs were partitioned in each reference chromosome.
Which of the two references (GRCh38 and CHM13) do more contigs map to? If there is a clear winner, why?

<details>
  <summary>Click me for the answer</summary>
CHM13 is a complete human genome (GRCh38 is almost complete, but not 100%), so we expect to be able to map more assembly contigs to it.
</details>

What chromosome do most contigs map to in GRCh38 and CHM13? and which chromosome has the least number of contigs mapped to GRCh38 and CHM13? Why?

<details>
  <summary>Click me for the answer</summary>
In GRCh38 the shor arms of the acrocentric chromosomes (chromosome 13, 14, 15, 21 and 22) are missing, so we struggle to align acrocentric contigs to it.
These short arms are available in CHM13 and indeed we are able to map lots of contigs against them.
</details>

It should be noted that also with `wfmash -N`, there can be cases with contigs fully mapping to different cromosomes. For example:

    grep 'HG00438#2#JAHBCA010000147.1' *paf
      HG00438.maternal.f1_assembly_v2_genbank.vs.ref.paf:HG00438#2#JAHBCA010000147.1  738336  0       738336  +       chm13#1#chr13 113566686       8582789 9321125 245     738336  23      id:f:99.4899    kc:f:0.057824
      HG00438.maternal.f1_assembly_v2_genbank.vs.ref.paf:HG00438#2#JAHBCA010000147.1  738336  0       738336  +       chm13#1#chr22 51324926        5299094 6037430 245     738336  23      id:f:99.4899    kc:f:0.057824

For which there is not enough information to determine which is the best chromosome to map against (acrocentric chromosomes are hard).
For these case, we just randomly take one result.

    DIR_BASE=~/human_pangenome_graphs
    ls $DIR_BASE/assemblies/*.f1_assembly_v2_genbank.fa | while read FASTA; do
      NAME=$(basename $FASTA .fa);
      echo $NAME

      PATH_PAF=$DIR_BASE/assemblies/partitioning/$NAME.vs.ref.paf
      cut -f 1,6 $PATH_PAF | sed -e 's/chm13#1#//g' -e 's/grch38#1#//g' | awk '{
          contig = $1;
      
          # If the contig is not already in the data array, add it
          if (!(contig in data)) {
              data[contig] = $2;
          }
      }
      END {
        # Output the result
        for (contig in data) {
          print contig "\t" data[contig];
        }
      }' > $DIR_BASE/assemblies/partitioning/$NAME.vs.ref.assignments.tsv
    done

Now we can subset assembly contigs by chromosome:

    DIR_BASE=~/human_pangenome_graphs
    ( seq 1 22; echo X; echo Y ) | while read i; do
      echo chr$i
      
      grep chr$i --no-filename -w $DIR_BASE/assemblies/partitioning/*.assignments.tsv | cut -f 1 > $DIR_BASE/assemblies/partitioning/chr$i.contigs.txt
    done

Then, we create a FASTA file for each chromosome, reference chromosomes included.
To save time and space, let's take only sequences from chromosome 20:

    DIR_BASE=~/human_pangenome_graphs
    ( echo 20 ) | while read i; do
      echo chr$i
      samtools faidx $DIR_BASE/assemblies/chm13+grch38.fa.gz chm13#1#chr$i grch38#1#chr$i > $DIR_BASE/assemblies/partitioning/chr$i.fa
      
      ls $DIR_BASE/assemblies/*.f1_assembly_v2_genbank.fa | while read FASTA; do
        NAME=$(basename $FASTA .fa);
        echo $NAME

        samtools faidx $FASTA $( comm -12 <(cut -f 1 $FASTA.fai | sort) <(sort $DIR_BASE/assemblies/partitioning/chr$i.contigs.txt) ) >> $DIR_BASE/assemblies/partitioning/chr$i.fa
      done

      bgzip -@ 16 $DIR_BASE/assemblies/partitioning/chr$i.fa
      samtools faidx $DIR_BASE/assemblies/partitioning/chr$i.fa.gz
    done

Check that everything went fine:

    head $DIR_BASE/assemblies/partitioning/chr20.fa.gz.fai | column -t
      chm13#1#chr20                66210255  15         60  61
      grch38#1#chr20               64444167  67313791   60  61
      HG00438#2#JAHBCA010000018.1  36199104  132832057  60  61
      HG00438#2#JAHBCA010000074.1  688514    169634509  60  61
      HG00438#2#JAHBCA010000089.1  29856295  170334528  60  61
      HG00438#2#JAHBCA010000131.1  448468    200688457  60  61
      HG00438#2#JAHBCA010000161.1  865084    201144429  60  61
      HG00438#2#JAHBCA010000193.1  130419    202023961  60  61
      HG00438#2#JAHBCA010000210.1  103918    202156583  60  61
      HG00438#1#JAHBCB010000023.1  26194192  202262262  60  61

### Pangenome graph building

Build the pangenome graph for chromosome 20.

    DIR_BASE=~/human_pangenome_graphs
    mkdir -p $DIR_BASE/graphs
    cd $DIR_BASE/graphs

    pggb -i $DIR_BASE/assemblies/partitioning/chr20.fa.gz -o $DIR_BASE/graphs/pggb.chr20 -p 98 -k 49 -D /scratch -t 32

This should take less than 1 hour.

**IMPORTANT**: The `-D` parameter in `pggb` is used to specify the directory used for temporary files.
This directory should be on a high-speed disk (like an SSD) to avoid severe slowdowns during graph construction.

Note that we are not specifying the number of haplotypes (`-n` parameter).
Indeed, `pggb` can automatically obtain this information thanks to the fact that the input sequences respect the PanSN naming.

Since human presents a low sequence divergence, we set the mapping identity (`-p` parameter) in `pggb` to `98`.
Additionally, we specify `-s 10k` to get a simpler and more linear graph structure, which is easier to work with.
Lower values lead to more sensitive mappings but to the possibility of having circular graphs due to the similarity of the telomeres.

The `-k` parameter is used to filter exact matches in the sequence alignments shorter than 49 bps.
Indeed, graph induction with `seqwish` often works better when we filter short matches out of the input alignments.
In practice, these often occur in regions of low alignment quality, which are typical of areas with large indels and structural variations.
Removing short matches can simplify the graph and remove spurious relationships caused by short repeated homologies.
However, this filter might lead to under-alignment, that is resolved in the graph normalization step with `smoothxg`.
By default, `pggb` uses `-k 23`.
With human data, higher values work well.

Use `odgi stats` to obtain the graph length, and the number of nodes, edges, and paths.
Do you think the resulting pangenome graph represents the input sequences well?
Check the length and the number of the input sequences to answer this question.

<details>
  <summary>Click me for the answer</summary>
chr20 is approximately 64-66 Mbp long in reference genomes.
There is no right value because it always depends on how much genetic variability we want to represent in the graph.
Lower `-s` and `-p` values will tend to represent more homologies, getting shorter, but more complex, graphs.
</details>

Take a look at the PNG files in the `pggb.chr20` folder. Is the layout of the graph roughly linear?

Generate another `odgi viz` visualization with

    cd $DIR_BASE/graphs/pggb.chr20
    odgi paths -i chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og -L | cut -f 1,2 -d '#' | uniq > prefixes.txt
    odgi viz -i chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og -o chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og.viz_multiqc.2.png -x 1500 -y 500 -a 10 -M prefixes.txt

What do you think is different between the `chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og.viz_multiqc.png` image 
and the newly generated image (`chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og.viz_multiqc.2.png`)?

<details>
  <summary>Click me for the answer</summary>

The contigs of the two additional haplotypes are visually collapsed into one path each.
</details>
</br>
