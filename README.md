# Reference Graph Pangenome Data Analysis Hackathon 2023

Material for the **Reference Graph Pangenome Data Analysis Hackathon 2023**, November 13-17 in Cape Town, South Africa.
The event was sponsored by [H3ABioNet](https://www.h3abionet.org/).

## Table of Contents

* [Reference-based pangenome graph building](#reference-based-pangenome-graph-building)
  * [Constructing and viewing your first graphs](#constructing-and-viewing-your-first-graphs)
* [Small pangenome graph building from sequence alignments](#pangenome-graph-building-from-sequence-alignments)
  * [Building HLA pangenome graphs](#building-hla-pangenome-graphs)
  * [Building LPA pangenome graphs](#building-lpa-pangenome-graphs)


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

Select the longest node (it is in the middle) by clicking on it, then click the `Paths...` button to see the paths that pass through the selected node.
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

Take a look at the files in the `out_DRB1_3123.1` folder. Visualize the graph with `Bandage`.

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

Choose another HLA gene from the `data` folder and explore how the statistics of the resulting graph change as` s` and `p` change.

### Building LPA pangenome graphs

[Lipoprotein(a) (LPA)](https://en.wikipedia.org/wiki/Lipoprotein(a)) is a low-density lipoprotein variant containing a protein called apolipoprotein(a).
Genetic and epidemiological studies have identified lipoprotein(a) as a risk factor for atherosclerosis and related diseases, such as coronary heart disease and stroke.

Try to make LPA pangenome graphs.
The input sequences are in `data/LPA/LPA.fa.gz`.
Sequences in this locus have a peculiarity: which one?
Hint: visualize the alignments and take a look at the graph layout (with `Bandage` and/or in the `.draw_multiqc.png` files).


