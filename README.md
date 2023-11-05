# Reference Graph Pangenome Data Analysis Hackathon 2023

Material for the **Reference Graph Pangenome Data Analysis Hackathon 2023**, November 13-17 in Cape Town, South Africa.
The event was sponsored by [H3ABioNet](https://www.h3abionet.org/).

## Reference-based pangenome graph building

### Learning objectives

- construct graphs using `vg construct`
- visualize graphs using `vg view` and `Bandage`
- convert graphs using `vg convert`

### Getting started

Make sure you have `vg` installed. It is already available on the course workstations.
In this exercise, you will use small toy examples from the `test` directory.
So make sure you have checked out `vg` repository:

    cd ~
	git clone https://github.com/vgteam/vg.git

Now create a directory to work on for this tutorial:

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
</br>

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
</br>

Try to generate a PNG image with `Bandage` (take a look at the `BandageNG-9eb84c2-x86_64.AppImage -h` output).
