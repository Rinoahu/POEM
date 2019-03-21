## Requirements
==============
Make sure that you have the following installed

1.  Python2.7 (Recommend [Anaconda](https://www.continuum.io/downloads#linux "https://www.continuum.io/downloads#linux" ) ) and Packages:
    1. Networkx
    2. Biopython (version >=1.68)
    3. Numpy (version >= 1.11.2)
    4. Keras (version >=1.1.2)
2.  [Diamond](https://github.com/bbuchfink/diamond "https://github.com/bbuchfink/diamond") (version >= 0.8.18) 
3.  [MetaGeneMark](http://exon.gatech.edu/Genemark/ "http://exon.gatech.edu/Genemark")
4.  [CD-hit](https://github.com/weizhongli/cdhit "https://github.com/weizhongli/cdhit")
5.  [IDBA-UD](https://github.com/loneknightpy/idba "https://github.com/loneknightpy/idba")



## Installation
===============

$ git clone https://github.com/Rinoahu/POEM

$ cd ./POEM

$ bash ./install.sh


## Usage
===============

For short reads:

    $ bash ./bin/run_poem.sh -f reads.fsa -a y

    reads.fsa is single fasta file. If the reads are paired-end files in fastq or fasta format, 
    use the fq2fa of IDBA_UD to convert them to  a single fasta file.

For genome/assembly:

    $ bash ./bin/run_poem.sh -f genome.fsa -a n

    genome.fsa is the genome/assembly fasta file.



## Example
===============

example directory contain a genome fasta file of E.coli, run the script named runme.sh to test the pipeline

$ cd ./example

$ bash ./runme.sh genome



## Output
===============

The output of POEM is a set files:

    1.  input.fsa:
        Contigs from IDBA-UD output

    2.  input.fsa_gmk_aa.fsa and input.fsa.gff:
        Proein sequences and gff file predicted by MetaGeneMark on contigs from step 1

    3.  input.fsa.cdhit and input.fsa.cdhit.clstr:
        CD-hit clustering results for protein sequences from step 2

    4.  input.fsa_aa.fsa:
        Filtered protein sequence according to the results from step 3

    5.  input.fsa.flt.blast8:
        -m 8 tabular file of protein sequences from step 4 against the COG database

    6.  input.fsa.flt.blast8.flt:
        Filtered blast results from step 5. Only hits which identity >= 80 are kept

    7.  input.fsa.cog:
        COG annotation for protein sequences from step 4

    8.  input.fsa.locus:
        A file lists loci for all genes.

    9.  input.fsa.adjacency and input.fsa.operon:
        Predicted operonic adjacency and full operon of step 8

    10. input.fsa.cog_adjacency:
        This file contain two part:
        1). The COG adjacency of operonic adjacency
        2). Operon's core COG and its closest real operon

    11. input.fsa.core_cog_adjacency:
        The core COG adjacency

    12. input.fsa.core_network.sif, input.fsa.core_node.attrs and input.fsa.core_edge.attrs:
        These files are used for cytoscape visualization of core COG adjacency from step 11. input.fsa.core_network.sif records relationship between nodes. input.fsa.core_node.attrs records colors and names of nodes. input.fsa.core_edge.attrs records weight of edge between each adjacent nodes.





