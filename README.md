## Requirement
==============
Make sure that you have the following installed

1.  Python2.7 (Recommend Anaconda) and Packages:
    1. Networkx
    2. Biopython (version >=1.68)
    3. Numpy (version >= 1.11.2)
    4. Keras (version >=1.1.2)
2.  Diamond (version >= 0.8.18)
3.  MetaGeneMark
4.  CD-hit
5.  IDBA-UD


## Installation
===============

$ git clone https://github.com/Rinoahu/POEM

$ cd ./POEM

$ bash ./install.sh

## Example
===============

example directory contain a genome fasta file of E.coli, run the script named runme.sh to test the pipeline

$ cd ./example

$ bash ./runme.sh genome



## Usage
===============

For short reads:

    $ bash ./bin/run_poem.sh -f reads.fsa -a y

    reads.fsa is single fasta file. If the reads are paired-end files in fastq or fasta format, 
    use the fq2fa of IDBA_UD to convert them to  a single fasta file.

For genome/assembly:

    $ bash ./bin/run_poem.sh -f genome.fsa -a n

    genome.fsa is the genome/assembly fasta file.


## Output
===============

POEM will create a directory named read.fasta_output to save the results. The results include serveral file:

    1.  input.fsa:
        Contig of IDBA-UD output

    2.  input.fsa_gmk_aa.fsa and input.fsa.gff:
        MetaGeneMark output of protein sequence and gff file on step 1

    3.  input.fsa.cdhit and input.fsa.cdhit.clstr:
        CD-hit outputs of step 2

    4.  input.fsa_aa.fsa:
        Filtered protein sequence according to the result of step 3

    5.  input.fsa.flt.blast8:
        m8 tabular output of step 4 against COG database

    6.  input.fsa.flt.blast8.flt:
        Filtered blast output of step 5. Only hits which identity >= 80 are retained

    7.  input.fsa.cog:
        COG annotation for protein sequence of step 4

    8.  input.fsa.locus:
        The file records the gene's locus on genome

    9.  input.fsa.adjacency and input.fsa.operon:
        Predicted operonic adjacency and full operon of step 8

    10. input.fsa.cog_adjacency:
        This file contain two part:
        1). The COG adjacency of operonic adjacency
        2). Operon's core COG and its closest real operon

    11. input.fsa.core_cog_adjacency:
        The core COG adjacency

    12. input.fsa.core_network.sif, input.fsa.core_node.attrs and input.fsa.core_edge.attrs:
        Network, node attribute and edge attribute extracted from step 11 for cytoscape visualization





