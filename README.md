lz-rlbwt-sparse: Run-Length Compressed Burrows-Wheeler transform with Sparse LZ77 suffix array sampling
===============
Author: Nicola Prezza (nicolapr@gmail.com)

From the paper: "Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza, and Mathieu Raffinot. Composite repetition-aware data structures"

### Brief description

The lz-rlbwt-sparse data structure is a fast index that takes advantage of text repetitions in order to reduce its memory footprint. Highly repetitive inputs produce longer BWT runs (RLE compressed using an Elias-Fano-encoded RLBWT data structure) and a smaller LZ77 parse. The lz-rlbwt-sparse exploits this fact and substitutes SA samples with components from a LZ77 index (basically, a sampled suffix tree and geometric range search data structures). 

As opposed to the original lz-rlbwt, this index also uses sparsification on the LZ77 parsing

### Download

To clone the repository, use the --recursive option:

> git clone --recursive http://github.com/nicolaprezza/lz-rlbwt-sparse

### Compile

The library has been tested under linux using gcc 4.9.2. You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

lz-rlbwt-sparse uses cmake to generate the Makefile. Create a build folder in the main lz-rlbwt-sparse folder:

> mkdir build

execute cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  lz-rlbwt-sparse-build input.txt

This command will create the lz-rlbwt-sparse index of the text file "input.txt" and will store it using as prefix "input.txt". Use option -o to specify a different basename for the index files. At the moment, building the index takes a lot of time since we use a slow but memory-efficient LZ77 parser.

Run

> lz-rlbwt-sparse-count index_basename patterns_file

or

> lz-rlbwt-sparse-locate index_basename patterns_file

to count/locate a set of patterns in pizza&chilli format (http://pizzachili.dcc.uchile.cl/experiments.html) using the lz-rlcsa index. At the moment, these executables do not store the output anywhere (we use it only for benchmarking purposes)
