# KmerCamel🐫
[![KmerCamel test](https://github.com/OndrejSladky/kmercamel/actions/workflows/ci.yml/badge.svg)](https://github.com/OndrejSladky/kmercamel/actions/)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Prerequisites](#prerequisites)
* [Getting started](#getting-started)
  * [Installation](#installation)
  * [Compression for k-mer set storage](#compression-for-k-mer-set-storage)
  * [k-mer set indexing](#k-mer-set-indexing)
* [Detailed instructions](#detailed-instructions)
  * [Arguments](#arguments)
  * [Converting k-mer set superstring representation to the (r)SPSS representations](#converting-k-mer-set-superstring-representation-to-the-rspss-representations)
* [How it works](#how-it-works)
* [How to test](#how-to-test)
* [Issues](#issues)
* [Changelog](#changelog)
* [Licence](#licence)
* [Contact](#contact)

<!-- vim-markdown-toc -->

## Introduction

KmerCamel🐫 is a tool for efficiently representing a set of k-mers a [masked superstring](https://doi.org/10.1101/2023.02.01.526717).

It is based on the following paper:

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

See [supplementary materials](https://github.com/karel-brinda/masked-superstrings-supplement) of the aforementioned paper for experimental results with KmerCamel🐫.

The computation of masked superstring using KmerCamel🐫 is done in two steps -
first a superstring is computed with its default mask and then its mask can be optimized.

The computation of the masked superstring works as follows. KmerCamel🐫 reads an input FASTA file (optionally `gzip`ed), retrieves the associated k-mers (with supported $k$ up to 127), and outputs
a fasta file with a single record - a masked-cased superstring, which is in the nucleotide alphabet with case of the letters determining the mask symbols.
KmerCamel🐫 implements two different algorithms for computing the superstring:
global greedy and local greedy. Global greedy produces more compact superstrings and therefore is the default option,
but local greedy requires less memory and hence can be more suitable in use cases where memory is the main limitation.

To compute masked superstrings takes about 4-6s / 1M k-mers, which means about 3h to compute masked superstrings for the human genome. The memory consumption on human genome is about 115 GB.

All algorithms can be used to either work in the unidirectional model or in the bidirectional model
(i.e. treat $k$-mer and its reverse complement as the same; in this case either of them appears in the result).

Additionally, KmerCamel🐫 can optimize the mask of the superstring via the `optimize`subcommand. The implemented mask optimization algorithms are the following:
- Minimize the number of 1s in the mask.
- Maximize the number of 1s in the mask.
- Minimize the number of runs of 1s in the mask.

## Prerequisites

* GCC
* Zlib
* GLPK (can be installed via `apt-get install libglpk-dev` on Ubuntu or `brew install glpk` on macOS)

## Getting started

### Installation

Download and compile KmerCamel🐫 by running the following commands:

```
git clone --recursive https://github.com/OndrejSladky/kmercamel
cd kmercamel && make
```

Alternatively, you can install KmerCamel from Bioconda:
```
   conda install kmercamel
```

### Compression for k-mer set storage

```
kmercamel ms -k 31 -o ms.msfa yourfile.fa             # Compute MS with the default mask
kmercamel msfa2ms -m mask.m -s superstring.s ms.msfa  # Extract superstring and mask
bzip2 --best mask.txt
xz -T1 -9 superstring.txt
```

For a super efficient compression of the superstring (often <2 bits / bp), you use some of the specialized tools based on statistical compression such as [GeCo3](https://github.com/cobilab/geco3) or [Jarvis3](https://github.com/cobilab/jarvis3).


### k-mer set indexing

Example with [FMSI](https://github.com/OndrejSladky/fmsi/activity?ref=main):
```
kmercamel ms -k 31 -o ms.msfa -M maxonemask.m yourfile.fa          # Compute MS and the maxone mask
kmercamel msfa2ms -m /dev/null -s superstring.s ms.msfa            # Extract superstring
kmercamel ms2msfa -m maxonemask.m -s superstring.s -o ms-opt.msfa  # Combine with maxone mask
fmsi index -p ms-opt.msfa                                          # Create a k-mer index
```

## Detailed instructions

Examples of computing masked superstrings (`ms` subcommand):
```
kmercamel ms -k 31 yourfile[.fa|.fa.gz] -o ms.msfa         # From a (gziped) fasta file, use "-" for stdin
kmercamel ms -k 31 -u yourfile.fa -o ms.msfa               # Treat k-mer and its reverse complement as distinct
kmercamel ms -k 31 -M maxonemask.m yourfile.fa -o ms.msfa  # Also store mask with maximum ones
kmercamel ms -k 31 -a local yourfile.fa -o ms.msfa         # Use local instead of global for lower memory footprint (likely worse result)
```

Examples of optimizing masks:
```
kmercamel optimize -a maxone -k 31 ms.msfa -o ms-opt.msfa    # Maximize the number of 1s in the mask
kmercamel optimize -a minone -k 31 ms.msfa -o ms-opt.msfa    # Minimize the number of 1s in the mask
kmercamel optimize -a minrun -k 31 ms.msfa -o ms-opt.msfa    # Minimize the number of runs of consecutive 1s in the mask.
```

Format conversions:
```
kmercamel ms2msfa -m dataset.m -s dataset.s -o dataset.msfa  # M and S -> mask-cased MS in msfa
kmercamel msfa2ms -m dataset.m -s dataset.s dataset.msfa     # Mask-cased MS -> M and S
kmercamel spss2msfa -k 31 -o dataset.msfa dataset.rspss      # rSPSS/general fasta to its corresponding MS
kmercamel msfa2spss -k 31 -o dataset.fa dataset.fa           # Splitting MS in msfa into rSPSS in fa
```

Compute lower bound on the minimum possible superstring length of a k-mer set:
```
./kmercamel lowerbound -p -k 31 yourfile.fa
```

To view all options for a particular subcommand, run `kmercamel <subcommand> -h`.

Additionally, KmerCamel🐫 experimentally implements both algorithms in their Aho-Corasick automaton versions. To use them, add `AC` to the algorithm name.
Note that they are much slower than the original versions, but they can handle arbitrarily large *k*s.

## How it works

For details about the algorithms and their implementation, see the [Code README](./src/README.md).

## How to test

To ensure correctness of the results, KmerCamel🐫 has two levels of tests - unit tests and file-specific integration tests.

For integration tests  install [jellyfish (v2)](https://github.com/gmarcais/Jellyfish)
and add it to PATH.

You can verify all the algorithms for `1 < k < 128` on a *S. pneumoniae* by running `make verify`.
To run it on another dataset, see the [verification script](./verify.py).

You can run the C++ unittests by `make cpptest`.

To run all the test, simply run `make test`.

## Issues

Please use [Github issues](https://github.com/OndrejSladky/kmercamel/issues).


## Changelog

See [Releases](https://github.com/OndrejSladky/kmercamel/releases).


## Licence

[MIT](https://github.com/OndrejSladky/kmercamel/blob/master/LICENSE.txt)


## Contact

[Ondrej Sladky](https://iuuk.mff.cuni.cz/~sladky/) \<ondra.sladky@gmail.com\>\

