Space-efficient variable-order Markov models
=========

Tools for building and using space-efficient representations of variable-order Markov models and of interpolated Markov models, that support a large number of context-selection criteria, scoring functions, probability smoothing methods, and interpolations, while taking less space than implementations based on tries, suffix trees, and suffix arrays. This software can handle multi-gigabyte training and query datasets.


Requirements
------------

* A modern, C++11 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
* The [cmake][cmake] build system.
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.


Installing and testing
------------

Install the [sdsl-lite library](https://github.com/simongog/sdsl-lite) at the project root:

```
git clone https://github.com/simongog/sdsl-lite
cd sdsl-lite
sh install.sh
cd ..
```

Install [BD_BWT_index](https://github.com/jnalanko/BD_BWT_index):

```
cd BD_BWT_index
cmake -DCMAKE_BUILD_TYPE=Release .
make
cd ..
```

or:

```
cd BD_BWT_index
cmake -DCMAKE_BUILD_TYPE=Debug . 
make
cd ..
```

Compile the rest with `make`:

```
make tests
make optimized
```

The `tests` executable runs the test suite, and takes approximately 10 minutes.



Building and rebuilding models
------------

The `build_model_optimized` executable builds a variable-order Markov model from a file. Example:

```
mkdir models
./build_model_optimized --reference-raw data.txt --entropy 0.2 --outputdir models --maxreps-pruning --rle
```

Full list of flags:

* `--reference-raw [file path]` Assumes that every byte in the input file is a character of a single input string, i.e. that the input file contains exactly one string with no header.

* `--reference-fasta [file path]` Assumes that the input file is in multi-FASTA format, i.e. that every line is either a FASTA header or it contains part of a string. Use this flag to train a Markov model from a set of strings rather than from a single string.
    
* `--outputdir [directory path]` Where to write the built model. This directory must exist before running! The model consists of a set of files such that the filename of each model file is prefixed by the filename. This means that if you build models for two files with the same filename into the same output directory, then the latter model will overwrite the former.
     
* `--maxreps-pruning` Enables maxrep pruning.
    
* `--rle` Enables run-length encoding.
    
* `--depth [integer depth]` Enables depth pruning to the given depth. **Also enables maxrep pruning**.
   
* `--entropy [float threshold]` Use entropy-style contexts with the given threshold.
   
* `--KL [float threshold]` Use Kullback–Leibler-style contexts with the given threshold.
    
* `--pnorm [integer p] [float threshold]` Use p-norm-style contexts with the given threshold.
    
* `--four-thresholds [float tau1] [float tau2] [float tau3] [float tau4]` Use the context formula with the four thresholds *tau1*, *tau2*, *tau3*, *tau4* from the paper "[A framework for space-efficient string kernels][KERNELSPAPER]".
    
* `--store-depths` Stores the string depth of every maximal repeat in the topology as binary integers into `outputdir + "/" + filename_prefix + ".string_depths"`. The binary representation has length that is just enough to store the largest depth.    The file is created even if the option is not enabled, but in that case it will be very small.

* `--context-stats` Computes statistics on the contexts. Writes two files into the model directory:
  * `stats.context_summary.txt`: number of context candidates and number of contexts.
  * `stats.depths_and_scores.txt`: one line for each context: `[string depth] [tree depth] [score(s)]`. The score(s) are:
    * In case of `--KL`, `--entropy` or `--pnorm`, the value that is compared against the threshold. 
    * In case of `--four-thresholds`, there are three values corresponding to equations 2,3 and 4 in the paper "[A framework for space-efficient string kernels][KERNELSPAPER]".

If a model had already been built for a string, program `reconstruct_optimized` rebuilds just its contexts part. Example:

```
./reconstruct_optimized --file data.txt --entropy 10 --dir models
```

Full list of flags:

* `--dir [directory path]` The directory where the model is.

* `--file` The filename of the input data from which the model was built. Only the
    filename, not the full path. That is, if the data is at `./foo/bar/data.txt`,
    give only `data.txt`. This is needed so that the code knows the prefix of the 
    model files.
    
* `--entropy [float threshold]` Use entropy-style contexts with the given threshold.
   
* `--KL [float threshold]` Use Kullback–Leibler-style contexts with the given threshold.
    
* `--pnorm [integer p] [float threshold]` Use p-norm-style contexts with the given threshold.
    
* `--four-thresholds [float tau1] [float tau2] [float tau3] [float tau4]` Use the context formula with the four thresholds, as described above.

* `--context-stats` Computes statistics on the contexts, as described above.



Computing the score of a query
---------

Program `score_string_optimized` scores a string given a model that had been previously built with `build_model_optimized`. It writes log-probabilities to `stdout`, one line per query, and it writes to `stderr` a progress report while it's running. Example:

```
./score_string_optimized --query-raw queries.txt --dir models --file data.txt --escapeprob 0.05
```

Full list of flags:

* `--query-raw [file path]` Assumes that every byte in the input file is a character of a single query string, i.e. that the input file contains exactly one string with no header.

* `--query-fasta [file path]` Assumes that the input file is in multi-FASTA format, i.e. that every line is either a FASTA header or it contains part of a query string. Use this flag to compute a distinct score for every string in a set.

Takes in the queries in fasta-format. Each read is one query. Not tested very well.
      
* `--dir [directory path]` Directory where the model is stored.
    
* `--file [filename]` Filename of the reference string. Only the filename, not the full path. 
    That is, if the data is at `./foo/bar/data.txt`, give only `data.txt`. This 
    is needed so that the code knows the prefix of the model files.

* `--escapeprob [float prob]` Escape probability used in scoring.

* `--recursive-fallback` Uses the recursive scoring method from the paper "[A framework for space-efficient string kernels][KERNELSPAPER]".

* `--lin-scoring` Uses the scoring method from the paper "[Probabilistic suffix array: efficient modeling and prediction of protein families][SAPAPER]".


[KERNELSPAPER]: https://link.springer.com/article/10.1007/s00453-017-0286-4 "A framework for space-efficient string kernels"
[SAPAPER]: https://academic.oup.com/bioinformatics/article/28/10/1314/211256 "Probabilistic suffix array: efficient modeling and prediction of protein families"
[PREZZA]: https://github.com/nicolaprezza/lz-rlbwt
[cmake]: http://www.cmake.org/ "CMake tool"

<!---
If there is a problem with some of the flags maybe I updated the flags but forgot
to update this documentation, or maybe I typoed something. In this case please check
the main-function in score_string.cpp to see what the flags really are and how they 
are parsed.    
    
If there is a problem with some of the flags maybe I updated the flags but forgot
to update this documentation, or maybe I typoed something. In this case please check
the main-function in build_model.cpp to see what the flags really are and how they 
are parsed.
-->


External resources used
---------

This repository includes:

* portions of code from [lz-rlbwt][PREZZA] by Nicola Prezza, for the run-length-encoded BWT;

* the [BD_BWT_index](https://github.com/jnalanko/BD_BWT_index) by Jarno Alanko, for the bidirectional BWT index used during construction.


Related software
---------


