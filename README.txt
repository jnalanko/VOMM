Compiling:

* Install sdsl-lite at the project root:

git clone https://github.com/simongog/sdsl-lite
cd sdsl-lite
sh install.sh
cd ..

* Install BD_BWT_index

cd BD_BWT_index
cmake -DCMAKE_BUILD_TYPE=Release .
make
cd ..

or

cd BD_BWT_index
cmake -DCMAKE_BUILD_TYPE=Debug . 
make
cd ..

* Compile rest with make

make tests
make optimized

---------

There are three executables:

* tests

Runs the test suite. Might take 10 minutes.

* build_model_optimized

Builds a variable order markov model from a file. Example:

mkdir models
./build_model_optimized --reference-raw data.txt --entropy 0.2 --outputdir models --maxreps-pruning --rle

Full list of flags:

--reference-fasta [file path]
    Takes the input file in the fasta-format. Not tested very well :)
    
--reference-raw [file path]
    Takes the input in a raw text file
    
--outputdir [directory path]
    Where to write the built model. This directory must exist before running!
    The model consists of a set of files such that the filename of each model file
    is prefixed by the filename. This means that if you build models for two 
    files with the same filename into the same output directory, then the latter 
    model will overwrite the former.
     
--maxreps-pruning
    Enables maxrep pruning
    
--rle
    Enables run length encoding
    
--depth [integer depth]
    Enables depth pruning to the given depth. ALSO ENABLES MAXREP PRUNING.
   
--entropy [float threshold]
    Use entropy-style contexts with the given threshold
   
--KL [float threshold]
    Use Kullback–Leibler-style contexts with the given threshold
    
--pnorm [integer p] [float threshold]
    Use p-norm-style contexts with the given threshold
    
--four-thresholds [float tau1] [float tau2] [float tau3] [float tau4]
    Use the context formula with the four thresholds tau1,tau2,tau3,tau4
    
--context-stats
    Computes statistics on the contexts. Writes three files 
    into the model directory:
    - stats.context_summary.txt: 
      Number of context candidates and number of contexts
    - stats.depths_and_scores.txt
      One line for each context: [string depth] [tree depth] [score(s)]
      The score(s) are:
       In case of --KL, --entropy or --pnorm, the values that is compared 
       against the threshold In case of --four-thresholds, there are three 
       values corresponding to equations 2,3 and 4 in the paper A 
       Framework for Space-Efficient String Kernels.

    
If there is a problem with some of the flags maybe I updated the flags but forgot
to update this documentation, or maybe I typoed something. In this case please check
the main-function in build_model.cpp to see what the flags really are and how they 
are parsed.

* score_string_optimized

Scores a string given a previously built model with build_model_optimized.

Example:

./score_string_optimized --query-raw queries.txt --dir models --file data.txt --escapeprob 0.05

Full list of flags:

--query-fasta [file path]
    Takes in the queries in fasta-format. Each read is one query. Not tested very well.
    
--query-raw [file path]
    Takes in single raw text file as one big query.
    
--dir [directory path]
    Directory where the model is stored
    
--file [filename]
    Filename of the reference string. This is needed to load the correct model
    because model files are prefixed by the filename of the reference string.

--escapeprob [float prob]
    Escape probability used in scoring.
    
If there is a problem with some of the flags maybe I updated the flags but forgot
to update this documentation, or maybe I typoed something. In this case please check
the main-function in score_string.cpp to see what the flags really are and how they 
are parsed.



