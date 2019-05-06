## Compiling
GNU g++-7 must be installed.
In order to compile Look4TRs, run the following commands

> cd Look4TRs

> mkdir bin

> cd bin

> cmake ..

> make

## Splitting fasta files
Look4TRs is designed to read in multiple fasta files where each file represents a differnt chromosome.
However, the fasta file format specifies that an entire genome can be placed into a single file.
To remedy this, we have provided a file under src/ called divideGenomeIntoChromosomesEfficient.pl.
divideGenomeIntoChromosomesEfficient.pl will divide the FASTA file into multiple files based on their headers.
To call this, use the following command

> perl divideGenomeIntoChromosomesEfficient.pl Input.fa OutputDirectory

Then, call Look4TRs using `OutputDirectory` as the input directory from which to read .fa files (using the --adr parameter)

## Running

To run Look4TRs on a genome, you'll need two folders consisting of the input directory and output directory.
The input directory should contain .fa files.
Run the following command to obtain microsatellites.

> Look4TRs --adr Input --out Output --default


In order to receive additional help, run one of the following commands

> Look4TRs

> Look4TRs --help

> Look4TRs -h



## Parameters
### REQUIRED OPTIONS
--adr <string>

    Takes the address of the directory of the sequences to be
    scanned. Each file must be a FASTA file and end with .fa

--out <string>

    Takes the address of the output prediction directory.

NOTE: After supplying the required options, either provide a single default
 option, or provide all of the parameter options!

### DEFAULT OPTIONS
--default        

    Runs the tool with the default parameters.

### PARAMETER OPTIONS
--trn <string>

    Takes the address of the directory of the training files.

--motif-location-file <file.bed>

    This optional will make Look4TRs train using repeats that are specified by a bed file.
    Doing this will mean that minm and maxm will not have to be initialized.
    WARNING: This file should correspond with the file used for training Look4TRs or the file indicated by --motif-fa-file!

--motif-fa-file <file.fa>

    This will determine the fa file used to generate repeats for training.
    If ommited, the training file (either specifed by the --trn flag or the largest one in the --adr directory) will be used.
    If included, a bed file (initialized by --motif-loc-file) must be included with it.



--min <integer>

    The minimum k-mer length analyzed (must be less than or equal
    to the maximum k-mer length).
    This is used by the scoring module.
    Default is 4.

--max <integer>

    The maximum k-mer length analyzed (must be greater than or equal
    to the minimum k-mer length).
    This is used by the scoring module.
    Default is 6.

--int <integer>

    The minimum allowed region size for predictions.
    (Must be greater than 0)
    Default is 20.


--win <integer>

    The number of multiples of --int to use. (must be positive)
    This will allow Look4TRs to auto-callibrate the scoring module.
    Default is 4.


--hmm-states-upper <integer>

    the upper bound for the number of hmm states used
    (must be positive, even, and higher than or equal to --hmm-states-lower)
    Default is 12.

--hmm-states-lower integer

    The lower bound for the number of HMM states used
    (must be postive, even, and lower than or equal to --hmm-states-upper)
    Default is 6.


--seg <integer>

    Determines the size of the training chromosome (must be positive)
    Default is 3,000,000


--mimn <integer>

    The minimum size of the motifs used to generate the training Chromosome.
    (must be at least 1 and less than or equal to maxm)
    Default is 1.

--maxm <interger>

    The maximum size of the motifs used to generate the training Chromosome.
    (must be at least 1 and bigger than or equal to minm)
    Default is 10.


--ord <integer>

    The markov order used to generate the training chromosome.
    (must not be negative)
    Default is 0.


--mtf <0|1>

     Enable ('1') or disable ('0') the motif discovery feature.
     Default is 1.

 --prn <integer>

     The number of times that the motif analyzer will train on motifs. Must be positive.
     Default is 2000.


--idn <integer>

    The minimum allowed identity score between the exact repeat and a candidate (Must be 0 or greater).
    Default is 50.




--save-data <directory>

    A directory where the -hmm, --glm, and --chmp files will be stored.

--hmm <file>

    The hmm file that is produced from a previous run that used --save-data.
    This HMM will be loaded, instead of training one.
    Note: If --hmm, --glm, or --chmp is used, all of them must be used!

--glm <file>

    The glm file that is produced from a previous run that used --save-data.
    This GLM will be loaded, instead of training one.
    Note: If --hmm, --glm, or --chmp is used, all of them must be used!

--chmp <file>

    The file that is produced from a previous run that used --save-data.
    This will be used for parameters used from training the HMM and GLM.
    Note: If --hmm, --glm, or --chmp is used, all of them must be used!


--get-testing-data <directory>

    The directory where the training and testing chromosomes will be stored.
    These chromosomes are the ones used to auto-callibrate Look4TRs.
    This will also store the locations of the tandem repeats


--hmm-states <integer>

    The numeber of HMM states (must be positive and even)


--thr <integer>

    Number of threads.
    Default is the total number of cores detect on the user's computer.
