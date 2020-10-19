# A computational framework for pattern detection 1 on unaligned sequences: An application on SARS-CoV-2 data
 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
An alignment-free method capable of processing and counting k-mers in a reasonable time, while evaluating multiple values of the k parameter concurrently.

## Installation

### Python Version
The code was developed in Python 2.7 version.

## Usage

### Input files
The current application supports only `fasta` files as input files.

### How to execute
1. In order to execute the application, there must be a unique `fasta` file inside the `data/` folder, which will be used as an input to the current k-mer analyzer toolkit.
2. Folders `Input/`, `Output/` and `ClusteringData/` need to be empty. Otherwise, the application will remove everything (file or subfolder) inside them.
3. Execute the python script `featuresExtraction.py` 

### Output files
Assuming that the input file is called `filename.fasta`:

1. Inside the `Input/` folder there are 3 files:
  * The first one is a `txt` file called `filename.txt` which contains the DNA sequences in a shuffled order. 
  * The second one is a `csv` file called `filename_sequenceIDs.csv` that contains the IDs of the sequences in the same shuffled order. 
  * The third one is also a `csv` file called `filename_sequencesIDs_unshuffled.csv` which contains the mapping of the sequences to the corresponding IDs in the original unshuffled order of the `fasta` file.
2. Inside the `Output/filename/` folder there are 3 `csv` files: 
   * File `output.csv` contains the list that is generated from the kmer-tree. Every row represents a k-mer. The first column is the k-mer itself, the second is the length of the k-mer, the third column its frequency (the number of times that was detected in the input data) and the fourth one is its evaluation in the tree. 
   * The two remaining files are associated with the sequences that every k-mer appears, as well as the number of times that each k-mer appears in every sequence occurs.
3. Inside the `ClusteringData/` directory there is a `csv` file called `clustData.csv` which is actually the data matrix that we aimed for. Every sequemce is being represented by a number k-mer based features. The value of every feature is the number of times each k-mer was detected in the current sequence.

## Data availability 

SARS-CoV-2 data have been downloaded from [NCBI SARS-CoV-2 Resources](https://www.ncbi.nlm.nih.gov/sars-cov-2/).

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the [MIT](https://opensource.org/licenses/MIT) License - see the [LICENSE](LICENSE) file for details
