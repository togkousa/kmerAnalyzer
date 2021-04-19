# A computational framework for pattern detection 1 on unaligned sequences: An application on SARS-CoV-2 data
 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
An alignment-free method capable of processing and counting k-mers in a reasonable time, while evaluating multiple values of the k parameter concurrently.

## Abstract
The rapid development of Bioinformatics after the first half of the 20th century has proved paramount for the decoding of biological sequences. Mining information from the chaos of biological data and separating species in a genetic basis has been widely studied, but nevertheless, still proves to be an extremely difficult problem to consider. Multiple algorithmic processes and techniques have been developed, in order to approach the problem multidimensionally. The current Thesis focuses on the genetic separation of organisms using metagenomic samples, based on the information extracted from DNA sequences. More specifically, we focus on the distribution of k-mers and the information that can be extracted from them. To that end, we have developed an algorithm capable of detecting characteristic k-mers of the sample for many different k-values in reasonable time intervals. We use these k-mers as characteristics for the DNA sequences, and the data are grouped using machine learning algorithms. The developed algorithm is part of an innovative and much promising approach both to the problem of separating organisms using metagenomic data, as well as for the study of changes in the distributions of k-mers, as the k-value is fluctuating within a range of values.


## Installation

### Python Version
The code was developed in Python 2.7 version.

### Required Python Packages
- [`pandas`](https://pandas.pydata.org/getting_started.html) 
- [`numpy`](https://numpy.org/install/)

## Usage

### Input files
The current application supports only `fasta` files as input files.

### How to execute
1. In order to execute the application, there must be a unique `fasta` file inside the `data/` folder, which will be used as an input to the current k-mer analyzer toolkit.
2. Folder `Output/` needS to be empty. Otherwise, the application will remove everything (file or subfolder) inside it. In case the folder doesn't exist, it wil be created automatically.
3. Specify the parameters inside `featuresExtraction.py` script in **lines 21-22**, `kmax` and `eval_factor`. `Eval_factor` parameter determines the strictness in the assessment of kmers of each length. Recommended values for `eval_factor` lie inside the interval [1,2]. For optimal results, it's highly recommended to select a value between [1.2, 1.5]. At any case, for values lower than 1, the application won't run properly
4. Execute the python script `featuresExtraction.py` 

### Output files
Assuming that the input file is called `filename.fasta`:

1. Inside the `Output/` directory there is a `csv` file called `clustData.csv` which is actually the data matrix that we aimed for. Every sequemce is being represented by a number k-mer based features. The value of every feature is the number of times each k-mer was detected in the current sequence.
2. Inside the `Output/filename/` sub-folder there are 3 `csv` files: 
   * File `output.csv` contains the list that is generated from the kmer-tree. Every row represents a k-mer. The first column is the k-mer itself, the second is the length of the k-mer, the third column its frequency (the number of times that was detected in the input data) and the fourth one is its evaluation in the tree. 
   * The two remaining files are associated with the sequences that every k-mer appears, as well as the number of times that each k-mer appears in every sequence occurs.

## Data availability 

SARS-CoV-2 data have been downloaded from [NCBI SARS-CoV-2 Resources](https://www.ncbi.nlm.nih.gov/sars-cov-2/).

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the [MIT](https://opensource.org/licenses/MIT) License - see the [LICENSE](LICENSE) file for details

