# A computational framework for pattern detection on unaligned sequences: An application on SARS-CoV-2 data
 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
An alignment-free method capable of processing and counting k-mers in a reasonable time, while evaluating multiple values of the k parameter concurrently.

## Installation

### Python Version
The code was developed in Python 2.7 version.

### Required Python Packages
- [`pandas`](https://pandas.pydata.org/getting_started.html) 
- [`numpy`](https://numpy.org/install/)

## Usage

### Input files
The current application supports only `.fasta` files as input files.

### How to execute
1. In order to execute the application, there must be a unique `fasta` file inside the `data/` folder, which will be used as an input to the current k-mer analyzer toolkit.
2. Folder `Output/` needs to be empty. Otherwise, the application will remove everything (file or subfolder) inside it. In case the folder doesn't exist, it wil be created automatically.
3. Specify the parameters inside `featuresExtraction.py` script in **lines 22-23**, `kmax` and `eval_factor`. `Eval_factor` parameter determines the strictness in the assessment of kmers of each length. Recommended values for `eval_factor` lie inside the interval [1,2]. For optimal results, it's highly recommended to select a value between [1.2, 1.5]. At any case, for values lower than 1, the application won't run properly
4. Execute the python script `featuresExtraction.py` 

### Output files
Assuming that the input file is called `filename.fasta`:

1. Inside the `Output/` directory there is a `.csv` file called `clustData.csv` which is actually the data matrix that we aimed for. Every sequemce is being represented by a number k-mer based features. The value of every feature is the number of times each k-mer was detected in the current sequence. There's also a CSV file called `heades_to_IDS.csv`, which maps the headers of each sequence from fasta input to code names ID-1, ID-2 etc.
2. Inside the `Output/filename/` sub-folder there are 3 `csv` files: 
   * File `output.csv` contains the list that is generated from the kmer-tree. Every row represents a k-mer. The first column is the k-mer itself, the second is the length of the k-mer, the third column its frequency (the number of times that was detected in the input data) and the fourth one is its evaluation in the tree. 
   * The two remaining files are associated with the sequences that every k-mer appears, and the number of times that each k-mer appears in every sequence occurs as well.

### Extra comments
- It's important to have a look at the lengths of the sequences, prior to executing kmerAnalyzer. For example, in the example dataset, the sequence with header `>ERR525627.984.1 984 length=31` has length 31, so probably its better to either exclude this sequence (data filtering) or examine lower k-values, e.g. up to 20. However, if we set `kmax = 35`, the code seems to work properly.
- While executing kmerAnalyzer, a folder called `input` is created isnide the project direcoty, containing some necessary files for the execution process. The folder is deleted at the end of the process.

## Data availability 

SARS-CoV-2 data have been downloaded from [NCBI SARS-CoV-2 Resources](https://www.ncbi.nlm.nih.gov/sars-cov-2/).

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the [MIT](https://opensource.org/licenses/MIT) License - see the [LICENSE](LICENSE) file for details
