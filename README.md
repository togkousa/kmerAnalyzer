 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Diploma-Thesis-Togkousidis-8920
This is the full code that was developed during my Diploma Thesis.

# Title
Design and Implementation of a k-mer based feature selection method for machine learning applications in bioinformatics

# Abstract
The rapid development of Bioinformatics after the first half of the 20th century has proved paramount for the decoding of biological sequences. Mining information from the chaos of biological data and separating species in a genetic basis has been widely studied, but nevertheless, still proves to be an extremely difficult problem to consider. Multiple algorithmic processes and techniques have been developed, in order to approach the problem multidimensionally. The current Thesis focuses on the genetic separation of organisms using metagenomic samples, based on the information extracted from DNA sequences. More specifically, we focus on the distribution of k-mers and the information that can be extracted from them. To that end, we have developed an algorithm capable of detecting characteristic k-mers of the sample for many different k-values in reasonable time intervals. We use these k-mers as characteristics for the DNA sequences, and the data are grouped using machine learning algorithms. The developed algorithm is part of an innovative and much promising approach both to the problem of separating organisms using metagenomic data, as well as for the study of changes in the distributions of k-mers, as the k-value is fluctuating within a range of values.

# Python Version
The code was developed in Python 2.7 version so you may need to intsall it before you run the application.

# Input files
The current application supports only .fasta files as input files.

# How to execute
1. In order to run the application, there must be 
2. Folders 'Input/', 'Output/' and 'ClusteringData/' need to be empty. Otherwise, the application will remove everything (file or subfolder) inside them.
3. Place the input .fasta file into folder 'data/'. (A test file is uploaded in the current github reposiroty called 'data/sample_dnaseq1.fasta')
4. Go to terminal and change directory to your current local repository of the application.
5. Execute the python script "featuresExtraction.py" -> (type the following command: $ python featuresExtraction.py)

Steps 3 and 4 can also be executed by using a python IDE.

# Outputs
Assuming that the input file is called 'filename.fasta':

1. Inside the 'Input/' folder there are 3 files. The first one is a .txt file called 'filename.txt' which contains the DNA sequences in a shuffled order. The second one is a .csv file called 'filename_sequenceIDs.csv' that contains the IDs of the sequences in the same shuffled order. The third one is also a .csv file called 'filename_sequencesIDs_unshuffled.csv' which contains the mapping of the sequences to the corresponding IDs in the original unshuffled order of the .fasta file.
2. Inside the 'Output/filename/' folder there are 3 .csv files. File 'output.csv' contains the list that is generated from the kmer-tree. Every row represents a k-mer. The first column is the k-mer itself, the second is the length of the k-mer, the third column its frequency (the number of times that was detected in the input data) and the fourth one is its evaluation in the tree. The two remaining files are associated with the sequences that every k-mer appears, as well as the number of times that each k-mer appears in every sequenceoccurs 
3. Inside the 'ClusteringData/' directory there is a .csv file called 'clustData.csv' which is actually the data matrix that we aimed for. Every sequemce is being represented by a number k-mer based features. The value of every feature is the number of times each k-mer was detected in the current sequence.

# Example
1. Clone the GitHub repository into a local directory.
2. Change directory to your current local repository.
3. Inside the 'data/' folder there is a small .fasta file called 'sample_dnaseq1.fasta'. The application will be executed with this file as an input.
4. Execute the python script 'featuresExtraction.py'
5. Output files should be pretty similar to those uploaded to the current repository. Small differences in the output files may be detected because of the shuffling process. However, these differences will be of minor importance, as the most significant k-mers will be accurately detected.


