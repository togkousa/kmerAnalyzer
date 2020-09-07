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
1. Folders '/Input', '/Output' and '/ClusteringData' need to be empty. Otherwise, the application will remove everything (file or subfolder) inside them.
2. Place the input .fasta file into folder '~/data'. (A test file is uploaded in the current github reposiroty called '/data/sample_dnaseq1.fasta')
3. Go to terminal and change directory to your current local repository of the application.
4. Execute the python script "featuresExtraction.py" -> (type the following command: $ python featuresExtraction.py)

Steps 3 and 4 can also be executed by using a python IDE.

# Outputs
Assuming that the input file is called 'filename.fasta'!
1.Inside the '/Input' folder there are 3 files. The first one is a .txt file called 'filename.txt' which contains the DNA sequences in a shuffled order. The second one is a .csv file called 'filename_sequenceIDs.csv' that contains the IDs of the sequences in the same shuffled order. The third one is also a .csv file called 'filename_sequencesIDs_unshuffled.csv' which contains the mapping of the sequences to the corresponding IDs in the original unshuffled order of the .fasta file.
2. Inside the 'Output/filename/' folder there are 3 .csv files. File 'output.csv' contains 
3.

# Example



