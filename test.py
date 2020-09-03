import time
from random import choice

mychars = ['A', 'C', 'G', 'T']

def kmer2path1(kmer):
    
    path = []
    
    for i in range(len(kmer)):
        
        if kmer[i] in mychars:
            path.append(mychars.index(kmer[i]))
        else:
            return -1
    
    return path

def kmer2path2(kmer):
    
    path = []
    
    for letter in kmer:
        
        try:
            path.append(mychars.index(letter))
        except ValueError:
            return -1
        
    return path

def kmerConstructor(length):

    kmer=""
    for count in range(length):
      kmer+=choice("CGTA")
    
    return kmer



if __name__ == "__main__":

    N = 20
    kmers = [kmerConstructor(40) for _ in range(5000000)]

    #   FIRST METHOD
    start = time.time()
    for kmer in kmers:
        path = kmer2path1(kmer)
    end = time.time()
    print("Method 1")
    print("Completed in  " + str(time.strftime('%H:%M:%S', time.gmtime(end-start))))
    print("Completed in  " + str(end - start) + " seconds")

    #   SECOND METHOD
    start = time.time()
    for kmer in kmers:
        path = kmer2path2(kmer)
    end = time.time()
    print("Method 1")
    print("Completed in  " + str(time.strftime('%H:%M:%S', time.gmtime(end-start))))
    print("Completed in  " + str(end - start) + " seconds")
    