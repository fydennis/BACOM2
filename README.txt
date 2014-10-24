 
   Table of Contents
   -----------------

   1. Enrivonment

   2. Input

   3. Output

   4. Contact Information
   
   
1. Enrivonment 

The BACOM2 program requires Matlab 2013a(64bit) or higher version. It has been successfully tested under Windows 7 enrivonment. RAM with 4G or higher volume is recommended.

The main program named "BACOM2.m" could be opened, edited and run by Matlab. 

2. Input

For current matlab function does not support direct annoated reading on CEL files from Affymetrix SNP6 array, the input needs to be extracted from the intermediate results from original BACOM (the code could be downloaded from https://code.google.com/p/bacom or https://code.google.com/p/aisaic ). For each pair of tumor-normal matched samples, the required input include two files from original BACOM's output files: one is a CSV file that starts with the tumor sample CEL file name and ends with "_outputdata.csv", in which records probeset IDs, their located chromosomes and positions, and intensity signals for both alleles of the two samples; the other is a CSV file that starts with the same CEL file name and ends with "_outputdetectionResult.csv", in which records the segmentation information. These files should be placed under the sub-directory named "Samples".
Another text file named "SampleList.txt" is also needed in the same directory of the program. It should list all the samples' IDs that are going to be analyized, one ID in a line.

3. Output
The program will generate a text file named "BACOM2_Results.txt" under the same directory. In each line, there will be a sample ID compined with its esitmated tumor purity and average ploidy.

4. Contact Information
For any question about the program, please contact yifu@vt.edu