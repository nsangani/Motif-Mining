# Motif-Mining
Script to preprocess CLIP-seq experiment from ENCODE

# Download the script from ENCODE
Select experiment of interest and downlaod the script : https://www.encodeproject.org/matrix/?type=Experiment
Follow the instructions in the script and download the data 

# Preprocessing 
Remove chr GL, KI, X, Y
Split the data into 80/20

# Script Details
clean_ENCODE_functions.py has all the functions needed to clean the files
seq_cleaning.py has code to generate scripts for bedtools among others to align your data with ENCODE
