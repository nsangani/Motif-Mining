# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 23:26:41 2022

@author: Sangani
"""

import pandas as pd
import numpy as np
import os

os.chdir('C:\\Users\\Sangani\\Desktop')
#path = os.getcwd()

motifs_col = ['a','b','c','d','e','f']
motifs = pd.read_csv('ENCODE_eclipData_script.txt', sep = '/', dtype=(str), header = None, names = motifs_col)
motifs.head()

filename = motifs['f'].tolist()
myfile = open('ENCODE_filenames.txt', 'w')
for file in filename:
    file = file[:-7]+'.bed'
    myfile.writelines("%s\r" % file)
myfile.close()

print('Total motif files:', len(filename)) # 99 motifs

# --------------------------------------------------------------------------------------------
# Processing ENCODE files and creating 80/20 split
n_80= 0
n_20 = 0
for motif_name in filename:
    #motif_filename = motif_name[:-7]+'.bed'
    
motif_filename = 'ENCFF005ZCI.bed'
print('Loading ' + motif_name)

# import ENCODE bed file and assign header
colname = ['chr', 'start', 'end', 'motifname', 'score','strand','signalValue', 'pValue', 'qValue', 'peak']
seq_motif = pd.read_csv(motif_filename, sep='\t', header = None, names = colname) #index_col=0
seq_motif.iloc[1: , :]
seq_motif['chr'].unique()

def remove_xy(seq_motif):
    # Remove [chr#_GL, chr#_KI, X, Y] from raw file
    chromosomes = np.arange(1,23)
    chr_Un = chromosomes.tolist()
    chr_Un.append('M')
    chr_Un.append('Un')
    for chr_num in chr_Un:
        rm_word = 'chr'+str(chr_num)+'_'
        seq_motif = seq_motif[~seq_motif['chr'].astype(str).str.contains(rm_word)]
       
    chr_sex = ['M','X','Y']
    for chr_num in chr_sex:
        rm_word = 'chr'+str(chr_num)
        seq_motif = seq_motif[~seq_motif['chr'].astype(str).str.contains(rm_word)]
    
    # Clean the motif names from '_K562_Rep01'
    seq_motif['motifname'] = seq_motif['motifname'].str.slice(0,-11)
    
    # Drop unwanted columns
    seq_motif.drop(['score','strand','signalValue','qValue','peak'], axis=1, inplace = True)
    
    
    curated_file = 'curated_'+motif_filename
    seq_motif.to_csv(curated_file, header=None, index=None) # cleaned file from chr X,Y,MT
    return print('Sex_chr removed and curated file saved')

seq_motif['motifname'].unique()


# Check the curated file 

with open('curated_ENCFF005ZCI.bed', 'r') as file:
    #l = file.readlines()
    head = [next(file) for x in range(5)]
    #print(repr(l))
    print(head)

# import the curated file as df
file = pd.read_csv('curated_ENCFF005ZCI.bed', names = ['chr', 'start', 'end', 'motifname', 'pValue'] )
file.head()

# for line in lines:
#     updated = line.strip("\t")
#     print(updated)
    
#22856679  22856687 chr

# read csv file 
 import csv
 with open(curated_file, 'r') as csv_file:
     csv_reader = csv.DictReader(csv_file)
     
     with open('final80', 'w') as new_file:
         fieldnames = []
        
         csv_writer = csv.DictWriter('final80', fieldname, delimiter = '\t')
         csv_writer.writeheader()
         
         next(csv_reader) #skip the header
         for line in csv_reader():
             csv.writer.writerow(line)
    
# Remove selected Chromosomes for 80/20 split
colname = 
df = pd.read_csv(curated_file, names = ['chr', 'start', 'end', 'motifname', 'pValue']) #index_col=0
df.head()
df['chr'].unique()

def remove_chr(df):
    """
    Input df with a chr column name starting with 'chr' 
    Note: Following chr will be removed: [17,18,19,20,21,22]
    Output chr row removes and updated the df
    """
    # Remove selected Chromosomes for 80 split
    chr_rm = [17,18,19,20,21,22]
    for rm_word in chr_rm:
        rm_word = 'chr'+str(rm_word)
        df = df[~df['chr'].astype(str).str.contains(rm_word)]
        final80 = 'final80_'+motif_filename
        df.to_csv(final80, header=None, index=None)
    return print(final80 + ' Done.')

remove_chr(df) #80% split data for ENCODE
n_80 += 1
print('final80 count: ' + n_80)

# check the saved final80 file
#df = pd.read_csv(final80, names = ['chr', 'start', 'end', 'motifname', 'pValue'])




def removed_chr(df):
    """
    Input df with a chr column name starting with 'chr' 
    Note: Following chr will be kept: [17,18,19,20,21,22]
    Output chr other than the listed will be removed and update the df
    """
    # Keep selected Chromosomes for 20 split
    chr_lst = [17,18,19,20,21,22]
    for keep_word in chr_lst:
        keep_word = 'chr'+ str(keep_word)
        #print(keep_word)
        df = df[df['chr'].str.contains(str(keep_word))]
        final20 = 'final20_'+motif_filename
        df.to_csv(final20, header=None, index=None)
    return print(final20 + ' Done.')

removed_chr(df) #20% data file
n_20 += 1
print('final20 count: ' + n_20)


## some rows have tabs, commas, quotation. Convert all to tab after generating the file
#  awk '{print $0}' final20_ENCFF853FGC.bed | sed 's/"//g' | tail

myfile = open('final80_tab.sh', 'w')
for motif_name in filename:
    motif_filename = motif_name[:-7]
    seqfile = r"sed -i 's/^I*$//' final80_"+motif_filename+'.bed'
    myfile.writelines("%s\r" % seqfile)
    seqfile = r"sed -i 's/^I*$//' final80_"+motif_filename+'.bed'
    myfile.writelines("%s\r" % seqfile)
myfile.close()
    
# ----------------------------------------------------------------------------------------------
# Create .sh file to get sequences from the coordinates & do tab correction

myfile = open('intersect_final80.sh', 'w')
for motif_name in filename:
    motif_filename = motif_name[:-7]
    seqfile = r"sed -i 's/^I*$//' final80_"+motif_filename+'.bed'
    myfile.writelines("%s\r" % seqfile)
    seqfile = 'bedtools getfasta -fi /N/project/NGS-JangaLab/Neel/Genome/hg38.fa -bed '+ 'final80_'+motif_filename+'.bed' + ' -tab -fo seq/'+ motif_filename +'.seq.final80.bed'
    myfile.writelines("%s\r" % seqfile)
    seqfile = r"sed 's/:/\t/g' seq/"+ motif_filename + r".seq.final80.bed | sed 's/-/\t/g' > seq/" + motif_filename + ".tab_final80.bed"
    myfile.writelines("%s\r" % seqfile)
    seqfile = r"bedtools intersect -wb -a final80_" + motif_filename+ ".bed -b seq/" + motif_filename + r".tab_final80.bed | awk -v OFS='\t' '{ print $1,$2,$3,$4,$9,$5}' > " + motif_filename + '.final80.bed'  
    myfile.writelines("%s\r" % seqfile)
myfile.close()

myfile = open('intersect_final20.sh', 'w')
for motif_name in filename:
    motif_filename = motif_name[:-7]
    seqfile = r"sed -i 's/^I*$//' final20_"+motif_filename+'.bed'
    myfile.writelines("%s\r" % seqfile)
    seqfile = 'bedtools getfasta -fi /N/project/NGS-JangaLab/Neel/Genome/hg38.fa -bed '+ 'final20_'+motif_filename+'.bed' + ' -tab -fo seq/'+ motif_filename +'.seq.final20.bed'
    myfile.writelines("%s\r" % seqfile)
    seqfile = r"sed 's/:/\t/g' seq/"+ motif_filename + r".seq.final20.bed | sed 's/-/\t/g' > seq/" + motif_filename + ".tab_final20.bed"
    myfile.writelines("%s\r" % seqfile)
    seqfile = r"bedtools intersect -wb -a final20_" + motif_filename+ ".bed -b seq/" + motif_filename + r".tab_final20.bed | awk -v OFS='\t' '{ print $1,$2,$3,$4,$9,$5}' > " + motif_filename + '.final20.bed'  
    myfile.writelines("%s\r" % seqfile)
myfile.close()


#----------------------------------------------------------------------------------------------
# Upper case to lower case

motifs_col = ['a','b','c','d','e','f']
motifs = pd.read_csv('ENCODE_eclipData_script.txt', sep = '/', dtype=(str), header = None, names = motifs_col)

file = motifs['f'].tolist()

for motif_name in file:
    print('Loading '+ motif_name)
    motif_filename = motif_name[:-7]+'.tab_seq.bed'
    colname = ['chr', 'start', 'end', 'seq','pValue']
    seq_motif = pd.read_csv(motif_filename, sep='\t', header = None, names = colname)
    seq_motif['seq'] = seq_motif['seq'].str.lower()
    
    filename = motif_name[:-7] +'.tab_final.bed'
    print('Writing '+ filename + '\n')
    seq_motif.to_csv(filename, header=None, index=None, sep='\t', mode='a')

#----------------------------------------------------------------------------------------------
# Create .sh to merge '*.tab_final' and '*_updated.bed' to combine Sequence and pValue info.

myfile = open('merge_seq_info.sh', 'w')
for motif_name in file:
    motif_filename = motif_name[:-7]
                    
    myfile.writelines("%s\r" % seqfile)
myfile.close()



for i in file:
    #seqfile = i[:-7]+'.seq.bed'
    seqfile = 'ENCFF005ZCI.seq.bed'
    seq_motif = pd.read_csv(seqfile, sep='\t')
    seq_motif.columns = ['a','b']
    remove_chr = ['chr11_','chr14_','chr15_','chr1_','chr22_','chrUn_']
    for rm in remove_chr:
        seq_motif = seq_motif[~seq_motif['a'].astype(str).str.contains(rm)]
    seq_motif['b'] = seq_motif['b'].str.lower()
    # seq_motif.head()
    seq_motif[0]
    seq_motif = seq_motif.rename(columns=seq_motif.iloc[0]).drop(seq_motif.index[0])
    seq_motif.to_csv(seqfile, sep = '\t', index = False)
    #n += 1 
    #print(n)

# convert : and - to tab
myfile = open('tab_correction.sh', 'w')
for i in file:
    seqfile = r"sed 's/:/\t/g' " + i[:-7]+ r".seq.bed | sed 's/-/\t/g' > " + i[:-7]+".tab_seq.bed"
    myfile.writelines("%s\r" % seqfile)
myfile.close()
    
#bedtools intersect for final files
#bedtools intersect -wb -a ENCFF005ZCI.bed -b seq/ENCFF005ZCI.tab_seq.bed | awk -v OFS='\t' '{print $1,$14,$4,$7,$8}'
myfile = open('intersect.sh', 'w')
for i in file:
    seqfile = r"bedtools intersect -wb -a " + i[:-7]+ r".bed -b seq/" + i[:-7]+ r".tab_seq.bed | awk -v OFS='\t' '{print $1,$2,$3,$14,$4,$7,$8}' > " + \
    i[:-7]+ r".final.bed"
    myfile.writelines("%s\r" % seqfile)
myfile.close()

# ---------------------------------------------------------------------------------------------------------
# NPOP K562 filter
import pandas as pd
import numpy as np
import os

os.chdir('C:\\Users\\Sangani\\Desktop\\Research\\K562 Motif Mining\\K562_NPOP')

K562_col = ['chr','start','end','id','score','strand','pValue']
df = pd.read_csv('K562_NPOP_merged_final.bed', sep = '\t', dtype=(str), header = None, names = K562_col)

def remove_sex_chr(df):
    """
    Input df with a chr column name starting with 'chr'
    Note following chr will be removed: ['KI','GL','MT','X','Y']
    Output removes KI and GL chr rows and update the df
    """
    chr_extra = ['KI','GL','MT','X','Y']
    for rm_word in chr_extra:
        rm_word = 'chr'+ rm_word
        df = df[~df['chr'].astype(str).str.contains(rm_word)]
        df.to_csv('curated_K562_file.bed', sep = '\t', index = False, header = None)
    return df['chr'].unique()

remove_sex_chr(df)

df = pd.read_csv('curated_K562_file.bed', sep = '\t', dtype=(str), header = None, names = K562_col)

def remove_chr(df):
    """
    Input df with a chr column name starting with 'chr' 
    Note: Following chr will be removed: [17,18,19,20,21,22]
    Output chr row removes and updated the df
    """
    # Remove selected Chromosomes for 80 split
    chr_rm = [17,18,19,20,21,22]
    for rm_word in chr_rm:
        rm_word = 'chr'+str(rm_word)
        df = df[~df['chr'].astype(str).str.contains(rm_word)]
    df.to_csv('final80_cleaned_K562_NPOP.bed', sep = '\t', index = False, header = None)
    return df['chr'].unique()

remove_chr(df)

df = pd.read_csv('curated_K562_file.bed', sep = '\t', dtype=(str), header = None, names = K562_col)
def removed_chr(df):
    """
    Input df with a chr column name starting with 'chr' 
    Note: Following chr will be kept: [17,18,19,20,21,22]
    Output chr other than the listed will be removed and update the df
    """
    # Keep selected Chromosomes for 20 split
    chr_lst = [17,18,19,20,21,22]
    for keep_word in chr_lst:
        keep_word = 'chr'+ str(keep_word)
        print(keep_word)
        df1 = df[df['chr'].str.contains(keep_word)]
        df1.to_csv('final20_cleaned_K562_NPOP.bed', mode = 'a', sep = '\t', index = False,header = None)
    return df['chr'].unique()

removed_chr(df)

df['chr'].unique()

