# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:52:59 2022

@author: Sangani
"""

# Motif_filename unknown variable needs to be defined prior.

def remove_xy(df, motif_filename):
    '''
    Takes a df and removes chrX,Y,MT, KI, GL + drops column - score, strand, signalValue, qValue, peak
    Output: Curated_filenmae.bed
    '''
    # Remove [chr#_GL, chr#_KI, X, Y] from raw file
    import numpy as np
    chromosomes = np.arange(1,23)
    chr_Un = chromosomes.tolist()
    chr_Un.append('M')
    chr_Un.append('Un')
    for chr_num in chr_Un:
        rm_word = 'chr'+str(chr_num)+'_'
        df = df[~df['chr'].astype(str).str.contains(rm_word)]
       
    chr_sex = ['M','X','Y']
    for chr_num in chr_sex:
        rm_word = 'chr'+str(chr_num)
        df = df[~df['chr'].astype(str).str.contains(rm_word)]
    
    # Clean the motif names from suffix '_K562_Rep01/2' 
    df['motifname'] = df['motifname'].str.slice(0,-11)
    
    # Drop unwanted columns
    df.drop(['score','strand','signalValue','qValue','peak'], axis=1, inplace = True)
    
    
    curated_file = 'curated_'+motif_filename
    df.to_csv(curated_file, header=None, index=None) # cleaned file from chr X,Y,MT
    return print('Sex_chr removed and curated file saved')

def ENCODE_split_80(df, motif_filename):
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
        df.to_csv(final80, header=None, index=None, sep = '\t')
    return print(final80 + ' Done.')

def ENCODE_split_20(df, motif_filename):
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
        new_df = df[df['chr'].str.contains(str(keep_word))]
        final20 = 'final20_'+motif_filename
        new_df.to_csv(final20, header=None, index=None, sep = '\t', mode = 'a')
    return print(final20 + ' Done.')

def getfasta(finalx, motif_filename):
    """
    Input file: final80 or final20 and motif_filename
    Output: Script for bedtools
    """
    seqfile = 'bedtools getfasta -fi /N/project/NGS-JangaLab/Neel/Genome/hg38.fa -bed '+ finalx+'_'+motif_filename+ ' -tab -fo seq/seq_' +finalx +'_'+ motif_filename
    updated_seqfile = r"sed 's/:/\t/g' seq/seq_" +finalx +'_'+ motif_filename + r" | sed 's/-/\t/g' > seq/tab_" + finalx + '_'+ motif_filename

    myfile = open('getfasta.sh', 'a+')
    myfile.writelines("{}\n".format(seqfile))
    myfile.writelines("{}\n".format(updated_seqfile))
    myfile.close()
    return 'Fasta file script retrieved!'

def intersect(finalx, motif_filename):
    file_interesect = 'bedtools intersect -wb -a '+ finalx +'_'+motif_filename +' -b seq/tab_'+finalx+'_'+motif_filename+ r" | awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$5}' > clean_" + finalx + '_' + motif_filename
    myfile = open('intersect.sh', 'a+')
    myfile.writelines("{}\n".format(file_interesect))
    myfile.close()
    return 'Added to intersect file'


def line_prepender(filename, line):
    """
    Input: filename and line added to the top of the file
    Output: Adds line to the top of the script
    """
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
        
def lower_case(filename, newfile):
    '''
    Input: file with seq column
    Ouput: lower cased column to the newfile
    '''
    import pandas as pd
    colname = ['chr', 'start', 'end', 'motif_name', 'seq','pValue']
    seq_motif = pd.read_csv(filename, sep='\t', header = None, names = colname)
    seq_motif['seq'] = seq_motif['seq'].str.lower()
    
    seq_motif.to_csv(newfile, header=None, index=None, sep='\t')
    return 'File converted to lower case'