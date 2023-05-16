#!/bin/bash

import pandas as pd
import clean_ENCODE_functions as f

motifs_col = ['a','b','c','d','e','f']
motifs = pd.read_csv('ENCODE_eclipData_script.txt', sep = '/', dtype=(str), header = None, names = motifs_col)
filename = motifs['f'].tolist()

n = 0
fasta = 0
for motif_name in filename:
    motif_filename = motif_name[:-7]+'.bed'
    print('Loading ' + motif_name)
    colname = ['chr', 'start', 'end', 'motifname', 'score','strand','signalValue', 'pValue', 'qValue', 'peak']
    seq_motif = pd.read_csv(motif_filename, sep='\t', header = None, names = colname)
    
    f.remove_xy(seq_motif, motif_filename)
    curated_file = 'curated_'+motif_filename
    curated_file = pd.read_csv(curated_file, names = ['chr', 'start', 'end', 'motifname', 'pValue'])
    
    f.ENCODE_split_80(curated_file, motif_filename)
    f.ENCODE_split_20(curated_file, motif_filename)
    n += 1
    print('Total ENCODE_motif_files done: ' + str(n))
    
    f.getfasta('final80', motif_filename)
    f.getfasta('final20', motif_filename)
    fasta+= 1
    print('Total fasta file script appended: '+ str(fasta) + '\n')
    
    f.intersect('final80', motif_filename)
    f.intersect('final20', motif_filename)
    
f.line_prepender('getfasta.sh', '#!/bin/bash')    
f.line_prepender('intersect.sh', '#!/bin/bash') 
