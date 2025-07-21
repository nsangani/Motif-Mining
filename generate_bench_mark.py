# To set seed random number in order to reproducable results in keras
from numpy.random import seed
seed(4)
#import tensorflow
#tensorflow.random.set_seed(1234)
########################################
import pandas as pd
from pandas import *
import numpy as np
import random
import pickle
import joblib
##################################################
# training a DescisionTreeClassifier 
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
##########clf = GaussianNB()
clf=LogisticRegression()  #LogisticRegression(C=1e4)
#clf=LogisticRegression(random_state=random.seed(1234))
#clf=DecisionTreeClassifier(max_depth = 2)
#clf=KNeighborsClassifier(n_neighbors=2)#(n_neighbors=3)default#n_neighbors=5
#clf=svm.SVC()
#clf=RandomForestClassifier() #RandomForestClassifier(n_estimators=30, max_depth=10, random_state = random.seed(1234))#random_state=0)

#classifier =svm.SVC(gamma='scale',C=1,probability=True)
####################################################
import plot_learning_curves as plc
from sklearn.preprocessing import MinMaxScaler #For feature normalization
scaler = MinMaxScaler()

from sklearn.model_selection import train_test_split


# Evaluate the model: Model Accuracy, how often is the classifier correct
from sklearn import metrics #Import scikit-learn metrics module for accuracy calculation
from sklearn.metrics import classification_report #for classifier evaluation
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score # for printing AUC
from sklearn.metrics import confusion_matrix
from sklearn.utils import shuffle

#read labeled encode data
encode = pd.read_csv("K562_ENCODE_20_clean.bed.txt",sep='\t',header=None)

print(encode.head())
print(encode.info())
print("encode.shape=",encode.shape)
encode.columns=['chr', 'start', 'end', 'RBP', 'seq', 'p_value']

print(len(list(encode['chr'])))
print("encode",len(list(encode['start'])))
print("encode",len(list(encode['end'])))
print("encode",len(list(encode['RBP'])))
print("encode",len(list(encode['seq'])))
print("uniuqe_encode_seq",len(list(encode['seq'].value_counts())))
print("encode",len(list(encode['p_value'])))

      
popseq=pd.read_csv('K562_NPOP.final20.bed.txt', sep='\t',header=None)    # To predict RBPs on popseq unlabled data K562_popseq_seq_pval.txt
  #subset=['start']

print(popseq.head())
print(popseq.info())
print("shape of popseq=",popseq.shape)

popseq.columns=['chr', 'start', 'end', 'seq', 'p_value']
print("popseq",len(list(popseq['chr'])))
print("popseq",len(list(popseq['start'])))
print("popseq",len(list(popseq['end'])))
print("popseq",len(list(popseq['seq'])))
print("uniuqe_popseq_seq",len(list(popseq['seq'].value_counts())))
print("popseq",len(list(popseq['p_value'])))


#label the data
start_list=list(set(popseq['start']).intersection(set(encode['start'])))

print("length of start intersection list",len(start_list))

end_list=list(set(encode['end']).intersection(set(popseq['end'])))
print("length of end intersection list@@@",len(end_list))



encode1=encode[encode['start'].isin(start_list)]
print("encode1.shape=",encode1.shape)
encode1 = encode1.drop_duplicates(subset=['start'])
print("encode1.shape_no_duplicates=",encode1.shape)

print(encode1.head())
#print("/****encode1_Shape=",encode1.shape)

popseq_1=popseq[popseq['start'].isin(start_list)]
popseq_1 = popseq_1.drop_duplicates(subset=['start'])
print("shape of popseq_1_no_duplicates=",popseq_1.shape)
print(popseq_1.head())
print("/****popseq_1_Shape=",popseq_1.shape)


popseq_01=popseq[~popseq['start'].isin(start_list)]

popseq_01 = popseq_01.drop_duplicates(subset=['start'])
print("shape of popseq_01_no_duplicates=",popseq_01.shape)
print(popseq_01.head())
print("/****popseq_01_Shape=",popseq_01.shape)


new_df1 =encode1.combine_first(popseq_1)
print("new_df1_shape=",new_df1.shape)
print("new_df1_columns=",new_df1.columns)
new_df1 = new_df1.dropna()
popseq_11=new_df1[new_df1['start'].isin(list(set(popseq['start'])))]
print(":::::",popseq_11.columns)
print("###/*****/popseq_11_Shape=",popseq_11.shape)



end_list=list(set(encode['end']).intersection(set(popseq['end'])))
print("length of end intersection list@@@",len(end_list))

encode2=encode[encode['end'].isin(end_list)]
print(encode2.head())
print("encode2_Shape=",encode2.shape)
encode2 = encode2.drop_duplicates(subset=['end'])
print("encode2.shape_no_duplicates=",encode2.shape)


popseq_2=popseq[popseq['end'].isin(end_list)]
print("length of popseq_2", len(popseq_2))
popseq_2 = popseq_2.drop_duplicates(subset=['end'])
print("shape of popseq_2_no_duplicates=",popseq_2.shape)

popseq_02=popseq[~popseq['end'].isin(end_list)]

popseq_02 = popseq_02.drop_duplicates(subset=['end'])
print("shape of popseq_02_no_duplicates=",popseq_02.shape)
print(popseq_02.head())
print("/****popseq_02_Shape=",popseq_02.shape)


print(popseq_2.shape)
new_df2 =encode2.combine_first(popseq_2)
print("new_df2_shape=",new_df2.shape)
print("new_df2_columns=",new_df2.columns)
new_df2 = new_df2.dropna()
print("new_df2_shape=",new_df2.shape)
popseq_22=new_df2[new_df2['end'].isin(list(set(popseq['end'])))]
print(popseq_22.columns)
print("###/*****/popseq_22_Shape=",popseq_22.shape)


new_df3=pd.concat([popseq_11,popseq_22])
print("*****/new_df3_shape=",new_df3.shape)


new_df4=pd.concat([popseq_01,popseq_02])
print("*****/new_df4_shape=",new_df4.shape)
new_df4= new_df4.sample(n=len(new_df3), replace=False) #try replace=false
print("*****/new_df4_shape=",new_df4.shape)

listofzeros=[0]*len(new_df4.index)
new_df4.insert(5,"RBP", listofzeros, True)

# Create DataFrame from multi and zero class examples


dataset = new_df4.append(new_df3, ignore_index=True)  #for multi-classification 
print(dataset.columns)
#multilabel_dataset.to_csv('labled_popseq_k562_b.csv')
dataset.to_csv('multi_class_popseq_k562_20_balanced.csv')




