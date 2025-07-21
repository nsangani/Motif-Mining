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
#######################################################################
from gensim.models import Word2Vec  #for getting kmer embedding

##################################################
# training a DescisionTreeClassifier 
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
##########clf = GaussianNB()
#clf=LogisticRegression()  #LogisticRegression(C=1e4)
#clf=LogisticRegression(random_state=random.seed(1234))
#clf=DecisionTreeClassifier(max_depth = 2)
#clf=KNeighborsClassifier(n_neighbors=2)#(n_neighbors=3)default#n_neighbors=5
#clf=svm.SVC()
clf=RandomForestClassifier() #RandomForestClassifier(n_estimators=30, max_depth=10, random_state = random.seed(1234))#random_state=0)

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

df=df=pd.read_csv('multi_class_popseq_k562_20_balanced.csv')
#df=dataset
df['chr']=df['chr'].astype('category')
df['chr'] = df['chr'].cat.codes
print(df['chr'].head())
features = ['chr','p_value'] 
print(df[features].head())
X = df[features]
print("@@@@@@@@@@@@@")
print(df['seq'].head())
print("^^^^^^^^^^^")
#insert onehot encoding of reference-kmer
#Onehot=pd.get_dummies(df['seq'], prefix='seq')
#X= pd.concat([X,Onehot],axis=1)
#X=Onehot

print("###############")
dff=df['RBP'].to_frame()# to convert series to df
print("###############")
print(type(dff))
#dff = dff.select_dtypes(exclude=[np.number])
print(dff.head())
#print("###############",dff.unique())
#dff.value_counts().to_csv('unique_RBPs.csv')
dff_series = dff.squeeze() # to convert df to series
print("$$$$$$$",dff_series.value_counts())
print(dff_series.dtype)
print("###############")
#convert from object to ctegory type
dff_series=dff_series.astype('category')
dff_series = dff_series.cat.codes
print(dff_series.head())
print("###############")

y = dff_series



print("#############",X.shape)
print("=============",y.shape)


##########################################################
#needed for sequence_embedding
seq_list_train=list(df['seq']) 
print(len(seq_list_train))
#########################################
#create the embedding model
processed_corpus=seq_list_train
print("&&&&&&&&",len(processed_corpus))
processed_corpus=[processed_corpus]
print("#########",len(processed_corpus[0]))
model = Word2Vec(processed_corpus,vector_size=20,min_count=1,window=3)
#model = Word2Vec(processed_corpus,size=100,min_count=1)#, window=5) #min_count=5 leads to an error
print(model)
#words = list(model.wv.vocab)#old version of genism see https://github.com/RaRe-Technologies/gensim/wiki/Migrating-from-Gensim-3.x-to-4
#words = list(model.wv) #raise error at https://stackoverflow.com/questions/72480289/how-to-handle-keyerrorfkey-key-not-present-wor2vec-with-gensim
#print(words)
##############################
#get kmer_embedding of train dataset
#################################################
df3=pd.DataFrame()
for i in range(len(seq_list_train)):
    x=seq_list_train[i]
    #print(type(x))
    ###print(model[x])
    df4=pd.DataFrame(model.wv[x]) #use model[x] in old version of gensim <4.00
    df4=df4.T
    # to append df4 at the end of df3 dataframe 
    df3 = pd.concat([df3,df4])
df3_train=df3
print(df3_train.head())
print("8888888888",df3_train.shape)    
######insert embedding of reference-kmer 
df3_train.reset_index(drop=True, inplace=True)      #To avoid the error at https://vispud.blogspot.com/2018/12/pandas-concat-valueerror-shape-of.html
X.reset_index(drop=True, inplace=True)
X= pd.concat([X,df3_train],axis=1)
#X=df3_train  #TO test the effect of embedding only
##############################

#print(X.head())

#scale training data
#X= scaler.fit_transform(X)
print(",,,,,,,,",X.shape)


############################################################################################
#train, test = train_test_split(df, test_size=0.2)   

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

print(clf)
print('type of y_train=', type(y_train))
print('type of y_test=', type(y_test))

model = clf.fit(X_train, y_train) 
#model = clf.fit(X_train, y_train.ravel()) 
y_pred = model.predict(X_test)
print('type of y_pred=', type(y_pred))
###############################################
#y_prob = clf.predict_proba(X_test.values.reshape(-1, 1))
#y_prob = y_prob[:,1]
# Evaluate the model: Model Accuracy, how often is the classifier correct
from sklearn import metrics #Import scikit-learn metrics module for accuracy calculation
from sklearn.metrics import classification_report #for classifier evaluation
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score # for printing AUC
from sklearn.metrics import confusion_matrix
# creating a confusion matrix 
cm = confusion_matrix(y_test, y_pred )
print(cm)

print("Accuracy:",metrics.accuracy_score(y_test, y_pred)*100)
 
print(classification_report(y_test, y_pred))
#auc=roc_auc_score(y_test.round(),y_pred,multi_class="ovr",average=None)
#auc = float("{0:.3f}".format(auc))
#print("AUC=",auc)

#true negatives c00, false negatives C10, true positives C11, and false positives C01 
#tn c00, fpC01, fnC10, tpC11 
print('CF=',confusion_matrix(y_test, y_pred))
l=confusion_matrix(y_test, y_pred)#https://towardsdatascience.com/accuracy-precision-recall-or-f1-331fb37c5cb9
print('TN',l.item((0, 0)))
print('FP',l.item((0, 1)))
print('FN',l.item((1, 0)))
print('TP',l.item((1, 1)))
#print(type(X_train), type(y_train))

print("@@@@@@@@@@@", clf)




