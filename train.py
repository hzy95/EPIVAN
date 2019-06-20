

# In[ ]:
import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
from model import get_model
import numpy as np
from keras.callbacks import Callback
from datetime import datetime
from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.model_selection import train_test_split

class roc_callback(Callback):
    def __init__(self, val_data,name):
        self.en = val_data[0]
        self.pr = val_data[1]
        self.y = val_data[2]
        self.name = name

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        y_pred = self.model.predict([self.en,self.pr])
        auc_val = roc_auc_score(self.y, y_pred)
        aupr_val = average_precision_score(self.y, y_pred)

        self.model.save_weights("./model/specificModel/%sModel%d.h5" % (self.name,epoch))


        print('\r auc_val: %s ' %str(round(auc_val, 4)), end=100 * ' ' + '\n')
        print('\r aupr_val: %s ' % str(round(aupr_val, 4)), end=100 * ' ' + '\n')
        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return





t1 = datetime.now().strftime('%Y-%m-%d-%H:%M:%S')


names = ['GM12878', 'HUVEC', 'HeLa-S3', 'IMR90', 'K562', 'NHEK','all','all-NHEK']
name=names[0]
#The data used here is the sequence processed by data_processing.py.
Data_dir='/home/hzy/data/%s/'%name
train=np.load(Data_dir+'%s_train.npz'%name)
#test=np.load(Data_dir+'%s_test.npz'%name)
X_en_tra,X_pr_tra,y_tra=train['X_en_tra'],train['X_pr_tra'],train['y_tra']
#X_en_tes,X_pr_tes,y_tes=test['X_en_tes'],test['X_pr_tes'],test['y_tes']


X_en_tra, X_en_val,X_pr_tra,X_pr_val, y_tra, y_val=train_test_split(
    X_en_tra,X_pr_tra,y_tra,test_size=0.05,stratify=y_tra,random_state=250)




model=None
model=get_model()
model.summary()
print ('Traing %s cell line specific model ...'%name)


back = roc_callback(val_data=[X_en_val, X_pr_val, y_val], name=name)
history=model.fit([X_en_tra, X_pr_tra], y_tra, validation_data=([X_en_val, X_pr_val], y_val), epochs=15, batch_size=64,
                  callbacks=[back])

t2 = datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
print("开始时间:"+t1+"结束时间："+t2)
















