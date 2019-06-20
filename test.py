

import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
from model import get_model
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score


models=['GM12878', 'HUVEC', 'HeLa-S3', 'IMR90', 'K562', 'NHEK']
m=models[0]
model=None
model=get_model()
model.load_weights("./model/specificModel/%sModel.h5" % m)

names = ['GM12878', 'HUVEC', 'HeLa-S3', 'IMR90', 'K562', 'NHEK']
for name in names:
    Data_dir='/home/hzy/data/%s/'%name
    test=np.load(Data_dir+'%s_test.npz'%name)
    X_en_tes,X_pr_tes,y_tes=test['X_en_tes'],test['X_pr_tes'],test['y_tes']

    print("****************Testing %s cell line specific model on %s cell line****************"%(m,name))
    y_pred = model.predict([X_en_tes,X_pr_tes])
    auc=roc_auc_score(y_tes, y_pred)
    aupr=average_precision_score(y_tes, y_pred)

    print("AUC : ", auc)
    print("AUPR : ", aupr)



