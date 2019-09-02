# EPIVAN(Promoter-Enhancer Interaction predictor with pre-trained Vector and Attention based neural Networks)
EPIVAN is a new deep learning method that enables predicting long-range EPIs using only genomic sequences.

The three contributions of this work are as follows: (1) We use DNA vectors pre-trained with human whole genome sequences to encode the enhancers and promoters. (2) We use the attention mechanism to boost the contribution of key features, thus improv-ing the performance of the model. (3) We build a general model, which has transfer ability and can be used to predict EPIs in various cell lines. 

# File Description 
- Data_Augmentation.R

  A tool of data augmentation provided by Mao et al. (2017). The details of the tool can be seen in https://github.com/wgmao/EPIANN.

  We used this tool to amplify the positive samples in the training set to 20 times to achieve class balance.

- sequence_processing.py

  Perform pre-processing of DNA sequences:

  1.	Convert the enhancer and promoter gene sequences into sequences consisting of words (6-mers), and if a word contains a ‘N’, the word is marked as ‘NULL’.
  2.	Construct a dictionary containing 4^6+1 words.
  3.	Convert each gene sequence into a sequence of word indexes according to the dictionary (each word has its own unique in-dex).
 
- embedding_matrix.npy

  The weight of the embedding layer converted from the pre-trained DNA vector provided by Ng (2017).

- train.py

  Perform model training.

  You can find the weight of the model mentioned in our paper under the directory model/.
  Directory| Content
  --- | ---
  model/specificModel/| the weight of EPIVAN-specific on each cell line.
  model/generalModel/| the weight of EPIVAN-general.
  model/retrainModel/| the weight of EPIVAN-best on each cell line.
  model/transferModel/| the weight of EPIVAN-general transferred to the new cell line.


- test.py

  Evaluate the performance of model.



References:

  Mao, W. et al. (2017) Modeling Enhancer-Promoter Interactions with Attention-Based Neural Networks. bioRxiv, 219667.

  Ng, P. (2017) dna2vec: Consistent vector representations of variable-length k-mers. arXiv:1701.06279.
