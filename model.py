
from keras.layers import *
from keras.models import *
from keras.optimizers import Adam
from keras.regularizers import l1, l2
import numpy as np
MAX_LEN_en = 3000
MAX_LEN_pr = 2000
NB_WORDS = 4097
EMBEDDING_DIM = 100
embedding_matrix = np.load('embedding_matrix.npy')

from keras import backend as K
from keras.engine.topology import Layer, InputSpec
from keras import initializers

class AttLayer(Layer):
    def __init__(self, attention_dim):
        # self.init = initializers.get('normal')
        self.init = initializers.RandomNormal(seed=10)
        self.supports_masking = True
        self.attention_dim = attention_dim
        super(AttLayer, self).__init__()

    def build(self, input_shape):
        assert len(input_shape) == 3
        self.W = K.variable(self.init((input_shape[-1], self.attention_dim)))
        self.b = K.variable(self.init((self.attention_dim, )))
        self.u = K.variable(self.init((self.attention_dim, 1)))
        self.trainable_weights = [self.W, self.b, self.u]
        super(AttLayer, self).build(input_shape)

    def compute_mask(self, inputs, mask=None):
        return mask

    def call(self, x, mask=None):
        # size of x :[batch_size, sel_len, attention_dim]
        # size of u :[batch_size, attention_dim]
        # uit = tanh(xW+b)
        uit = K.tanh(K.bias_add(K.dot(x, self.W), self.b))
        ait = K.dot(uit, self.u)
        ait = K.squeeze(ait, -1)

        ait = K.exp(ait)

        if mask is not None:
            # Cast the mask to floatX to avoid float64 upcasting in theano
            ait *= K.cast(mask, K.floatx())
        ait /= K.cast(K.sum(ait, axis=1, keepdims=True) + K.epsilon(), K.floatx())
        ait = K.expand_dims(ait)
        weighted_input = x * ait
        output = K.sum(weighted_input, axis=1)

        return output

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[-1])





def get_model():
    enhancers=Input(shape=(MAX_LEN_en,))
    promoters=Input(shape=(MAX_LEN_pr,))

    emb_en=Embedding(NB_WORDS,EMBEDDING_DIM,weights=[embedding_matrix],trainable=True)(enhancers)
    emb_pr=Embedding(NB_WORDS,EMBEDDING_DIM,weights=[embedding_matrix],trainable=True)(promoters)

    enhancer_conv_layer = Conv1D(filters = 64,kernel_size = 40,padding = "valid",activation='relu')(emb_en)
    enhancer_max_pool_layer = MaxPooling1D(pool_size = 20, strides = 20)(enhancer_conv_layer)

    promoter_conv_layer = Conv1D(filters = 64,kernel_size = 40,padding = "valid",activation='relu')(emb_pr)
    promoter_max_pool_layer = MaxPooling1D(pool_size = 20, strides = 20)(promoter_conv_layer)

    merge_layer=Concatenate(axis=1)([enhancer_max_pool_layer, promoter_max_pool_layer])
    bn=BatchNormalization()(merge_layer)
    dt=Dropout(0.5)(bn)

    l_gru = Bidirectional(GRU(50, return_sequences=True))(dt)
    l_att = AttLayer(50)(l_gru)

    preds = Dense(1, activation='sigmoid')(l_att)


    model=Model([enhancers,promoters],preds)
    model.compile(loss='binary_crossentropy',optimizer='adam')
    return model