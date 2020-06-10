
import os, sys
import tensorflow as tf;
import numpy as np

class Repeat_CNN_model:
    def __init__(self, _num_labels, _height, _width, _depth, _gpu_cpu, _num_filter_in_first=64):
       self.prepare_sess = False;
       self.mdtype = 'compat'
       self.mdtype = 'old'
       self.num_filter_in_first = _num_filter_in_first
       self.gpu_cpu = _gpu_cpu

       self.height = _height;
       self.width = _width;
       self.depth = _depth;
       self.num_labels = _num_labels;
       self.X = tf.compat.v1.placeholder("float", [None, _height, _width, _depth], name='X') if self.mdtype=='compat' else tf.placeholder("float", [None, _height, _width, _depth], name='X'); 
       self.Y = tf.compat.v1.placeholder("float", [None, _num_labels], name='Y') if self.mdtype=='compat' else tf.placeholder("float", [None, _num_labels], name='Y'); # _num_labels = 2
       self.learning_rate = tf.compat.v1.placeholder(tf.float32, name='learning_rate') if self.mdtype=='compat' else tf.placeholder(tf.float32, name='learning_rate');
       self.dropout = tf.compat.v1.placeholder(tf.float32, name='dropout') if self.mdtype=='compat' else tf.placeholder(tf.float32, name='dropout');
       
    def weight_variable(self, _shape):
       _initial = tf.random.truncated_normal(_shape, stddev=0.1) if self.mdtype=='compat' else tf.truncated_normal(_shape, stddev=0.1)
       return tf.Variable(_initial)
    def bias_variable(self, _shape):
       _initial = tf.constant(0.1, shape=_shape)
       return tf.Variable(_initial)
    def conv_max_pool(self, _x, _W, _b, _filter):
       c_x = tf.nn.conv2d(_x, _W, strides=[1, 1, 1, 1], padding='SAME')
       c_x = tf.nn.bias_add(c_x, _b);
       c_x = tf.nn.relu(c_x);
       return tf.nn.max_pool(c_x, ksize=[1, 2, _filter, 1], strides=[1, 2, _filter, 1], padding='SAME');
       #return tf.nn.max_pool(c_x, ksize=[1, 2, 1, 1], strides=[1, 2, 1, 1], padding='SAME');
       

    def create_model(self):
       if self.num_filter_in_first==2: num_filter_c =    [  2,  4,  8,  32];
       elif self.num_filter_in_first==4: num_filter_c =  [  4,  8, 16,  64];
       elif self.num_filter_in_first==8: num_filter_c =  [  8, 16, 32, 128];
       elif self.num_filter_in_first==16: num_filter_c = [ 16, 32, 64, 256];
       elif self.num_filter_in_first==32: num_filter_c = [ 32, 64,128, 512];
       elif self.num_filter_in_first==64: num_filter_c = [ 64,128,256, 1024];
       elif self.num_filter_in_first==128: num_filter_c =[128,256,512, 1024]; #[128,256,384, 1024]
       else: num_filter_c = [256,512,768, 1024];
       #else: num_filter_c = [256,512,1024, 1024];
       #num_filter_c = [128,256,512, 2048]
       self.weights = {\
                      'w_c0': self.weight_variable([4, 1, self.depth, num_filter_c[0]]), \
                      'w_c1': self.weight_variable([3, 1, num_filter_c[0], num_filter_c[1]]), \
                      'w_c2': self.weight_variable([2, 1, num_filter_c[1], num_filter_c[2]]) if self.width<7 else self.weight_variable([2, 2, num_filter_c[1], num_filter_c[2]]), \
                      'w_f0': self.weight_variable([7*self.width*num_filter_c[2], num_filter_c[3]]) if self.width<7 else self.weight_variable([7*int(self.width/2+0.5)*num_filter_c[2], num_filter_c[3]]), \
                      'wout': self.weight_variable([num_filter_c[3], self.num_labels])\
                      }
       self.baises = {\
                      'b_c0': self.bias_variable([num_filter_c[0]]), \
                      'b_c1': self.bias_variable([num_filter_c[1]]), \
                      'b_c2': self.bias_variable([num_filter_c[2]]), \
                      'b_f0': self.bias_variable([num_filter_c[3]]), \
                      'bout': self.bias_variable([self.num_labels])\
                     }
       self.cn0 = self.conv_max_pool(self.X, self.weights['w_c0'], self.baises['b_c0'], 1);
       print (self.cn0.get_shape())
       self.cn1 = self.conv_max_pool(self.cn0, self.weights['w_c1'], self.baises['b_c1'], 1);
       print (self.cn1.get_shape())
       self.cn2 = self.conv_max_pool(self.cn1, self.weights['w_c2'], self.baises['b_c2'], 1 if self.width<7 else 2);
       print (self.cn2.get_shape())

       self.cn_flat = tf.reshape(self.cn2, [-1, self.weights['w_f0'].get_shape().as_list()[0]]);
       print (self.cn_flat.get_shape())
       self.fc1 = tf.nn.relu(tf.add(tf.matmul(self.cn_flat, self.weights['w_f0']), self.baises['b_f0']))
       print (self.fc1.get_shape())
       self.fcd = tf.nn.dropout(self.fc1, rate=1-self.dropout) if self.mdtype=='compat' else tf.nn.dropout(self.fc1, 1-self.dropout);

       self.logits = tf.add(tf.matmul(self.fcd, self.weights['wout']),self.baises['bout'], name='logits');
 
    def creat_optimizer(self):
       self.prediction = tf.nn.softmax(self.logits, name='prediction');
       self.pred_labels = tf.argmax(self.prediction,1, name='pred_labels') 

       self.loss_op = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=self.logits, labels=self.Y));
       # for optimizer   
       self.optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=self.learning_rate) if self.mdtype=='compat' else tf.train.AdamOptimizer(learning_rate=self.learning_rate);
       self.train_op = self.optimizer.minimize(self.loss_op);
    
    def create_metric(self):
       self.correct_pred = tf.equal(tf.argmax(self.prediction, 1), tf.argmax(self.Y, 1));
       self.accuracy = tf.reduce_mean(tf.cast(self.correct_pred, tf.float32));

       # auc, precision, and recall.
       #self.auc_op = tf.metrics.auc(self.Y, self.prediction) 
       self.auc_op = tf.compat.v1.metrics.auc(self.Y, self.prediction) if self.mdtype=='compat' else tf.metrics.auc(self.Y, self.prediction)

       #self.precision = tf.metrics.precision(tf.argmax(self.Y, 1), tf.argmax(self.prediction, 1))
       self.precision = tf.compat.v1.metrics.precision(tf.argmax(self.Y, 1), tf.argmax(self.prediction, 1)) if self.mdtype=='compat' else tf.metrics.precision(tf.argmax(self.Y, 1), tf.argmax(self.prediction, 1))

       #self.recall = tf.metrics.recall(tf.argmax(self.Y, 1), tf.argmax(self.prediction, 1))
       self.recall = tf.compat.v1.metrics.recall(tf.argmax(self.Y, 1), tf.argmax(self.prediction, 1)) if self.mdtype=='compat' else tf.metrics.recall(tf.argmax(self.Y, 1), tf.argmax(self.prediction, 1))

    def prepare_run(self):
       self.create_model();
       self.creat_optimizer();
       self.create_metric();

       # initialization of variables
       self.init = tf.compat.v1.global_variables_initializer() if self.mdtype=='compat' else tf.global_variables_initializer();
       self.init_l = tf.compat.v1.local_variables_initializer() if self.mdtype=='compat' else tf.local_variables_initializer();
       self.saver = tf.compat.v1.train.Saver() if self.mdtype=='compat' else tf.train.Saver();
      
       if self.gpu_cpu in ['CPU','cpu']:
          self.config = tf.ConfigProto(device_count={"CPU": 32})
       else: 
          self.config = tf.compat.v1.ConfigProto() if self.mdtype=='compat' else tf.ConfigProto();
          self.config.gpu_options.allow_growth = True
       self.sess = tf.compat.v1.Session(config=self.config) if self.mdtype=='compat' else tf.Session(config=self.config)
       self.prepare_sess = True;

       self.sess.run(self.init);
       self.sess.run(self.init_l);

    def __del__(self):
       if self.prepare_sess: self.sess.close();


if __name__=='__main__':
   c_r3cnn = Repeat_CNN_model(4, 50, 3, 3);

   c_r3cnn.prepare_run();


