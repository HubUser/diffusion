from multiprocessing import Pool
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
    
# void get_diffusion(float energy, float ratio, float kp, float* coefA, float* coefB);

lib = ctypes.cdll.LoadLibrary('../diff/libdiffusion.so')
get_diffusion = lib.get_diffusion
get_diffusion.restype = None
get_diffusion.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, 
        ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
        ndpointer(ctypes.c_float, flags="C_CONTIGUOUS")]
coefA = np.zeros((4, 4), dtype=np.float32)
coefB = np.zeros((4, 4), dtype=np.float32)
en = [ 100.0, 1020.0, 10.0, 1.0, 0.1 ];
s = [ 5.6, 4.5, 5.0, 7.0, 3.5, 4.7, 4.68 ];

# [kp][lat]
coefB[0, 0] = 0.002280709100887
coefB[0, 1] = 0.000007740005458
coefB[0, 2] = 0.000000043666208
coefB[0, 3] = 0.000000011799739

coefB[1, 0] = -0.000006133966053
coefB[1, 1] = 0.000004421397080
coefB[1, 2] = 0.000000357794391
coefB[1, 3] = -0.000000001681957

coefB[2, 0] = 0.000003074870165
coefB[2, 1] = 0.000000691596767
coefB[2, 2] = -0.000000052153258
coefB[2, 3] = 0.000000000679970

coefB[3, 0] = -0.000000065216319
coefB[3, 1] = -0.000000019691972
coefB[3, 2] = 0.000000001242614
coefB[3, 3] = -0.000000000019907

coefA[0, 0] = -2.962275028228760
coefA[0, 1] = 0.005432456266135
coefA[0, 2] = 0.000125454942463
coefA[0, 3] = -0.000001982993581

coefA[1, 0] =  -0.022489553317428
coefA[1, 1] = 0.025674106553197
coefA[1, 2] = -0.000471087114420
coefA[1, 3] = 0.000001404062914

coefA[2, 0] = 0.004958000034094
coefA[2, 1] = -0.001571639673784
coefA[2, 2] = 0.000026058032745
coefA[2, 3] = -0.000000054605785

coefA[3, 0] = -0.000074383700849
coefA[3, 1] = 0.000023944545319
coefA[3, 2] = -0.000000390142418
coefA[3, 3] = 0.000000000685272
get_diffusion(1., 4.5, 5., coefA, coefB);
