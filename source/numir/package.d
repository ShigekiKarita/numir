/**
   License: $(LINK2 http://boost.org/LICENSE_1_0.txt, Boost License 1.0).

   Authors: $(LINK2 http://shigekikarita.github.io, karita)
*/

module numir;

/*
  this library is motivated by
  https://github.com/torch/torch7/wiki/Torch-for-Numpy-users
*/

/*

  TODO:
  1. bring random ndarray (normal, uniform, binomial, multinomial) from d-svm
  1. implement Selection, Calculation, Comparison operations
  1. check linalg operations in mir
     https://github.com/ShigekiKarita/numir/commit/08de747b7b51ffd03e1cf0d7e35c83fe287fc20d
  1. add more linalg operations (LU, SVD)
  1. CUDA array from d-nvrtc

*/

public import numir.core;
public import numir.random;
