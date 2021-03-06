/*
 * Copyright (c) 2012-2013 Haitham Hassanieh, Piotr Indyk, Dina Katabi,
 *   Eric Price, Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H
 
  
  void get_expermient_vs_N_parameters(int N, bool WITH_COMB, double &Bcst_loc, double &Bcst_est,  double &Comb_cst,  int &loc_loops, 
                                      int &est_loops,  int &threshold_loops, int &comb_loops,
                                      double &tolerance_loc, double &tolerance_est);

  void get_expermient_vs_K_parameters(int K, bool WITH_COMB, double &Bcst_loc, double &Bcst_est,  double &Comb_cst,  int &loc_loops,
                                      int &est_loops,  int &threshold_loops, int &comb_loops,
                                      double &tolerance_loc, double &tolerance_est);


#endif
