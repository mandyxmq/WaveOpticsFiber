
//  MoM_ob.cpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//

#include "MoM_ob.hpp"
#include <math.h>
#include <fstream>

const double consb = 1.781/2.0;
namespace sp = sp_bessel;

MoM_ob::~MoM_ob(){}

// dense assembly
void MoM_ob::assembly(){
  double xp,xp1,len1,len2,dis;
  Eigen::Vector3d p1, p2, R;
  double cons = consb * w.newk;
  std::complex<double> cons1 = consb * w.newk1;
  std::complex<double> temp1, temp2, temp3, temp4,
    singularity1, singularity2, singularity3, singularity4;
  Eigen::Vector3d nor, c, tan1, tan2;
  double d, pow1 ,pow2, constant1, constant2;
  int b1, b2;

  for (int i =0; i<Nscatter*numel; ++i){
    len1 = el[i].getlen();
    tan1 = el[i].gettan();
    for (int m = 0; m<2; ++m){
      b1 = el[i].getbasis(m);
      pow1 = pow(-1, m+1)/len1;
      for (int p = 0; p<quadrature; ++p){
        p1 = el[i].getqp(p);
        nor = el[i].getnormal(p);
        constant1 = len1/2*weights[p];
        for (int j = 0; j<Nscatter*numel; ++j){
          len2 = el[j].getlen();
          tan2 = el[j].gettan();
          for (int n = 0; n < 2; ++n){
            b2 = el[j].getbasis(n);
            pow2 = pow(-1,n+1)/len2;
            if (i==j){ // on the same element, need to handle the Green's function singularity
              if (n==0){
                xp= (p1-nodes[el[i].getn2()]).norm();
              }else{
                xp = (p1-nodes[el[i].getn1()]).norm();
              }
              xp1 = (p1-nodes[el[i].getn1()]).norm();
              double factor;
              singularity1 = len2/2-2/M_PI*(xp*xp/(2*len2)*log(xp/(len2-xp))
                                            + len2/2*log(cons*(len2-xp))-xp/2-len2/4)*cunit;
              singularity2= len2-cunit*2.0/M_PI*(len2*log(cons*2.*(len2-xp1)/(2.f*exp(1.f)))
                                               + xp1*log(xp1/(len2-xp1)));
              singularity3 = len2/2-2.0/M_PI*(xp*xp/(2*len2)*log(xp/(len2-xp))
                                             + len2/2*log(cons1*(len2-xp))-xp/2-len2/4)*cunit;
              singularity4 = len2-2.0/M_PI*(len2*log(2.0*cons1*(len2-xp1)/(2.0*exp(1.0)))
                                          + xp1*log(xp1/(len2-xp1)))*cunit;
              assembly_same_Di(b1, b2, f(m,p), f(n,p), g(m,p), g(n,p), tan1, tan2, nor,
                               pow1, pow2, constant1,singularity1, singularity2, singularity3, singularity4);
            }else{ // different elements
              for (int q = 0; q<quadrature; ++q){
                p2 = el[j].getqp(q);
                dis = (p1-p2).norm();
                R = (p1-p2)/dis;
                constant2 = len2/2*weights[q];
                c = tan2.cross(R);
                d = nor.dot(R);

                computehankel(w.newk*dis, w.newk1*dis, temp1, temp2, temp3, temp4);
                assembly_diff_Di(b1, b2, f(m,p), f(n,q), g(m,p), g(n,q), constant1, constant2,
                                 pow1, pow2, tan1, tan2, R, temp1, temp2, temp3, temp4);
              }
            }
          }
        }
      }
    }
  }
}
