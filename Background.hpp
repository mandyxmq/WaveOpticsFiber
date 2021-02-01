//
//  Background.hpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//
#ifndef Background_hpp
#define Background_hpp

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <complex>

class Background{
public:
  double c0 = 299792458;
  double mu0 = 4e-7 * M_PI;
  double eps0 = 8.8541878e-12;
  double Z0 = sqrt(mu0/eps0);

  // parameters for scatterer
  std::complex<double> c1;
  double mu1;
  std::complex<double> eps1;
  std::complex<double> Z1;
    
  Background():c1(0), mu1(0), eps1(0), Z1(0){};
  Background(double mur, std::complex<double> epsr){
    this->mu1 = mu0 * mur;
    this->eps1 = eps0 * epsr;
    this->c1 = (double)1./sqrt(mu1*eps1);
    this->Z1 = sqrt(mu1/eps1);
  }
    
};
#endif /* Background_hpp */
