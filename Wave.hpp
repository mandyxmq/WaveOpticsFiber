//
//  Wave.hpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//

#ifndef Wave_hpp
#define Wave_hpp

#include <stdio.h>
#include "Background.hpp"
class Wave{
public:
  double phi_i;       // incident angle
  char   mode;       // 'M' for TM mode, 'E' for TE mode
  double freq;
  double wavelength;  // wavelength = c0/freq;
  double k0;          // k0 = 2*pi/wavelength
  double newk;
  double kz;
  double omega;

  // parameters for scatterer
  std::complex<double> wavelength1;
  std::complex<double> k1;
  std::complex<double> newk1;
    
public:
  Wave():phi_i(0), mode('\0'), freq(0), wavelength(0), k0(0), omega(0), wavelength1(0), k1(0){};
  Wave(Background b, double phi_i, char mode,double freq, double angle){
    this->phi_i = phi_i;
    this->mode = mode;
    this->freq = freq;
    omega = 2*M_PI*freq;
    
    this->wavelength = b.c0 / freq;
    this->k0 = (double)2*M_PI / wavelength;
    this->newk = k0 * sin(angle);
    this->kz = k0 * cos(angle);
    
    if (std::real(b.c1)!=0||std::imag(b.c1)!=0){
      this->wavelength1 = b.c1 / freq;
      this->k1 = 2 * M_PI / wavelength1;
      this->newk1 = std::sqrt(k1*k1-k0*k0*cos(angle)*cos(angle));
    }
  };
  ~Wave(){};
};

#endif /* Wave_hpp */
