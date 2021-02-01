//
//  Cylinder.hpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//

#ifndef MoM_ob_hpp
#define MoM_ob_hpp

#include <stdio.h>
#include <vector>
#include <complex>
#include "Background.hpp"
#include "Wave.hpp"
#include "Element.hpp"
#include "Solver.hpp"

class MoM_ob: public Solver{
private:
 
public:
  int nTheta;

  MoM_ob():Solver(){};

  MoM_ob(int numel, double radius, double radius2, int quadrature, double phi_i, char mode,
         double freq, double mur, std::complex<double> epsr, double angle):
    Solver(numel, radius, radius2, quadrature, phi_i, mode, freq, mur, epsr, angle){};

  MoM_ob(std::vector<Eigen::Vector3d> nodes, int quadrature, double phi_i, char mode,
         double freq, double mur, std::complex<double> epsr, double angle):
    Solver(nodes, quadrature, phi_i, mode, freq, mur, epsr, angle){}

  ~MoM_ob();
  void assembly();
};

#endif /* MoM_ob_hpp */

