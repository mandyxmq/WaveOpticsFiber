//
//  Solver.hpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//

#ifndef Solver_hpp
#define Solver_hpp

#include <stdio.h>
#include "Background.hpp"
#include "Wave.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include <complex>
#include <assert.h>
#include <iomanip>
#include <complex_bessel.h>

#define TABLE_SIZE 29971000

extern const std::complex<double> cunit;

class Solver{
protected:
  // Scatter parameters
  double radius;      // radius of cylinder or for an arbitrary 2d object: longest distance to the center
  double radius2;     // short axis radius for ellipse
  std::vector<Eigen::Vector3d> centers;   // vector or centers
  std::vector<Eigen::Vector3d> nodes;     // Vector of nodes
  std::vector<Element> el;                // Vector of elements
  std::vector<Basis> ba;                  // Vector of basis

  std::vector<double> xvec;
  std::vector<double> weights;
  char mode;
  char type;
  double angle;

public:
  Wave w;
  int Nscatter;       // number of scatters;
  int numel;          // number of elements/ling segments for circle discretization
  int quadrature;     // number of quadrature points for each line segment
  Eigen::MatrixXd f,g, f_in, g_in;
  Background b;
  Eigen::VectorXcd vvec;
  Eigen::VectorXcd sol;
  Eigen::MatrixXcd Z;

  Solver(){};

  // Constructor for elliptical cross-section
  Solver(int numel, double radius, double radius2, int quadrature,
         double phi_i, char mode, double freq, double mur, std::complex<double> epsr, double angle);

  // Constructor for arbitrary cross-sectional shape
  Solver(std::vector<Eigen::Vector3d> nodes, int quadrature, double phi_i, char mode, double freq, double mur, std::complex<double> epsr, double angle);

  double CalculateR();
  void CreateNodes();
  void NumericPar();
  void CreateEl();
  void CreateBasis();
  void CommonParameter(char mode, double mur, std::complex<double> epsr, double phi_i, double freq, int numel, double angle);

  virtual ~Solver();

  void print_elqp() const;
  void print_elnormal() const;
  void print_eltan() const;

  Eigen::VectorXcd get_sol() const;

  void updatewave(double phi_i1,char mode1, double freq1, double angle1);
  void incident();
  void incident_helper(int b1, const Eigen::Vector3d& c, const Eigen::Vector3d& d2,
                       const Eigen::Vector3d& tan, std::complex<double> constant, double eb, double hb, Eigen::Vector3d nor);

  // MoM, FMM, MLFMA have different implementation for assemby(near assembly)
  virtual void assembly() = 0;
  void computehankel(std::complex<double> arg1, std::complex<double> arg2, std::complex<double> &temp1, std::complex<double> &temp2, std::complex<double> &temp3, std::complex<double> &temp4);
  void assembly_same_Di(int b1, int b2, double f1, double f2, double g1, double g2,
                        Eigen::Vector3d tan1, Eigen::Vector3d tan2, Eigen::Vector3d nor, double pow1, double pow2,
                        double constant, std::complex<double> singularity1, std::complex<double> singularity2,
                        std::complex<double> singularity3, std::complex<double> singularity4);
  void assembly_diff_Di(int b1, int b2, double f1, double f2, double g1, double g2,
                        double constant1, double constant2, double pow1, double pow2,
                        Eigen::Vector3d tan1, Eigen::Vector3d tan2, Eigen::Vector3d R,
                        std::complex<double> temp1, std::complex<double> temp2, std::complex<double> temp3, std::complex<double> temp4);

  void BSDF(int nTheta, double Dis, Eigen::VectorXf& sigma);
  double energy();
};

#endif /* Solver_hpp */
