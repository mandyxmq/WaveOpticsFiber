//
//  Solver.cpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//
#include "Solver.hpp"
#include <math.h>

#define realnum 29971
extern double hankelre[TABLE_SIZE];
extern double hankelim[TABLE_SIZE];
extern double dhankelre[TABLE_SIZE];
extern double dhankelim[TABLE_SIZE];

const std::complex<double> cunit(0,1);
namespace sp = sp_bessel;

// ellipse, radius is longer axis, radius2 is shorter axis
Solver::Solver(int numel, double radius, double radius2, int quadrature, double phi_i, char mode,
               double freq, double mur, std::complex<double> epsr, double angle){
    this->type = 'e';
    this->Nscatter = 1;
    this->numel = numel;
    this->radius = radius;
    this->radius2 = radius2;
    this->centers.push_back(Eigen::Vector3d(0,0,0));
    CreateNodes();

    this->angle = angle;
    this->quadrature = quadrature;
    NumericPar();
    CreateBasis();
    CommonParameter(mode, mur, epsr, phi_i, freq, numel, angle);

    vvec = Eigen::VectorXcd::Zero(4*numel);
    this->Z = Eigen::MatrixXcd::Zero(4*numel, 4*numel);
}

// Arbitrary geometry
Solver::Solver(std::vector<Eigen::Vector3d> nodes, int quadrature, double phi_i, char mode,
               double freq, double mur, std::complex<double> epsr, double angle){
    this->type = 'e'; // set type to ellipse for now
    this->Nscatter = 1;
    this->numel = nodes.size();
    assert(numel>0  && "Invalid nodes input!");
    this->centers.push_back(Eigen::Vector3d(0,0,0));
    this->nodes = nodes;
    // assume center is origin find the longest distance and set it as radius.
    this->radius = CalculateR();
    this->radius2 = radius;

    this->angle = angle;
    this->quadrature = quadrature;
    NumericPar();
    CreateBasis();
    CommonParameter(mode, mur, epsr, phi_i, freq, numel, angle);

    vvec = Eigen::VectorXcd::Zero(4*numel);
    this->Z = Eigen::MatrixXcd::Zero(4*numel, 4*numel);
}

double Solver::CalculateR(){
    double r = nodes[0].norm();
    for (int i=1; i<numel; ++i){
        if (nodes[i].norm()>r){
            r = nodes[i].norm();
        }
    }
    return r;
}

void Solver::CreateNodes(){
    int numel1 = this->numel;
    int n = this->Nscatter;
    // ellipse
    // circle cross section
    // Create Nodes
    double unit = 2*M_PI/numel1;
    this->nodes = std::vector<Eigen::Vector3d>(numel1*n);
    for (int i=0; i<n; ++i){
      Eigen::Vector3d c = centers[i];
      for (int j =0; j < numel1; ++j){
        this->nodes[i*numel1+j] = Eigen::Vector3d(c(0)+radius*cos(j*unit),c(1)+radius2*sin(j*unit),0);
      }
    }
}

void Solver::NumericPar(){
    if (quadrature == 2){
        xvec = {-sqrt(1./3),sqrt(1./3)};
        weights = {1,1};
        CreateEl();
    }
    else if(quadrature == 3){
        xvec = {0, -sqrt(3./5), sqrt(3./5)};
        weights = {8./9, 5./9, 5./9};
        CreateEl();
    }
    else if (quadrature ==4){
        xvec = {sqrt(3./7-2./7*sqrt(6/5)), -sqrt(3./7-2./7*sqrt(6/5)),
            sqrt(3./7+2./7*sqrt(6/5)), -sqrt(3./7+2./7*sqrt(6/5))};
        weights = {(18.+sqrt(30))/36, (18.+sqrt(30))/36, (18.-sqrt(30))/36, (18.-sqrt(30))/36};
        CreateEl();
    }
}

void Solver::CreateEl(){
    int numel1 = this->numel;
    for (int j=0; j<this->Nscatter; ++j){
        int pre = j*numel1;
        el.push_back(Element(pre+0, pre+1, nodes[pre+0], nodes[pre+1], {pre+numel1-1, pre+0},xvec));
        for (int i = 1; i < numel1-1; ++i){
            el.push_back(Element(pre+i, pre+i+1, nodes[pre+i],nodes[pre+i+1],{pre+i-1,pre+i},xvec));
        }
        el.push_back(Element(pre+numel1-1,pre+0, nodes[pre+numel-1], nodes[pre+0], {pre+numel1-2,pre+numel1-1},xvec));
    }
}

void Solver::CreateBasis(){
    int n = this->Nscatter;
    // Create Basis
    ba = std::vector<Basis>(numel*n);
    for (int j=0; j<n; ++j){
        for (int i=0; i<numel; ++i){
            ba[j*numel+i] = Basis(j*numel+i,j*numel+((i+1)%numel));
        }
    }
}


void Solver::CommonParameter(char mode, double mur, std::complex<double> epsr, double phi_i,
                             double freq, int numel, double angle){
    this->mode = mode;
    b = Background(mur, epsr);
    w = Wave(b,phi_i,mode,freq, angle);
    // initalize basis matrix size
    f = Eigen::MatrixXd::Zero(2, quadrature);
    for (int i = 0; i<quadrature; ++i){
        f(0,i) = (1-xvec[i])/2;
        f(1,i) = (1+xvec[i])/2;
    }
    g = f;
}

Solver::~Solver(){}

void Solver::print_elqp() const{
    for (int i=0; i<el.size(); ++i){
        el[i].print_qp();
    }
}

void Solver::print_elnormal() const{
    for (int i=0; i<el.size(); ++i){
        el[i].print_normal();
    }
}

void Solver::print_eltan() const{
    for (int i=0; i<el.size(); ++i){
        std::cout<<el[i].gettan()<<std::endl;
    }
}

Eigen::VectorXcd Solver::get_sol() const{
    return sol;
}


Eigen::Vector3cd crossproduct(Eigen::Vector3cd a, Eigen::Vector3cd b){
  Eigen::Vector3cd result;
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = - (a[0] * b[2] - a[2] * b[0]);
  result[2] = a[0] * b[1] - a[1] * b[0];

  return result; 
}

void Solver::updatewave(double phi_i1, char mode1, double freq1, double angle1){
  mode = mode1;
  w = Wave(b, phi_i1, mode1, freq1, angle1);
}


void Solver::incident_helper(int b1, const Eigen::Vector3d& d2, const Eigen::Vector3d& d3,
                             const Eigen::Vector3d& tan, std::complex<double> constant, double eb, double hb, Eigen::Vector3d nor){
  int n = this->Nscatter*this->numel;
  double tmp = 1./(b.Z0);
  if (mode =='M'){
    //Evvec_z
    vvec(n+b1) += eb*constant*d2.dot(Eigen::Vector3d(0,0,1));
    //Hvvec_t
    vvec(2*n+b1) += tmp * hb*constant*d3.dot(tan);
    //Evvec_t
    vvec(b1) += eb*constant*d2.dot(tan);
    //Hvvec_z
    vvec(3*n+b1) += tmp * hb*constant*d3.dot(Eigen::Vector3d(0,0,1));
  }else{
    //Evvec_z
    vvec(n+b1) += eb*constant*-d3.dot(Eigen::Vector3d(0,0,1));
    //Hvvec_t
    vvec(2*n+b1) += tmp * hb*constant*d2.dot(tan);
    //Evvec_t
    vvec(b1) += eb*constant*-d3.dot(tan);
    //Hvvec_z
    vvec(3*n+b1) += tmp * hb*constant*d2.dot(Eigen::Vector3d(0,0,1));
  }
}

void Solver::incident(){
  // Evvec is Evvec_z; Hvvec is Hvvec_t; Evvec_t, Hvvec_z
  // d1 is wave vector direction; d2 is E direction; d3 is H direction (TM)
  int n = this->Nscatter*this->numel;
  vvec.setZero();

  Eigen::Vector3d d1(sin(angle)*cos(w.phi_i),sin(angle)*sin(w.phi_i),cos(angle));
  Eigen::Vector3d d2(-cos(angle)*cos(w.phi_i),-cos(angle)*sin(w.phi_i),sin(angle));
  Eigen::Vector3d d3 = d2.cross(d1);
  std::complex<double> constant;
  Eigen::Vector3d tan,nor;
  double len,eb,hb;
  int b1;
  for (int i = 0; i<Nscatter*numel; ++i){
    tan = el[i].gettan();
    len = el[i].getlen();
    for (int j=0; j<2; ++j){
      b1 = el[i].getbasis(j);
      for (int k=0; k<quadrature; ++k){
        double tmp = len/2.*weights[k];
        constant = tmp * exp(cunit*(w.k0)*((el[i].getqp(k)).dot(d1)));
        nor = el[i].getnormal(k);
        eb = f(j,k);
        hb = g(j,k);
        incident_helper(b1, d2, d3, tan, constant, eb, hb, nor);
      }
    }
  }
}

void Solver::computehankel(std::complex<double> arg1, std::complex<double> arg2, std::complex<double> &temp1, std::complex<double> &temp2, std::complex<double> &temp3, std::complex<double> &temp4){
  double hankelunit = 0.1; // imageunit = 0.1 too
  int floor, ceil, imagfloor, imagceil;
  int num = 29971;
  double weight, imagweight, A, B, hre, him, dhre, dhim;

  // for small arguement, call library directory
  if (std::real(arg1)<=3){
    temp1 = sp::hankelH2(0,arg1);
    temp2 = sp::hankelH2(1,arg1);
  }else{
    double re = std::real(arg1);
    floor = (int)std::floor((re - 3)/hankelunit);
    ceil = (int)std::ceil((re - 3)/hankelunit);
    weight = (re - 3)/hankelunit-floor;

    hre = (weight * hankelre[ceil] + (1-weight) * hankelre[floor]);
    him = (weight * hankelim[ceil] + (1-weight) * hankelim[floor]);
    temp1.real(hre);
    temp1.imag(him);
    temp1 *= std::exp(-cunit*arg1);

    dhre = (weight * dhankelre[ceil] + (1-weight) * dhankelre[floor]);
    dhim = (weight * dhankelim[ceil] + (1-weight) * dhankelim[floor]);
    temp2.real(dhre);
    temp2.imag(dhim);
    temp2 *= std::exp(-cunit*arg1);
  }

  if (std::real(arg2)<=3){
    temp3 = sp::hankelH2(0,arg2);
    temp4 = sp::hankelH2(1,arg2);
  }else{
    double im = std::imag(arg2);
    imagfloor = (int)std::floor(-im/hankelunit);
    imagceil = (int)std::ceil(-im/hankelunit);
    imagweight = -im / hankelunit - imagfloor;

    double re = std::real(arg2);
    floor = (int)std::floor((re - 3)/hankelunit);
    ceil = (int)std::ceil((re - 3)/hankelunit);
    weight = (re - 3)/hankelunit-floor;

    A = (1-weight) * hankelre[imagfloor*realnum + floor] + weight * hankelre[imagfloor*realnum + ceil];
    B = (1-weight) * hankelre[imagceil*realnum + floor] + weight * hankelre[imagceil*realnum + ceil];
    hre = (1-imagweight) * A + imagweight * B;

    A = (1-weight) * hankelim[imagfloor*realnum + floor] + weight * hankelim[imagfloor*realnum + ceil];
    B = (1-weight) * hankelim[imagceil*realnum + floor] + weight * hankelim[imagceil*realnum + ceil];
    him = (1-imagweight) * A + imagweight * B;

    temp3.real(hre);
    temp3.imag(him);
    temp3 *= std::exp(-cunit*arg2);

    A = (1-weight) * dhankelre[imagfloor*realnum + floor] + weight * dhankelre[imagfloor*realnum + ceil];
    B = (1-weight) * dhankelre[imagceil*realnum + floor] + weight * dhankelre[imagceil*realnum + ceil];
    dhre = (1-imagweight) * A + imagweight * B;

    A = (1-weight) * dhankelim[imagfloor*realnum + floor] + weight * dhankelim[imagfloor*realnum + ceil];
    B = (1-weight) * dhankelim[imagceil*realnum + floor] + weight * dhankelim[imagceil*realnum + ceil];
    dhim = (1-imagweight) * A + imagweight * B;

    temp4.real(dhre);
    temp4.imag(dhim);
    temp4 *= std::exp(-cunit*arg2);

  }
}



void Solver::assembly_same_Di(int b1, int b2, double f1, double f2, double g1, double g2,
                              Eigen::Vector3d tan1, Eigen::Vector3d tan2, Eigen::Vector3d nor, double pow1,
                              double pow2, double constant,std::complex<double> singularity1, std::complex<double> singularity2,
                              std::complex<double> singularity3, std::complex<double> singularity4){
  int n = this->Nscatter*this->numel;
  // ZZZ
  // outside
  std::complex<double> common1 = constant*f1;
  std::complex<double> outside = w.omega*b.mu0/4*common1*singularity1 - w.kz*w.kz/(4.0*w.omega*b.eps0) * common1*singularity1;
  std::complex<double> inside = w.omega*b.mu1/4*common1*singularity3 - w.kz*w.kz/(4.0*w.omega*b.eps1) * common1*singularity3;
  Z(n+b1, n+b2) += outside;
  Z(n+b1, n+b2) += inside;
  Z(3*n+b1, 3*n + b2) += b.eps0/b.mu0*outside;
  Z(3*n+b1, 3*n + b2) += b.eps1/b.mu1*inside;

  // // ZZT use KZT sum to zero
  // // ZTZ = - KTZ (so take the negative of KTZ), sum to zeros

  // tan dot normal component
  // outside
  std::complex<double> temp = f1*tan1.dot(Eigen::Vector3d(0,0,1).cross(g2*tan2));
  outside = w.kz/4.0 * constant*temp*singularity2;
  // inside
  inside = w.kz/4.0 * constant*temp*singularity4;
  Z(b1,2*n+b2) += outside;
  Z(b1,2*n+b2) += inside;
  Z(2*n+b1,b2) += -outside;
  Z(2*n+b1,b2) += -inside;

  // ZTT
  // outside
  // first part of the formula
  common1 = w.omega*constant*g1;
  std::complex<double> common2 = - 1./(4*w.omega)*constant*pow1*pow2;
  outside = b.mu0/4.0 * common1*singularity1*tan1.dot(tan2) + 1.0/b.eps0*common2*singularity2;
  Z(b1,b2) += outside;
  Z(2*n+b1, 2*n+b2) += b.eps0/b.mu0 * outside;
  // inside
  //first part of the formula
  inside = b.mu1/4.0 * common1*singularity3*tan1.dot(tan2) + 1.0/b.eps1*common2*singularity4;
  Z(b1,b2) += inside;
  Z(2*n+b1, 2*n+b2) += b.eps1/b.mu1 * inside;

  //oblqiue specific
  std::complex<double> common3 = -cunit * w.kz / (4.0 * w.omega) * constant * pow1;
  outside = 1.0/b.eps0 * common3 * singularity1;
  inside = 1.0/b.eps1 * common3 * singularity3;
  Z(b1, n+b2) += outside;
  Z(b1, n+b2) += inside;
  Z(2*n+b1, 3*n+b2) += b.eps0/b.mu0*outside;
  Z(2*n+b1, 3*n+b2) += b.eps1/b.mu1*inside;

  common3 = cunit * w.kz / (4.0 * w.omega) * constant * pow2;
  outside = 1.0/b.eps0 * common3 * singularity2 * f1;
  inside = 1.0/b.eps1 * common3 * singularity4 * f1;
  Z(n+b1, b2) += outside;
  Z(n+b1, b2) += inside;
  Z(3*n+b1, 2*n+b2) += b.eps0/b.mu0*outside;
  Z(3*n+b1, 2*n+b2) += b.eps1/b.mu1*inside;
}

void Solver::assembly_diff_Di(int b1, int b2, double f1, double f2, double g1, double g2, double constant1, double constant2, double pow1, double pow2, Eigen::Vector3d tan1, Eigen::Vector3d tan2, Eigen::Vector3d R,
                              std::complex<double> temp1, std::complex<double> temp2, std::complex<double> temp3, std::complex<double> temp4){
  int n = this->Nscatter*this->numel;
  // ZZZ
  // outside
  std::complex<double> common1 = constant1*f1*constant2*f2;
  std::complex<double> outside = (w.omega)*(b.mu0)/4*common1*temp1 - w.kz*w.kz/(4.0*w.omega*b.eps0) * common1*temp1;
  std::complex<double> inside = (w.omega)*(b.mu1)/4*common1*temp3 - w.kz*w.kz/(4.0*w.omega*b.eps1) * common1*temp3;
  Z(n+b1, n+b2) += outside;
  Z(n+b1, n+b2) += inside;
  Z(3*n+b1, 3*n + b2) += b.eps0/b.mu0*outside;
  Z(3*n+b1, 3*n + b2) += b.eps1/b.mu1*inside;

  // ZZT use KZT
  // gradient of H0 gives k*sin(theta) in front
  // outside
  Eigen::Vector3d c = R.cross(g2*tan2);
  common1 = cunit/4.0*constant1*f1*constant2;
  Z(n+b1, 2*n+b2) += common1 * w.newk * temp2 * c(2);
  Z(n+b1, 2*n+b2) += common1 * w.newk1 * temp4 * c(2);
  Z(3*n+b1, b2) += -common1 * w.newk * temp2 * c(2);
  Z(3*n+b1, b2) += -common1 * w.newk1 * temp4 * c(2);

  // ZTZ = - KTZ (so take the negative of KTZ)
  // gradient of H0 gives k*sin(theta) in front
  // outside
  double cdot = tan1.dot(R.cross(Eigen::Vector3d(0,0,f2)));
  Z(b1, 3*n+b2) += common1 * w.newk * cdot * temp2;
  Z(b1, 3*n+b2) += common1 * w.newk1 * cdot * temp4;
  Z(2*n+b1, n+b2) += -common1 * w.newk * cdot * temp2;
  Z(2*n+b1, n+b2) += -common1 * w.newk1 * cdot * temp4;

  // tan dot normal component
  // outside
  cdot = tan1.dot(Eigen::Vector3d(0,0,1).cross(g2*tan2));
  outside = w.kz/4.0*constant1*f1*constant2*cdot*temp1;
  inside = w.kz/4.0*constant1*f1*constant2*cdot*temp3;
  Z(b1,2*n+b2) += outside;
  Z(b1,2*n+b2) += inside;
  Z(2*n+b1,b2) += -outside;
  Z(2*n+b1,b2) += -inside;

  // ZTT
  // outside
  common1 = constant1*g1*constant2*g2;
  std::complex<double> common2 = - 1.0/(4.0*w.omega)*constant1*constant2*pow1*pow2;
  outside = w.omega*b.mu0/4.0*common1*tan1.dot(tan2)*temp1 + 1.0/b.eps0*common2*temp1;
  Z(b1,b2) += outside;
  Z(2*n+b1, 2*n+b2) += b.eps0/b.mu0 * outside;

  inside = w.omega*b.mu1/4.0*common1*tan1.dot(tan2)*temp3 + 1.0/b.eps1*common2*temp3;
  Z(b1,b2) += inside;
  Z(2*n+b1, 2*n+b2) += b.eps1/b.mu1 * inside;

  //oblqiue specific
  // gradient
  common1 = -cunit * w.kz / (4.0 * w.omega) * constant1 * constant2 * pow1 * f2;
  outside = 1.0/b.eps0 * common1 * temp1 ;
  inside = 1.0/b.eps1 * common1 * temp3 ;
  Z(b1, n+b2) += outside;
  Z(b1, n+b2) += inside;
  Z(2*n+b1, 3*n+b2) += b.eps0/b.mu0*outside;
  Z(2*n+b1, 3*n+b2) += b.eps1/b.mu1*inside;

  common1 = cunit * w.kz / (4.0 * w.omega) * constant1 * constant2 * f1 * pow2;
  outside = 1.0/b.eps0 * common1 * temp1;
  inside = 1.0/b.eps1 * common1 * temp3;
  Z(n+b1, b2) += outside;
  Z(n+b1, b2) += inside;
  Z(3*n+b1, 2*n+b2) += b.eps0/b.mu0*outside;
  Z(3*n+b1, 2*n+b2) += b.eps1/b.mu1*inside;
}


// calculate energy balance
double Solver::energy(){
  double energy = 0;
  int n = this->Nscatter*this->numel;
  Eigen::Vector3d tanvec, nor;
  std::complex<double> J_coef_z, J_coef_t, M_coef_t, M_coef_z;
  Eigen::Vector3cd Jz, Mz, Jt, Mt, Et, Ht;
  std::vector<Eigen::Vector3cd>J(2), M(2);
  J[0] = Eigen::Vector3cd::Zero(3);
  J[1] = Eigen::Vector3cd::Zero(3);
  M[0] = Eigen::Vector3cd::Zero(3);
  M[1] = Eigen::Vector3cd::Zero(3);
  for(int j = 0; j < numel; ++j){
    tanvec = el[j].gettan();
    double len = el[j].getlen();
    J[0].setZero();
    J[1].setZero();
    M[0].setZero();
    M[1].setZero();
    for (int k = 0; k<2; ++k){
      int b1 = el[j].getbasis(k);
      J_coef_z = sol(n+b1);
      J_coef_t = sol(b1);
      M_coef_t = sol(2*n+b1);
      M_coef_z = sol(3*n+b1);

      for(int p = 0; p<quadrature; ++p){
        Jz = J_coef_z * f(k,p) * Eigen::Vector3d(0,0,1);
        Mz = M_coef_z * g(k,p) * Eigen::Vector3d(0,0,1);
        Jt = J_coef_t * f(k,p) * tanvec;
        Mt = M_coef_t * g(k,p) * tanvec;

        J[k] += Jz + Jt;
        M[k] += Mz + Mt;
      }
    }
    nor = Eigen::Vector3d{tanvec[1], -tanvec[0], 0};

    Eigen::Vector3cd poynting = crossproduct(J[0].conjugate(), M[0]).real() + crossproduct(J[1].conjugate(), M[1]).real();
    energy += std::real(nor.dot(poynting)) * len/2.0 * weights[0];
  }
  energy *= b.Z0 / (2*radius*sin(angle));;
  return energy;
};

// calculate scattered field from current;
void Solver::BSDF(int nTheta, double Dis, Eigen::VectorXf& sigma){
  // TM
  Eigen::Vector3cd d1(sin(angle)*cos(w.phi_i),sin(angle)*sin(w.phi_i),cos(angle));
  Eigen::Vector3cd d2(-cos(angle)*cos(w.phi_i),-cos(angle)*sin(w.phi_i),sin(angle));
  Eigen::Vector3cd d3 = d2.cross(d1);

  std::complex<double> hankel, J_coef_z, J_coef_t, M_coef_t, M_coef_z;
  Eigen::Vector3cd Jt, Mt, dhankel, tmp, E_s, H_s;
  Eigen::Vector3d p1, diffvec, R;
  double len1, pow1, dis, common;

  int n = this->Nscatter*this->numel;

  double the, x, y;
  Eigen::Vector3d tanvec, nor;
  for (int i = 0; i < nTheta; ++i){
    the = i * 2 * M_PI / nTheta;
    x = Dis*cos(the);
    y = Dis*sin(the);
    Eigen::Vector3d dir(cos(the), sin(the),0);

    Eigen::Vector3cd E(0, 0, 0);
    Eigen::Vector3cd H(0, 0, 0);

    Eigen::Vector3cd propdir(sin(angle)*dir(0), sin(angle)*dir(1), -cos(angle));

    for(int j = 0; j < numel; ++j){
      tanvec = el[j].gettan();
      len1 = el[j].getlen();
      for (int k = 0; k<2; ++k){
        //for (int k = 0; k<1; ++k){
        int b1 = el[j].getbasis(k);
        pow1 = pow(-1,k+1)/len1;

        J_coef_z = sol(n+b1);
        J_coef_t = sol(b1);
        M_coef_t = sol(2*n+b1);
        M_coef_z = sol(3*n+b1);

        for(int p = 0; p<quadrature; ++p){
          p1 = el[j].getqp(p);
          Jt = J_coef_t * f(k,p) * tanvec;
          Mt = M_coef_t * g(k,p) * tanvec;

          nor = el[j].getnormal(p);
          diffvec = Eigen::Vector3d(x,y,0) - p1;
          dis = diffvec.norm();
          R = diffvec / dis;

          hankel = std::sqrt(2.0*cunit/(M_PI*w.newk*dis)) * exp(-cunit*w.newk*dis);

          common = el[j].getlen()/2.0 * weights[p];
          tmp = -w.omega*b.mu0/4.0 *  common * Eigen::Vector3d(0,0,1) * hankel;
          E += J_coef_z *f(k,p) * tmp;
          H += b.eps0/b.mu0 * M_coef_z *g(k,p) * tmp;

          tmp = -w.omega*b.mu0/4.0 * common * tanvec * hankel;
          E += J_coef_t * f(k,p) * tmp;
          H += b.eps0/b.mu0 * M_coef_t * g(k,p) * tmp;
        }
      }
    }
    // substract the component that is projected on to the propagation direction
    E_s = E - propdir.dot(E) * propdir;
    H_s = H - propdir.dot(H) * propdir;
    E_s += crossproduct(-propdir, H_s) * b.Z0;
    H_s = crossproduct(propdir, E_s) / b.Z0;

    // scattered
    //sigma(i) = std::real(Eigen::Vector3d(dir(0), dir(1), 0).dot(wavevec)) / sin(angle);
    sigma(i) = (float) E_s.squaredNorm();
  }
  sigma *=  (float) (Dis / (2*radius));
};
