//
//  BSDFgen.cpp
//  2Dsolve
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include "MoM_ob.hpp"

#include <string.h>
#include <iomanip>

#include <cstdio>
#include <ctime>
#include <cstdlib>

#include "omp.h"
#include <Eigen/Dense>
#include <sys/time.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <bits/stdc++.h>

double hankelre[TABLE_SIZE];
double hankelim[TABLE_SIZE];
double dhankelre[TABLE_SIZE];
double dhankelim[TABLE_SIZE];

// save scattering disttribution, pdf, cdf as 32bit floats
void postprocess(bool lossless, Eigen::VectorXf& sigma, double energy, int phionum, int vectindex, float unit, Eigen::VectorXf& vect, Eigen::VectorXf& pdf, Eigen::VectorXf& cdf){
  double expected;
  if (lossless)
    expected = 1.f;
  else
    expected = (float) (energy + 1);

  if (expected<=0){
    vect.segment(vectindex, phionum).setZero();
    pdf.segment(vectindex, phionum).array() += (float) (1.0 / (2 * M_PI));
  }else{
    sigma *= expected / (sigma.sum() * unit);
    vect.segment(vectindex, phionum) = sigma;
    float normalization = sigma.sum() * unit;

    if (std::abs(normalization)<1e-6)
      pdf.segment(vectindex, phionum).array() = (float) (1.0 / (2 * M_PI));
    else
      pdf.segment(vectindex, phionum) = sigma / normalization;

  }
  for (int index = 0; index<phionum-1; ++index)
    cdf(vectindex+index) = pdf.segment(vectindex, index+1).sum() * unit;
  cdf(vectindex+phionum-1) = 1;
}

std::vector<Eigen::Vector3d> readgeometry(std::string xfile, std::string yfile){
  std::string line;
  std::ifstream myfile (xfile);
  std::vector<double> xvec;
  std::vector<double> yvec;
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
        {
          xvec.push_back(std::stod(line));
        }
      myfile.close();
    }

  std::ifstream myfile2 (yfile);
  if (myfile2.is_open())
    {
      while ( getline (myfile2,line) )
        {
          yvec.push_back(std::stod(line));
        }
      myfile2.close();
    }
  else std::cout << "Unable to open file";

  std::vector<Eigen::Vector3d> nodes(xvec.size());
  for (int i = 0; i < xvec.size(); ++i){
    nodes[i] = Eigen::Vector3d{xvec[i], yvec[i], 0};
  }

  return nodes;
}

void readiorimag(std::vector<double> &kval, int num, std::string filename){
  std::cout<<"iorfile "<<filename<<std::endl;
  std::ifstream myfile(filename, std::ios::in|std::ios::binary);
  myfile.read((char*)&kval[0], num * sizeof(double));
}

int main(int argc, const char * argv[]) {
  std::string example = argv[1];
  std::string output = example+"/";
  std::cout<<"output directory "<<output<<std::endl;

  std::stringstream ellipsestr(argv[2]);
  bool ellipse;
  ellipsestr >> ellipse;

  double radius1, radius2, etare, etaim;
  int numel, phiinum, phionum, thetanum, lambdanum;
  std::vector<double> nval, kval;
  std::string filename, filename2, filename3;
  std::vector<Eigen::Vector3d> nodes;
  bool lossless = true;
  if (ellipse){

    std::cout<<"------Simulating an elliptical cross-section------"<<std::endl;

    std::string radius1str = argv[3]; // semi-major radius, in micron
    std::string::size_type sz;     // alias of size_t
    radius1 = std::stof (radius1str,&sz);

    std::string radius2str = argv[4]; // semi-minor radius, in micron
    radius2 = std::stof (radius2str,&sz);
    std::cout<<"radius1 "<<radius1<<"um radius2 "<<radius2<<"um"<<std::endl;
    radius1 *= 1e-6;
    radius2 *= 1e-6;

    std::stringstream elstr(argv[5]); // number of elements for the boundary of the cross-section
    elstr >> numel;
    std::cout<<"number of elements "<<numel<<std::endl;

    std::stringstream phionumstr(argv[6]); // number of outgoing phi directions
    phionumstr >> phionum;
    if (radius1==radius2){
      phiinum = 1;
    }else{
      phiinum = (int) phionum / 4;
    }
    std::stringstream thetanumstr(argv[7]); // number of theta directions
    thetanumstr >> thetanum;
    std::cout<<"phiinum "<<phiinum<<" phionum "<<phionum<<" thetanum "<<thetanum<<std::endl;

    std::stringstream lambdanumstr(argv[8]); // number of wavelength
    lambdanumstr >> lambdanum;
    std::cout<<"lambdanum "<<lambdanum<<std::endl;

    std::stringstream etarestr(argv[9]); // index of refraction of the fiber (real part)
    etarestr >> etare;

    std::stringstream etaimstr(argv[10]); // index of refraction of the fiber (imaginary part)
    etaimstr >> etaim;
    std::cout<<"etare "<<etare<<" etaim "<<etaim<<std::endl;

    if (etaim!=0)
      lossless = false;

  }else{

    std::cout<<"------Simulating a non-elliptical cross-section------"<<std::endl;

    std::string xfile = argv[3]; // filename of x coordinates of nodes 
    std::string yfile = argv[4]; // filename of y coordinates of nodes
    std::cout<<"x coord file: "<<xfile<<"; y coord file: "<<yfile<<std::endl;
    nodes = readgeometry(xfile, yfile); // create nodes using the input files

    radius1 = -1; // setting the inital to be a negative number
    for (int i = 0; i < nodes.size(); ++i){
      double curradius = nodes[i].norm();
      if (radius1 < curradius)
        radius1 = curradius;
    }
    std::cout<<"radius "<<radius1<<std::endl;
    std::cout<<"number of elements "<<nodes.size()<<std::endl;

    std::stringstream phionumstr(argv[5]); // number of outgoing phi directions
    phionumstr >> phionum;
    phiinum = phionum;

    std::stringstream thetanumstr(argv[6]); // number of theta directions
    thetanumstr >> thetanum;
    std::cout<<"phiinum "<<phiinum<<" phionum "<<phionum<<" thetanum "<<thetanum<<std::endl;

    std::stringstream lambdanumstr(argv[7]); // number of wavelength
    lambdanumstr >> lambdanum;
    std::cout<<"lambdanum "<<lambdanum<<std::endl;

    std::stringstream etarestr(argv[8]); // index of refraction of the fiber (real part)
    etarestr >> etare;

    std::string iorfile = argv[9];
    kval.resize(lambdanum);
    readiorimag(kval, lambdanum, iorfile); // read in wavelength dependent imaginary part of ior

    lossless = false;
  }

  // wavelength parameter
  int lambdastart = 400;
  int lambdaend = 700;

  // read in hankel tables
  std::string hankeldir = "../hankeldouble/";
  filename = hankeldir + "hankelre.binary";
  std::ifstream myfile(filename, std::ios::in|std::ios::binary);
  myfile.read((char*)&hankelre[0], TABLE_SIZE * sizeof(double));

  filename = hankeldir + "hankelim.binary";
  std::ifstream myfile2(filename, std::ios::in|std::ios::binary);
  myfile2.read((char*)&hankelim[0], TABLE_SIZE * sizeof(double));

  filename = hankeldir + "dhankelre.binary";
  std::ifstream myfile3(filename, std::ios::in|std::ios::binary);
  myfile3.read((char*)&dhankelre[0], TABLE_SIZE * sizeof(double));

  filename = hankeldir + "dhankelim.binary";
  std::ifstream myfile4(filename, std::ios::in|std::ios::binary);
  myfile4.read((char*)&dhankelim[0], TABLE_SIZE * sizeof(double));

  int quadrature = 2;
  double mur = 1.0;
  double Dis = 1000*radius1;
  char mode = 'M';
  double freq;
  std::complex<double> eta, epsr;

  struct timeval start, end;
  gettimeofday(&start, NULL);

  int nb_samples = thetanum * phiinum * phionum;
  std::cout<<"nb_samples "<<nb_samples<<std::endl; 
  int nb_pdf = thetanum * phiinum;
  float unit = (float) 2*M_PI / phionum;

  Eigen::VectorXf vect(nb_samples);
  Eigen::VectorXf pdf(nb_samples);
  Eigen::VectorXf cdf(nb_samples);
  Eigen::VectorXf cstot(nb_pdf * lambdanum);

  Eigen::VectorXf sigma;
  double energy;
  float cs;

  char dir_array[output.length()];
  strcpy(dir_array, output.c_str());
  mkdir(dir_array, 0777);
  for (int i = 0; i < lambdanum; ++i){
    std::cout<<"lambdaindex "<<i<<std::endl;
    double lambda = (lambdastart + (double)(lambdaend-lambdastart)/(double)lambdanum * i)*1e-9;

    // index of refraction calculation
    if (ellipse)
      eta = etare - etaim * cunit;
    else
      eta = 1.55 - kval[i] * cunit;
    epsr = eta * eta;

    freq = 299792458.0/double(lambda);
    MoM_ob m1;
    double theta;
    Eigen::PartialPivLU<Eigen::MatrixXcd> luvar;
    #pragma omp parallel for private (m1, theta, luvar)
    for (int j = 0 ; j < thetanum; ++j){
      theta = M_PI / 2 - M_PI / 2 * (double) j / (double) thetanum;

      if (ellipse)
        m1 = MoM_ob(numel, radius1, radius2, quadrature, 0, mode, freq, mur, epsr, theta);
      else
        m1 = MoM_ob(nodes, quadrature, 0, mode, freq, mur, epsr, theta);

      m1.assembly();
      luvar.compute(m1.Z);
      for (int k = 0; k < phiinum; ++k){
        double phi_i;
        if (phiinum==360)
          phi_i = (double) k / (double) phiinum * M_PI * 2;
        else if (phiinum==90)
          phi_i = (double) k / (double) phiinum * M_PI / 2;
        else
          phi_i = 0;
        m1.updatewave(phi_i, 'M', freq, theta);
        m1.incident();
        m1.sol = luvar.solve(m1.vvec);
        Eigen::VectorXf sigma1(phionum), sigma2(phionum);
        m1.BSDF(phionum, Dis, sigma1);
        double energy1, energy2;
        if (lossless)
          energy1 = 0;
        else
          energy1= m1.energy();

        m1.updatewave(phi_i, 'E', freq, theta);
        m1.incident();
        m1.sol = luvar.solve(m1.vvec);
        m1.BSDF(phionum, Dis, sigma2);
        if (lossless)
          energy2 = 0;
        else
          energy2 = m1.energy();

        sigma = (sigma1 + sigma2) * 0.5f;
        cs = (float) sigma.sum() / phionum * 2.f * M_PI;
        energy = (energy1 + energy2) / 2;

        cstot(i*thetanum*phiinum + j*phiinum + k) = (float) (cs - energy);

        // process dist, cs, energy info for each wavelength and save to disk
        int vectindex = j*phiinum*phionum+k*phionum;
        postprocess(lossless, sigma, energy, phionum, vectindex, unit, vect, pdf, cdf);
      }
    }
    // write scattering distribution, pdf, cdf to disk
    filename = output + "TEM_"+std::to_string(i)+".binary";
    std::ofstream out(filename, std::ios::out|std::ios::binary|std::ios_base::app);
    out.write((char *) &vect(0), sizeof(float)*nb_samples);

    filename2 = output + "TEM_"+std::to_string(i)+"_pdf.binary";
    std::ofstream out2(filename2, std::ios::out|std::ios::binary|std::ios_base::app);
    out2.write((char *) &pdf(0), sizeof(float)*nb_samples);

    filename3 = output + "TEM_"+std::to_string(i)+"_cdf.binary";
    std::ofstream out3(filename3, std::ios::out|std::ios::binary|std::ios_base::app);
    out3.write((char *) &cdf(0), sizeof(float)*nb_samples);
  }
  // compute cross-section ratio and write to disc
  float bound = 6;
  float maxcs = std::min(bound, cstot.maxCoeff());
  std::cout<<"maxcs "<<maxcs<<std::endl;
  Eigen::VectorXf ratio = cstot / maxcs;
  // clamp ratio
  for (int i =0; i < nb_pdf * lambdanum; ++i){
    if (ratio(i) > 1)
      ratio(i) = 1;
    if(ratio(i) < 0){
      ratio(i) = 0;
      std::cout<<"ratio(i)<0, i "<<i<<" ratio(i) "<<ratio(i)<<std::endl;
    }
  }
  // output ratio
  filename = output + "ratio.binary";
  std::ofstream out(filename, std::ios::out|std::ios::binary|std::ios_base::app);
  out.write((char *) &ratio(0), sizeof(float)*nb_pdf*lambdanum);

  gettimeofday(&end, NULL);

  double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
                  end.tv_usec - start.tv_usec) / 1.e6;
  std::cout<<"time used: "<<delta<<std::endl;

  return 0;
}
