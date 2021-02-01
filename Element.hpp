//
//  Element.hpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//

#ifndef Element_hpp
#define Element_hpp

#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>

class Element{
    
public:    
  float length;
  std::vector<int> basis;
  Eigen::Vector3d tanvec;         // tangent vector is a 3d vector
  std::vector<Eigen::Vector3d> quadpoints;
  std::vector<Eigen::Vector3d> normal;         // each normal is a 3d vector

  int node1;
  int node2;
  Element(){};
  Element(int no1, int no2, Eigen::Vector3d n1, Eigen::Vector3d n2, std::vector<int> b, std::vector<double> xvec){
    node1 = no1;
    node2 = no2;
    length = (n1-n2).norm();
    basis = b;
    tanvec = (n2-n1)/(n1-n2).norm();
    // use xvec to generate quadpoints and normal
    for (int i = 0; i < xvec.size(); ++i){
      Eigen::Vector3d temp = n1+(1+xvec[i])/2*(n2-n1);
      quadpoints.push_back(temp);
      normal.push_back(Eigen::Vector3d{tanvec[1],-tanvec[0],0});
    }
  };
  ~Element(){};
    
  void print_qp() const{
    for (int i = 0; i <(this->quadpoints).size();++i){
      std::cout<<(this->quadpoints)[i]<<std::endl;
    }
  };
  void print_normal() const{
    for (int i = 0; i <normal.size();++i){
      std::cout<<normal[i]<<std::endl;
    }
  };
  void print_b() const{
    std::cout<<"basis 1"<<basis[0]<<std::endl;
    std::cout<<"basis 2"<<basis[1]<<std::endl;
  };
    
  // getters
  int getn1() const { return node1; };
  int getn2() const { return node2; };
  float getlen() const { return length; };
  int getbasis(int x) const { return basis[x]; };
  Eigen::Vector3d gettan() const { return tanvec; };
  Eigen::Vector3d getqp(int x) const { return quadpoints[x]; };
  Eigen::Vector3d getnormal(int x) const { return normal[x]; };
    
};

#endif /* Element_hpp */
