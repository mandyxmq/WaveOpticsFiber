//
//  Basis.hpp
//  2Dsolver
//
//  Created by Mandy Xia.
//  Copyright © 2019 MandyXia. All rights reserved.
//

#ifndef Basis_hpp
#define Basis_hpp

#include <stdio.h>

class Basis{
public:
    std::vector<int> el;
    Basis():el({}){};
    Basis(int e1, int e2){
        el = {e1, e2};
    };
};

#endif /* Basis_hpp */
