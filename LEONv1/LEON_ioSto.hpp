//
//  LEON_ioSto.hpp
//  LEONv1
//
//  Created by Shuotao Diao on 4/6/22.
//

#ifndef LEON_ioSto_hpp
#define LEON_ioSto_hpp

#include <stdio.h>
#include "LEON_dataStructure.hpp"

// including be, bi, Ce and Ci
secondStageRHSmap readStochasticMap(const std::string& stochasticPath);

// merge randomVector
secondStageRHSpoint merge_randomVector(const dataPoint& be_point, const dataPoint& bi_point, const dataPoint& Ce_point, const dataPoint& Ci_point);

#endif /* LEON_ioSto_hpp */
