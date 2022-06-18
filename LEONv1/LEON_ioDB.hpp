//
//  LEON_ioDB.hpp
//  LEONv1
//
//  Created by Shuotao Diao on 4/6/22.
//

#ifndef LEON_ioDB_hpp
#define LEON_ioDB_hpp

#include <stdio.h>

#include <fstream>
#include <sstream>

#include "LEON_dataStructure.hpp"

std::vector<std::vector<dataPoint>> readNonparametricDB(std::string readPath); // read database from a text file

void printNonparametricDB(const std::vector<std::vector<dataPoint>>& dataPointDB);
// print dataPoint
void printDataPoint(const dataPoint& dataPoint01);
void inputDBTest(); // test on input functions of nonparametric DB

#endif /* LEON_ioDB_hpp */
