//
//  LEON_utils.hpp
//  LEONv1
//
//  Created by Shuotao Diao on 4/6/22.
//

#ifndef LEON_utils_hpp
#define LEON_utils_hpp

#include <stdio.h>

#include "LEON_solver.hpp"

double twoStageLP_secondStageCost(const std::vector<double>& x, standardTwoStageParameters& model_parameters, const secondStageRHSpoint& rhs, const secondStageRHSmap& RHSmap);

validationResult twoStageLP_validation_outputResultsV2(const std::string& folder_path, const std::vector<double>& x_candidate);

void twoStageLP_empirical_cost(const std::string& folder_path);

void interface_leon(const std::string& folder_path,
                    const std::string& validation_folder_path,
                    const std::string& method,
                    const std::vector<double>& x_init,
                    int max_outer_it,
                    double Dx,
                    double max_subgradient,
                    int m,
                    int N);

#endif /* LEON_utils_hpp */
