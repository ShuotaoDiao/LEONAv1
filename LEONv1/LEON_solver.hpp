//
//  LEON_solver.hpp
//  LEONv1
//
//  Created by Shuotao Diao on 4/6/22.
//

#ifndef LEON_solver_hpp
#define LEON_solver_hpp

#include <stdio.h>
#include <stdlib.h> // rand
#include <ctime>
#include <cmath>
#include <stdexcept>

#include "LEON_dataStructure.hpp"
#include "LEON_ioDB.hpp"
#include "LEON_ioModel.hpp"
#include "LEON_ioSto.hpp"

// obtain dual multiplers of the second stage, given x (first stage decision variable)
dualMultipliers twoStageLP_secondStageDual(const std::vector<double>& x, standardTwoStageParameters& model_parameters, const secondStageRHSpoint& rhs, const secondStageRHSmap& RHSmap);

// functions for generating feasibility cut
// find extreme ray
dualMultipliers twoStageLP_secondStageExtremRay(const std::vector<double>& x, standardTwoStageParameters& model_parameters, const secondStageRHSpoint& rhs, const secondStageRHSmap& RHSmap);
// construct feasibility cut
feasibilityCut twoStageLP_feasibilityCutGeneration(const dualMultipliers& extremeRay, standardTwoStageParameters& model_parameters, const secondStageRHSpoint& rhs, const secondStageRHSmap& RHSmap);

// projection for the first stage in the two stage linear programming
std::vector<double> twoStageLP_projection(const std::vector<double>& x, standardTwoStageParameters& model_parameters);


// LEON solver
std::vector<double> leon_solver(const std::string& folder_path,
            int max_iterates,
            const std::vector<double>& x_init,
            double Dx,
            int m,
            double max_subgradient);

#endif /* LEON_solver_hpp */
