//
//  main.cpp
//  LEONv1
//
//  Created by Shuotao Diao on 4/6/22.
//

#include <iostream>
#include <ilcplex/ilocplex.h>

#include "LEON_utils.hpp"

void qp_dual_test() {
    std::cout << "Test on basic solver performance!\n";
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x(env,2,-5,5,ILOFLOAT);
    mod.add(x);
    IloExpr expr_obj(env);
    expr_obj = x[0] * x[0] + x[1] * x[1];
    IloObjective obj = IloMinimize(env, expr_obj);
    mod.add(obj);
    // constrant
    IloRangeArray constraints(env);
    IloExpr expr_con(env);
    expr_con = x[0] + x[1];
    constraints.add(expr_con == 1);
    // extra constraints
    constraints.add(x[0] + x[1] <= 2);
    mod.add(constraints);
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.solve();
    // output solutions
    std::cout << "Optimal value: " << cplex.getObjValue() << std::endl;
    std::cout << "x[0] = " << cplex.getValue(x[0]) << std::endl;
    std::cout << "x[1] = " << cplex.getValue(x[1]) << std::endl;
    // dual variables
    IloNumArray duals(env);
    cplex.getDuals(duals, constraints);
    std::cout << "duals[0] = " << duals[0] << std::endl;
    std::cout << "duals[1] = " << duals[1] << std::endl;
    env.end();
}

void bk19() {
    std::string folderPath = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment11/replication6";
    std::string method = "EpanechnikovKernel2";
    std::string trainDB_path = folderPath + "/" + method; // path
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment11/trueDist";
    int m = 1;
    double max_subgradient = 190;
    double Dx = 40;
    int N = 50;
    std::vector<double> x_init(4,0.0);
    double Dx2 = 2.0 * Dx;
    int maxOuterLoops[] = {10,15};
    for (int idx = 0; idx < 1; ++idx) {
        interface_leon(folderPath, validationFolder_path, method, x_init, maxOuterLoops[idx], Dx2, max_subgradient, m, N);
    }
}

void bk19_2(int caseNumber, std::string method) {
    std::string folderPath = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment11/replication" + std::to_string(caseNumber);
    //std::string method = "EpanechnikovKernel2";
    std::string trainDB_path = folderPath + "/" + method; // path
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment11/trueDist";
    int m = 1;
    double max_subgradient = 190;
    double Dx = 40;
    int N = 50;
    std::vector<double> x_init(4,0.0);
    double Dx2 = 2.0 * Dx;
    int maxOuterLoops[] = {13};//{10,15};
    for (int idx = 0; idx < 1; ++idx) {
        interface_leon(folderPath, validationFolder_path, method, x_init, maxOuterLoops[idx], Dx2, max_subgradient, m, N);
    }
} // end bk_19_2

void baa99_large() {
    std::string folderPath = "/Users/sonny/Documents/numericalExperiment/baa99/experiment8/case5";
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/baa99/experiment8/kNNValidation2";
    std::string method = "kNN2";
    int m = 1;
    double max_subgradient = 1;
    double Dx = 1;
    int N = 50;
    std::vector<double> x_init(50,0.0);
    int maxOuterLoops[] = {12,20};
    for (int idx = 0; idx < 2; ++idx) {
        interface_leon(folderPath, validationFolder_path, method, x_init, maxOuterLoops[idx], Dx, max_subgradient, m, N);
    }
} // end baa99_large

void baa99_large2(int caseNumber, std::string method) {
    std::string folderPath = "/Users/sonny/Documents/numericalExperiment/baa99/experiment8/case" + std::to_string(caseNumber);
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/baa99/experiment8/kNNValidation2";
    int m = 1;
    double max_subgradient = 1;
    double Dx = 1;
    int N = 50;
    std::vector<double> x_init(50,0.0);
    int maxOuterLoops[] = {12,20};
    for (int idx = 0; idx < 2; ++idx) {
        interface_leon(folderPath, validationFolder_path, method, x_init, maxOuterLoops[idx], Dx, max_subgradient, m, N);
    }
} // end baa99_large

int main(int argc, const char * argv[]) {
    //bk19();
    std::string method = "GaussianKernel2"; //naiveKernel2, EpanechnikovKernel2, quarticKernel2, GaussianKernel2
    for (int caseNumber = 6; caseNumber < 16; ++caseNumber) {
        bk19_2(caseNumber, method);
    }
    //baa99_large();
    /*
    std::string method = "naiveKernel";
    for (int caseNumber = 1; caseNumber < 11; ++caseNumber) {
        baa99_large2(caseNumber, method);
    }
     */
    return 0;
}
