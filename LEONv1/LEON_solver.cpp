//
//  LEON_solver.cpp
//  LEONv1
//
//  Created by Shuotao Diao on 4/6/22.
//

#include "LEON_solver.hpp"

// declare global variables (need to be defined in the source file)
double SOLVER_PRECISION_LOWER = -1e-6;
double SOLVER_PRECISION_UPPER = 1e-6;

double SOLVER_INF = 1e10;
double DUAL_PRECISION_LOWER = -1e-2;
double DUAL_PRECEISION_UPPER = 1e-2;

// compare duals
bool if_vec_equal(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        return false;
    }
    // use L1 distance to measure the distance between two vectors
    double diff = 0;
    for (int index = 0; index < vec1.size(); ++index) {
        diff += std::abs(vec1[index] - vec2[index]);
        if (diff > DUAL_PRECEISION_UPPER) { // difference is significant
            //std::cout << "Debug check if two vector are equal.\n";
            //std::cout << "index: " << index;
            //std::cout << "  diff: " << diff << std::endl;
            return false;
        }
    }
     
    // coarse way
    /*
    for (int index = 0; index < vec1.size(); ++index) {
        double diff = std::abs(vec1[index] - vec2[index]);
        if (diff > DUAL_PRECEISION_UPPER) { // difference is significant
            std::cout << "Debug check if two vector are equal.\n";
            std::cout << "index: " << index;
            std::cout << "  diff: " << diff << std::endl;
            std::cout << "vec1[index] = " << vec1[index] << std::endl;
            std::cout << "vec2[index] = " << vec2[index] << std::endl;
            return false;
        }
    }
    */
    return true;
}

bool if_duals_equal(const dualMultipliers& dual1, const dualMultipliers& dual2){
    if (if_vec_equal(dual1.equality, dual2.equality)) {
        if (if_vec_equal(dual1.inequality, dual2.inequality)) {
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
    return true;
}

bool if_new_dual(const std::vector<dualMultipliers>& duals, const std::vector<int>& indices, const dualMultipliers& candidate_dual) {
    if (indices.size() < 1) {
        return true;
    }
    for (int index = 0; index < indices.size(); ++index) {
        if (if_duals_equal(duals[indices[index]], candidate_dual)) { // if candidate_dual is equal to some dual in the duals
            //std::cout << "Two duals are equal.\n";
            //print(duals[index].inequality);
            //print(candidate_dual.inequality);
            return false;
        }
    }
    return true;
}

bool if_new_dual(const std::vector<dualMultipliers>& duals, const dualMultipliers& candidate_dual) {
    if (duals.size() < 1) {
        return true; // new dual is found
    }
    for (int index = 0; index < duals.size(); ++index) {
        if (if_duals_equal(duals[index], candidate_dual)) {
            std::cout << "Dual is not new.\n";
            return false; // candidate dual already exists
        }
    }
    return true; // new dual is found
}


// obtain dual multiplers of the second stage by solving the primal, given x (first stage decision variable)
dualMultipliers twoStageLP_secondStageDual(const std::vector<double>& x, standardTwoStageParameters& model_parameters, const secondStageRHSpoint& rhs, const secondStageRHSmap& RHSmap) {
    dualMultipliers pi; // initialize dual multipliers
    long y_size = model_parameters.D.num_col;
    // set up the model
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray y(env,y_size,0,IloInfinity,ILOFLOAT); // second stage decision variables
    mod.add(y);
    // objective function
    IloExpr expr_obj(env);
    for (auto it = model_parameters.d.vec.begin(); it != model_parameters.d.vec.end(); ++it) {
        expr_obj += (it -> second) * y[it -> first];
    }
    IloObjective obj = IloMinimize(env,expr_obj);
    mod.add(obj);
    // equality constraints Dy + Cx = e
    IloRangeArray constraintsEquality(env);
    std::vector<IloExpr> exprs_eq;
    for (int index_eq = 0; index_eq < model_parameters.D.num_row; ++index_eq) {
        IloExpr expr(env);
        exprs_eq.push_back(expr);
    }
    // coefficients before y
    for (std::map<std::pair<int,int>, double>::iterator it = model_parameters.D.mat.begin(); it != model_parameters.D.mat.end(); ++it) {
        // get the location of the entry
        exprs_eq[(it -> first).first] += (it -> second) * y[(it -> first).second];
    }
    // coefficients before x (deterministic part)
    for (std::map<std::pair<int,int>, double>::iterator it = model_parameters.C.mat.begin(); it != model_parameters.C.mat.end(); ++it) {
        // get the location of the entry
        exprs_eq[(it -> first).first] += (it -> second) * x[(it -> first).second];
    }
    // right hand side (deterministic part)
    for (auto it = model_parameters.e.vec.begin(); it != model_parameters.e.vec.end(); ++it) {
        exprs_eq[it -> first] -= (it -> second);
    }
    // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
    for (int idx_Ce = 0; idx_Ce < rhs.Ce.size(); ++idx_Ce) {
        exprs_eq[RHSmap.Ce_map[idx_Ce].first] += rhs.Ce[idx_Ce] * x[RHSmap.Ce_map[idx_Ce].second];
    }
    // coefficients before x (stochastic part) inequality (location is behind equality constraints)
    for (int idx_Ci = 0; idx_Ci < rhs.Ci.size(); ++idx_Ci) {
        exprs_eq[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq] += rhs.Ci[idx_Ci] * x[RHSmap.Ci_map[idx_Ci].second];
    }
    // right hand side (stochastic part) equality be_(i) equality
    for (int idx_be = 0; idx_be < rhs.be.size(); ++idx_be) {
        exprs_eq[RHSmap.be_map[idx_be]] -= rhs.be[idx_be];
    }
    // right hand side (stochastic part) equality bi_(i) inequality
    for (int idx_bi = 0; idx_bi < rhs.bi.size(); ++idx_bi) {
        exprs_eq[RHSmap.bi_map[idx_bi] + model_parameters.num_eq] -= rhs.bi[idx_bi];
    }
    // add the equality constraints
    for (int index_eq = 0; index_eq < model_parameters.D.num_row; ++index_eq) {
        constraintsEquality.add(exprs_eq[index_eq] == 0);
    }
    mod.add(constraintsEquality);
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    IloBool solvable_flag = cplex.solve();
    IloNumArray dual_equality(env);
    if (solvable_flag == IloTrue) {
        pi.feasible_flag = true; // tell the subproblem is feasible for a given x, first stage decision variable
        cplex.getDuals(dual_equality,constraintsEquality);
        for (int index_eq = 0; index_eq < model_parameters.D.num_row; ++index_eq) {
            double pi_temp = dual_equality[index_eq]; // move y to the right hand side
            pi.equality.push_back(pi_temp);
        }
    }
    else {
        pi.feasible_flag = false; // tell the subproblem is infeasible for given x
    }
    env.end();
    return pi;
}


// functions for generating feasibility cut
dualMultipliers twoStageLP_secondStageExtremRay(const std::vector<double>& x, standardTwoStageParameters& model_parameters, const secondStageRHSpoint& rhs, const secondStageRHSmap& RHSmap) {
    dualMultipliers extremeRay;
    extremeRay.feasible_flag = false;
    // obtain the sizes of input parameters
    long y_size = model_parameters.D.num_col;
    // set up the model
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray y(env,y_size,0,IloInfinity,ILOFLOAT); // second stage decision variables
    mod.add(y);
    // objective function
    IloExpr expr_obj(env);
    for (auto it = model_parameters.d.vec.begin(); it != model_parameters.d.vec.end(); ++it) {
        expr_obj += (it -> second) * y[it -> first];
    }
    IloObjective obj = IloMinimize(env,expr_obj);
    mod.add(obj);
    // equality constraints
    IloRangeArray constraintsEquality(env);
    std::vector<IloExpr> exprs_eq;
    for (int index_eq = 0; index_eq < model_parameters.D.num_row; ++index_eq) {
        IloExpr expr(env);
        exprs_eq.push_back(expr);
    }
    // coefficients before y
    for (std::map<std::pair<int,int>, double>::iterator it = model_parameters.D.mat.begin(); it != model_parameters.D.mat.end(); ++it) {
        // get the location of the entry
        exprs_eq[(it -> first).first] += (it -> second) * y[(it -> first).second];
    }
    // coefficients before x (deterministic part)
    for (std::map<std::pair<int,int>, double>::iterator it = model_parameters.C.mat.begin(); it != model_parameters.C.mat.end(); ++it) {
        // get the location of the entry
        exprs_eq[(it -> first).first] += (it -> second) * x[(it -> first).second];
    }
    // right hand side (deterministic part)
    for (auto it = model_parameters.e.vec.begin(); it != model_parameters.e.vec.end(); ++it) {
        exprs_eq[it -> first] -= (it -> second);
    }
    // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
    for (int idx_Ce = 0; idx_Ce < rhs.Ce.size(); ++idx_Ce) {
        exprs_eq[RHSmap.Ce_map[idx_Ce].first] += rhs.Ce[idx_Ce] * x[RHSmap.Ce_map[idx_Ce].second];
    }
    // coefficients before x (stochastic part) inequality (location is behind equality constraints)
    for (int idx_Ci = 0; idx_Ci < rhs.Ci.size(); ++idx_Ci) {
        exprs_eq[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq] += rhs.Ci[idx_Ci] * x[RHSmap.Ci_map[idx_Ci].second];
    }
    // right hand sie (stochastic part) equality be_(i) equality
    for (int idx_be = 0; idx_be < rhs.be.size(); ++idx_be) {
        exprs_eq[RHSmap.be_map[idx_be]] -= rhs.be[idx_be];
    }
    // right hand sie (stochastic part) equality bi_(i) inequality
    for (int idx_bi = 0; idx_bi < rhs.bi.size(); ++idx_bi) {
        exprs_eq[RHSmap.bi_map[idx_bi] + model_parameters.num_eq] -= rhs.bi[idx_bi];
    }
    // add the equality constraints
    for (int index_eq = 0; index_eq < model_parameters.D.num_row; ++index_eq) {
        constraintsEquality.add(exprs_eq[index_eq] == 0);
    }
    mod.add(constraintsEquality);
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setParam(IloCplex::PreInd, false); // need to turn off presolve in order to get dual extreme rays
    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // use dual simplex optimizer
    cplex.solve(); // solve the problem
    IloNumArray extremeRay_eq(env);
    cplex.dualFarkas(constraintsEquality, extremeRay_eq);
    for (int index_eq = 0; index_eq < model_parameters.D.num_row; ++index_eq) {
        double pi_temp = extremeRay_eq[index_eq]; // move y to the right hand side
        extremeRay.equality.push_back(pi_temp);
    }
    return extremeRay;
}


// construct feasibility cut
feasibilityCut twoStageLP_feasibilityCutGeneration(const dualMultipliers& extremeRay, standardTwoStageParameters& model_parameters, const secondStageRHSpoint& rhs, const secondStageRHSmap& RHSmap) {
    feasibilityCut cut_scenario;
    // intercept deterministic part  A_newRow x \leq b_newRow, so the negative is taken
    cut_scenario.b_newRow = (-1.0) * (model_parameters.e * extremeRay.equality);
    // stochastic part
    // equality part
    for (int idx_eq = 0; idx_eq < rhs.be.size(); ++idx_eq) {
        cut_scenario.b_newRow -= extremeRay.equality[RHSmap.be_map[idx_eq]] * rhs.be[idx_eq];
    }
    // inequality part (before standardizing) inequality constraint is after the equality constraints
    for (int idx_ineq = 0; idx_ineq < rhs.bi.size(); ++idx_ineq) {
        cut_scenario.b_newRow -= extremeRay.equality[RHSmap.bi_map[idx_ineq] + model_parameters.num_eq] * rhs.bi[idx_ineq];
    }
    // slope
    // deterministic part
    cut_scenario.A_newRow = (-1.0) * (extremeRay.equality * model_parameters.C);
    // stochastic part
    // equality
    for (int idx_Ce = 0; idx_Ce < rhs.Ce.size(); ++idx_Ce) {
        cut_scenario.A_newRow[RHSmap.Ce_map[idx_Ce].second] += -1.0 * rhs.Ce[idx_Ce] * extremeRay.equality[RHSmap.Ce_map[idx_Ce].first];
    }
    // inequality before standardizing
    for (int idx_Ci = 0; idx_Ci < rhs.Ci.size(); ++idx_Ci) {
        cut_scenario.A_newRow[RHSmap.Ce_map[idx_Ci].second] += -1.0 * rhs.Ci[idx_Ci] * extremeRay.equality[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq];
    }
    return cut_scenario;
}

// projection for the first stage in the two stage linear programming
std::vector<double> twoStageLP_projection(const std::vector<double>& x, standardTwoStageParameters& model_parameters) {
    std::vector<double> x_projected(x.size(),0.0);
    // solve a quadratic programming
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x_temp(env,x.size(),-IloInfinity,IloInfinity,ILOFLOAT);
    mod.add(x_temp);
    IloExpr expr_obj(env);
    for (int x_index = 0; x_index < x.size(); ++x_index) {
        expr_obj += x_temp[x_index] * x_temp[x_index] - 2.0 * x_temp[x_index] * x[x_index];
    }
    IloObjective obj = IloMinimize(env,expr_obj); // objective function
    mod.add(obj);
    // constraints
    std::vector<IloExpr> exprs;
    for (int index_cons = 0; index_cons < model_parameters.A.num_row; ++index_cons) {
        IloExpr expr(env);
        exprs.push_back(expr);
    }
    for (std::map<std::pair<int,int>, double>::iterator it = model_parameters.A.mat.begin(); it != model_parameters.A.mat.end(); ++it) {
        // get the location of the entry
        exprs[(it -> first).first] += (it -> second) * x_temp[(it -> first).second];
    }
    // right hand side
    for (auto it = model_parameters.b.vec.begin(); it != model_parameters.b.vec.end(); ++it) {
        exprs[it -> first] -= (it -> second);
    }
    // add constraints
    for (int index_cons = 0; index_cons < model_parameters.A.num_row; ++index_cons) {
        mod.add(exprs[index_cons] <= 0);
    }
    // create cplex environment
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    cplex.solve();
    // obtain the projected point
    for (int x_index = 0; x_index < model_parameters.A.num_col; ++x_index) {
        x_projected[x_index] = cplex.getValue(x_temp[x_index]);
        //std::cout << cplex.getValue(x_temp[x_index]) << std::endl;
    }
    env.end();
    return x_projected;
}


// LEON solver
std::vector<double> leon_solver(const std::string& folder_path,
            int max_iterates,
            const std::vector<double>& x_init,
            double Dx,
            int m,
            double max_subgradient) {
    // INITIALIZATION
    bool flag_be; // tell if be stochastic is generated
    bool flag_bi; // tell if bi stochastic is generated
    bool flag_Ce; // tell if Ce stochastic is generated
    bool flag_Ci; // tell if Ci stochastic is generated
    std::vector<secondStageRHSpoint> RHS_dataset;
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string Ce_DB_path = folder_path + "/Ce_DB.txt";
    std::string Ci_DB_path = folder_path + "/Ci_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(leon).txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* Ce_DB_path_const = Ce_DB_path.c_str();
    const char* Ci_DB_path_const = Ci_DB_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    std::ifstream readFile_Ce(Ce_DB_path_const);
    std::ifstream readFile_Ci(Ci_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    std::vector<std::vector<dataPoint>> Ce_DB;
    std::vector<std::vector<dataPoint>> Ci_DB;
    
    // create model structure
    standardTwoStageParameters model_parameters;
    // create sto object
    secondStageRHSmap RHSmap;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
        flag_be = true;
    }
    else {
        readFile_be.close(); // close the file
        flag_be = false;
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
        flag_bi = true;
    }
    else {
        readFile_bi.close(); // close the file
        flag_bi = false;
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read Ce
    if (readFile_Ce.is_open()) {
        std::cout << "Ce_DB stochastic part is found." << std::endl;
        readFile_Ce.close(); // close the file
        // Ce database
        Ce_DB = readNonparametricDB(Ce_DB_path);
        flag_Ce = true;
    }
    else {
        readFile_Ce.close(); // close the file
        flag_Ce = false;
        std::cout << "Ce_DB is not found!" << std::endl;
    }
    // read Ci
    if (readFile_Ci.is_open()) {
        std::cout << "Ci_DB stochastic part is found." << std::endl;
        readFile_Ci.close(); // close the file
        // Ci database
        Ci_DB = readNonparametricDB(Ci_DB_path);
        flag_Ci = true;
    }
    else {
        readFile_Ci.close(); // close the file
        flag_Ci = false;
        std::cout << "Ci_DB is not found!" << std::endl;
    }
    
    // read model file
    model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    RHSmap = readStochasticMap(sto_path);
    // initialize feasibility cut collection
    std::vector<feasibilityCut> feasibility_cuts;
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write initial setup
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "LEON (v1) is initialized\n";
    writeFile << "LEON (v1) is initialized\n";
    std::cout << "Algorithmic Parameters\n";
    writeFile << "Algorithmic Parameters\n";
    std::cout << "Dx, m, max_subgradient\n";
    writeFile << "Dx, m, max_subgradient\n";
    std::cout << Dx << ", " << m << ", " << max_subgradient << std::endl;
    writeFile << Dx << ", " << m << ", " << max_subgradient << std::endl;
    writeFile << "Initial solution: ";
    long x_init_size = x_init.size();
    for (int index = 0; index < x_init_size - 1; ++index) {
        writeFile << x_init[index] << ", ";
    }
    writeFile << x_init[x_init_size - 1] << "\n";
    writeFile << "Max number of iterations: " << max_iterates << "\n";
    std::cout << "Problem Complexity\n";
    writeFile << "Problem Complexity\n";
    std::cout << "A_num_row, A_num_col\n";
    writeFile << "A_num_row, A_num_col\n";
    std::cout << model_parameters.A.num_row << ", " << model_parameters.A.num_col << std::endl;
    writeFile << model_parameters.A.num_row << ", " << model_parameters.A.num_col << std::endl;
    std::cout << "D_num_row, D_num_col (after converting into standard form)\n";
    writeFile << "D_num_row, D_num_col (after converting into standard form)\n";
    std::cout << model_parameters.D.num_row << ", " << model_parameters.D.num_col << std::endl;
    writeFile << model_parameters.D.num_row << ", " << model_parameters.D.num_col << std::endl;
    // set up initial solution
    long x_size = model_parameters.c.num_entry;
    //long A_rowsize = model_parameters.A.num_row;
    //long A_colsize = model_parameters.A.num_col;
    // initial variable
    std::vector<double> x_new(x_size,0.0);
    std::vector<double> x_leon(x_size,0.0);
    // projection
    std::vector<double> x_old = twoStageLP_projection(x_init, model_parameters);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        x_new[x_index] = x_old[x_index];
        x_leon[x_index] = x_old[x_index];
    }
    // main loop
    // write main loop
    writeFile << "Main loop\n";
    // iterate index
    int dataset_idx= 0;
    // main loop (outer loop)
    for (int outer_it = 0; outer_it < max_iterates; ++outer_it) {
        // initialization reset x_leon
        for (int x_index = 0; x_index < x_size; ++x_index) {
            x_leon[x_index] = 0;
        }
        // max number iterates in the inner loop
        int maxInner_iterates = (outer_it + 1) * m;
        std::cout << "Outer Loop: Iteration " << outer_it + 1 << std::endl;
        writeFile << "##################################################\n";
        writeFile << "Outer Loop: Iteration " << outer_it + 1 << "\n";
        // step size
        double stepsize = Dx / (max_subgradient * sqrt((double) maxInner_iterates));
        for (int inner_it = 0; inner_it < maxInner_iterates; ++inner_it) {
            // compute quasigradient
            feasibilityCut feasibilityCut_scenario;
            bool feasibility_flag = true;
            std::vector<double> pi_AGG(model_parameters.D.num_row,0.0); // aggregate sum of duals
            long dataset_size = 0;
            if (flag_bi == true) {
                dataset_size = bi_DB[dataset_idx].size();
            }
            else if (flag_be == true) {
                dataset_size = be_DB[dataset_idx].size();
            }
            else if (flag_Ci == true) {
                dataset_size = Ci_DB[dataset_idx].size();
            }
            else if (flag_Ce == true) {
                dataset_size = Ce_DB[dataset_idx].size();
            }
            else {
                std::invalid_argument("ERROR: Database is missing.");
            }
            // write number of scenarios
            writeFile << "------------------------------------------------\n";
            writeFile << "Inner loop: Iteration  " << inner_it + 1 << "\n";
            writeFile << "Number of scenarios: " << dataset_size << "\n";
            std::vector<double> nsqg(x_size,0.0); // NSQG
            // loop through all scenarios in this dataset
            for (int scenario_index = 0; scenario_index < dataset_size; ++scenario_index) {
                std::cout << "scenario_index: " << scenario_index << std::endl;
                dataPoint bi_datapoint;
                if (flag_bi == true) {
                    bi_datapoint = bi_DB[dataset_idx][scenario_index];
                }
                dataPoint be_datapoint;
                if (flag_be == true) {
                    be_datapoint = be_DB[dataset_idx][scenario_index];
                }
                dataPoint Ci_datapoint;
                if (flag_Ci == true) {
                    Ci_datapoint = Ci_DB[dataset_idx][scenario_index];
                }
                dataPoint Ce_datapoint;
                if (flag_Ce == true) {
                    Ce_datapoint = Ce_DB[dataset_idx][scenario_index];
                }
                secondStageRHSpoint RHS_datapoint = merge_randomVector(be_datapoint, bi_datapoint, Ce_datapoint, Ci_datapoint);
                
                dualMultipliers dualsTemp = twoStageLP_secondStageDual(x_old, model_parameters, RHS_datapoint, RHSmap); // calculate duals
                
                if (dualsTemp.feasible_flag == false) { // if second stage subproblem is infeasible
                    // need to generate feasibility cut
                    // obtain extrem ray of the subproblem
                    // write one second stage problem is infeasible
                    writeFile << "Scenario " << scenario_index << " is infeasible, generate extreme ray\n";
                    std::invalid_argument("ERROR: Subproblem is infeasible.");
                }
                else {
                    // determinictic C
                    std::vector<double> pi_C = dualsTemp.equality * model_parameters.C;
                    std::vector<double> minus_pi_C = (-1.0) * pi_C;
                    // stochastic C
                    // equality
                    for (int idx_Ce = 0; idx_Ce < RHS_datapoint.Ce.size(); ++idx_Ce) {
                        //pi_C[RHSmap.Ce_map[idx_Ce].second] += RHS_dataset[kNNSet[index]].Ce[idx_Ce] * explored_duals[dual_index].equality[RHSmap.Ce_map[idx_Ce].first];
                        minus_pi_C[RHSmap.Ce_map[idx_Ce].second] -= RHS_datapoint.Ce[idx_Ce] * dualsTemp.equality[RHSmap.Ce_map[idx_Ce].first];
                    }
                    // inequality before standardizing
                    for (int idx_Ci = 0; idx_Ci < RHS_datapoint.Ci.size(); ++idx_Ci) {
                        //pi_C[RHSmap.Ci_map[idx_Ci].second] += RHS_dataset[kNNSet[index]].Ci[idx_Ci] * explored_duals[dual_index].equality[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq];
                        minus_pi_C[RHSmap.Ci_map[idx_Ci].second] -= RHS_datapoint.Ci[idx_Ci] * dualsTemp.equality[RHSmap.Ci_map[idx_Ci].first + model_parameters.num_eq];
                    }
                    nsqg = nsqg + RHS_datapoint.weight * minus_pi_C;
                } // end checking feasibility
            } // end reading all the scenarios in the dataset
            // first stage cost
            for (auto it = model_parameters.c.vec.begin(); it != model_parameters.c.vec.end(); ++it) {
                nsqg[it -> first] += (it -> second);
            }
            // write stepzie
            writeFile << "Stepsize: " << stepsize << "\n";
            // write quasigradient
            writeFile << "NSQG: ";
            for (int Gx_index = 0; Gx_index < x_size - 1; ++Gx_index) {
                writeFile << nsqg[Gx_index] << ", ";
            }
            writeFile << nsqg[x_size - 1] << "\n";
            // update estimate
            for (int x_index = 0; x_index < x_size; ++x_index) {
                x_old[x_index] = x_old[x_index] - nsqg[x_index] * stepsize;
            }
            // projection
            x_new = twoStageLP_projection(x_old,model_parameters);
            std::cout << "Inner Iteration: " << inner_it + 1 << std::endl;
            std::cout << "x: (";
            // update x_old and go to the next iteration
            for (int x_index = 0; x_index < x_size; ++x_index) {
                x_old[x_index] = x_new[x_index];
                std::cout << x_new[x_index] << " ";
            }
            std::cout << ")" <<std::endl;
            std::cout << "**************************" << std::endl;
            writeFile << "New estimated solution in the subsequence: ";
            for (int x_index = 0; x_index < x_size - 1; ++x_index) {
                writeFile << x_new[x_index] << ", ";
                // calculate x_robust (take the average of the subsequence of estimates)
                x_leon[x_index] +=  1.0 / ((double) maxInner_iterates) * x_new[x_index];
            }
            writeFile << x_new[x_size - 1] << "\n";
            x_leon[x_size - 1] +=  1.0 / ((double) maxInner_iterates) * x_new[x_size - 1];
            writeFile << "------------------------------------------------\n";
            dataset_idx++; // increment dataset index
        } // end inner loop
        // write new estimate
        std::cout <<  "New estimated solution by LEON: ";
        writeFile << "New estimated solution by LEON: ";
        for (int x_index = 0; x_index < x_size - 1; ++x_index) {
            std::cout << x_leon[x_index] << ", ";
            writeFile << x_leon[x_index] << ", ";
        }
        std::cout << x_leon[x_size - 1] << "\n";
        writeFile << x_leon[x_size - 1] << "\n";
        writeFile << "##################################################\n";
    } // end outer loop
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close(); // close writeFile
    return x_leon;
} // end LEON solver
