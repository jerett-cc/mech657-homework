#ifndef SpatialLinearSystemCons_H
#define SpatialLinearSystemCons_H


#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"
#include "../Eigen/Sparse"
#include "QuasiEuler.h"
#include "Geometry.h"

// needs to calculate the matrix, calculate the RHS, and solve the system for deltaQ

class Solver{

  public:
    std::vector<Eigen::Vector3d> delta_Q;
    void setupSystem();//does this need to include a parameter Data Class?
    void solveSystem();//updates delta_Q with the result of the system solution.

  private:
    Eigen::MatrixXd A;
    Eigen::VectorXd b;

    void setup_A();
    void setup_b();

};



























#endif
