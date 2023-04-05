#ifndef SpatialLinearSystemCons_H
#define SpatialLinearSystemCons_H

#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"
#include "../Eigen/Sparse"
#include "QuasiEuler.h"
#include "Geometry.h"
#include <iomanip>
#include <cmath>

/**
 * this class should take the data structure, the problem (which has the pressure sensor and
 * dissipation calculation) and produce a matrix A, the rhs b, and can solve it with the iterative method
 * and the diagonal form
 */
class Solver{
  public:
    std::vector<Eigen::Vector3d> delta_Q;//solution vector
    void setupSystem();
    void solveSystem();
    Solver(ProblemData *data, QuasiEuler *problem, double t_end = 0)
      : data(data), problem(problem)
      , stencil_rows(data->q.size()*3)
      , stencil_cols(3*(data->q.size()+4))
    {
      delta_Q = data->q;//we will overwrite the values in this. when we compute delta_Q
      //this could introduce bugs. be aware?
      A = Eigen::MatrixXd::Identity(data->getQVect().size(), data->getQVect().size());
      b = Eigen::VectorXd::Zero(data->getQVect().size(), 1);
      assert(A.cols() == data->q.size() * 3);//FIXME: will this assertion be tru in more than one dim?
      end_time = t_end;
    }
    double L2Error();
//  private: //FIXME: make some of these private??
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    int stencil_cols;
    int stencil_rows;
    double end_time;
    ProblemData *data;
    QuasiEuler *problem;

    void setup_A();
    void calculateAndAddSpatialMatrix();
    void calculateAndAddL();
    void calculateAndAddS();
    void reinit();

    void setup_b();
    void calcDe();
    void calcDx();
    void calcSx();
};

/**
 * returns the little el2 norm of deltaQ the change in Q for the given time step
 */
double Solver::L2Error(){
  double tmp = 0;
  for (int i = 0; i < delta_Q.size(); ++i)
  {
    tmp += delta_Q[i].squaredNorm();
  }
  return std::sqrt(tmp);
}

void Solver::setupSystem(){
  problem->calculateSensorContributions();
  setup_A();
  setup_b();
}

/**
 *this function should solve for delta Q, update Q by adding it to Q, and
 *then reinit the matrix A and the vector b to their initial states, so that
 *we can recompute in a next iteration*/
void Solver::solveSystem(){
//  std::cout << "A is " << A << std::endl;
//  std::cout << "b is " << b << std::endl;
  Eigen::VectorXd delta_q = A.colPivHouseholderQr().solve(b);
//  std::cout << "step solution is \n" << delta_q << std::endl;
  for (int i=0; i<delta_Q.size(); ++i)
  {
    for (int j=0; j<3; ++j)
    {
      delta_Q[i](j) = delta_q(3*i+j);
    }
  }
  *data += delta_Q;
  reinit();
}

/**
 *this function calculates the contribution of En and Dx and adds it to b */
void Solver::setup_b(){
  calcDe();
  calcDx();
  calcSx();
}

/**
 *this function calculates An and L and adds them to A. */
void Solver::setup_A(){
  calculateAndAddSpatialMatrix();
  calculateAndAddL();
  calculateAndAddS();
}

/**
 * this function calculates all the local flux jacobians and does an second order approx
 * to the total flux jacobian in space. It stitches each together in the system
 * matrix A.
 */
void Solver::calculateAndAddSpatialMatrix(){
  int local_matrix_size = data->q[0].size();
  Eigen::MatrixXd Ai(local_matrix_size, local_matrix_size);
  Eigen::MatrixXd Ai_n(local_matrix_size, local_matrix_size);
  Eigen::MatrixXd TOSHOW = Eigen::MatrixXd::Zero(data->getQVect().size(), data->getQVect().size());
  for (int i = 0; i < data->q.size()-1; ++i)//loop through all nodes but last one
  {
    Ai = problem->calculateLocalFluxJacobian(i);//calculate flux at node i
    Ai_n = problem->calculateLocalFluxJacobian(i+1);//calculate flux at node i+1
    int sub_index = local_matrix_size * (i);
    assert(sub_index+local_matrix_size + 3<=A.cols() && "indexing outside of bounds in the matrix");
    TOSHOW.block<3,3>(sub_index, sub_index+local_matrix_size) += 0.5/data->dx *Ai_n;
    TOSHOW.block<3,3>(sub_index + local_matrix_size, sub_index) += -0.5/data->dx*Ai;
  }
  //std::cout << "dxA" << std::endl << TOSHOW << std::endl;
  A =  A + problem->dt*TOSHOW;
}

/**
 * this function calculates the linearization to the dissipation term Dx
 * TODO there is probably a bug here.
 */
void Solver::calculateAndAddL(){
  Eigen::MatrixXd stencil_high_order(stencil_rows, stencil_cols);
  Eigen::MatrixXd stencil_low_order(stencil_rows, stencil_cols);
  stencil_high_order = Eigen::MatrixXd::Zero(stencil_rows, stencil_cols);
  stencil_low_order = Eigen::MatrixXd::Zero(stencil_rows, stencil_cols);

  Eigen::MatrixXd result(A.rows(), A.cols());
  result = Eigen::MatrixXd::Zero(A.rows(), A.cols());
  assert(result.cols()==A.cols()
      && result.rows() == A.rows()
      && "matrices for dissipation and system not same dimensions");


//  std::cout << "identity matrix: \n" << identity << std::endl;

//  std::cout << "result matrix: \n" << result << std::endl;
//  std::cout << "Solving high order stencil." << std::endl;
  for (int i = 0; i < data->q.size(); ++i)//loop through each node
  {
    //this sets up the vector of coefficients {lambda-1/2, lambda1/2} at each node point
    //noting that we use a constant extrapolation of these at the boundary
    std::vector<double> gamma(2);
    for (int j = 0; j < 2; ++j)
    {
      if (i-1 <0)
      {
        gamma[j] = problem->calc_lambda4_half(i);
      }
      else if (i+1 >= data->q.size())
      {
        gamma[j] = problem->calc_lambda4_half(i-1);
      }
      else
      {
        gamma[j] = problem->calc_lambda4_half(i-1+j);
      }
    }
    stencil_high_order.block<3,3>(i*3, i*3)     = -gamma[0]*(identity);
    stencil_high_order.block<3,3>(i*3, (i+1)*3) = (gamma[1] + 3*gamma[0])*(identity);
    stencil_high_order.block<3,3>(i*3, (i+2)*3) = -(3*gamma[0] + 3*gamma[1])*(identity);
    stencil_high_order.block<3,3>(i*3, (i+3)*3) = (3*gamma[1] + gamma[0])*(identity);
    stencil_high_order.block<3,3>(i*3, (i+4)*3) = -gamma[1]*(identity);
  }
  //solve for low order stencil
//  std::cout << "Calculating low order stencil." << std::endl;
  for (int i = 0; i < data->q.size(); ++i)//loop through each node
  {
    std::vector<double> gamma(2);
    for (int j = 0; j < 2; ++j)
    {
      std::vector<double> gamma(2);
      for (int j = 0; j < 2; ++j)
      {
        if (i-1 <0)
        {
          gamma[j] = problem->calc_lambda4_half(i);
        }
        else if (i+1 >= data->q.size())
        {
          gamma[j] = problem->calc_lambda4_half(i-1);
        }
        else
        {
          gamma[j] = problem->calc_lambda4_half(i-1+j);
        }
      }
    }
    stencil_low_order.block<3,3>(i*3, (i+1)*3) = gamma[0]*(identity);
    stencil_low_order.block<3,3>(i*3, (i+2)*3) = -(gamma[1] - gamma[0])*(identity);
    stencil_low_order.block<3,3>(i*3, (i+3)*3) = gamma[1]*(identity);
  }
  std::cout << "done" << std::endl;
  stencil_high_order += stencil_low_order;
  // std::cout << "Stencil before extraction: \n" << stencil_high_order << std::endl;
  //extract correct matrix
  for (int i=6; i < stencil_cols-6; ++i)
  {
    result.col(i - 6) = stencil_high_order.col(i);
  }
  // std::cout << "Matrix dissipation: \n" << 1.0/ data->dx * result << std::endl;
  result = problem->dt / data->dx * result;
  //std::cout << "matrix dissipation from e4 and e2 contribution\n" << result  <<std::endl;
  A = A - result;
}

/**
 * this function calculates dxE for each node and place that in the vector b
 */
void Solver::calcDe(){
  Eigen::VectorXd DxE = Eigen::VectorXd::Zero(data->getQVect().size(), 1);
  int j = 0;
  for(int i = 0; i<data->q.size(); ++i)
  {
    Eigen::Vector3d Ei = data->E(i+1) - data->E(i-1);
    double e1 = Ei(0);
    double e2 = Ei(1);
    double e3 = Ei(2);
    DxE(j) = e1;
    j++;
    DxE(j) = e2;
    j++;
    DxE(j) = e3;
    j++;
  }
  b = b - problem->dt*DxE/(2.0*data->dx);
}

/**
 *This function calculates the artificial dissipation for the RHS of the quasi Euler system and
 *adds it to member variable b.
 * */
void Solver::calcDx(){
  Eigen::VectorXd Dx = Eigen::VectorXd::Zero(data->getQVect().size());
  std::vector<Eigen::Vector3d> tmp_h(data->q.size());
  std::vector<Eigen::Vector3d> tmp_l(data->q.size());

  for (int i = 0; i < data->q.size(); ++i)//loop through all nodes
  {
    int sensor_node = i;
    double two_gamma_plus, two_gamma_minus, four_gamma_plus, four_gamma_minus;
    //special case for first node
    if (i-1 < 0)
    {
      two_gamma_plus  = problem->calc_lambda2_half(sensor_node);
      two_gamma_minus = two_gamma_plus;

      four_gamma_plus = problem->calc_lambda4_half(sensor_node);
      four_gamma_minus = four_gamma_plus;
    }
    else if(i+1 >= data->q.size())
    {
      two_gamma_minus = problem->calc_lambda2_half(sensor_node-1);
      two_gamma_plus = two_gamma_minus;

      four_gamma_minus = problem->calc_lambda4_half(sensor_node-1);
      four_gamma_plus = four_gamma_minus;
    }
    else
    {
      two_gamma_plus  = problem->calc_lambda2_half(sensor_node);
      two_gamma_minus = problem->calc_lambda2_half(sensor_node-1);

      four_gamma_plus = problem->calc_lambda4_half(sensor_node);
      four_gamma_minus = problem->calc_lambda4_half(sensor_node-1);
    }
    tmp_l[i] = two_gamma_plus* problem->lowOrderDifferencing(i)
            - two_gamma_minus * problem->lowOrderDifferencing(i-1);
    tmp_h[i] = four_gamma_plus* problem->highOrderDifferencing(i)
            - four_gamma_minus * problem->highOrderDifferencing(i-1);
    int ier = 3*i;
    Dx(ier  ) =  1/data->dx * (tmp_l[i](0) - tmp_h[i](0));
    Dx(ier+1) =  1/data->dx * (tmp_l[i](1) - tmp_h[i](1));
    Dx(ier+2) =  1/data->dx * (tmp_l[i](2) - tmp_h[i](2));
  }
  b = b + problem->dt * Dx;
}

void Solver::calcSx(){
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(data->getQVect().size());
  for(unsigned int i = 0; i < data->q.size(); ++i)
  {
    tmp(3*i) = 0;
    tmp(3*i+1) = data->Pressure(i) * data->Sprime(data->X(i));
    tmp(3*i+2) = 0;
  }
  b = b + problem->dt * tmp;
}

/**
 *This function calculates the dS/dQ from the quasi-euler equations and adds it to member A.
 * */
void Solver::calculateAndAddS(){
  int local_matrix_size = data->q[0].size();
  Eigen::MatrixXd local_S_matrix;
  Eigen::MatrixXd Sjac = Eigen::MatrixXd::Zero(data->getQVect().size(), data->getQVect().size());
  local_S_matrix = Eigen::MatrixXd::Zero(local_matrix_size, local_matrix_size);
  for(unsigned int i = 0; i < data->q.size(); ++i)
  {
    int ier = 3*i;
    local_S_matrix(1,0) = (data->parameter.gamma-1)/2 * std::pow(data->get_q2(i)/data->get_q1(i),2);
    local_S_matrix(1,1) = -(data->parameter.gamma-1) * data->get_q2(i)/data->get_q1(i);
    local_S_matrix(1,2) = (data->parameter.gamma-1.0);
    Sjac.block<3,3>(ier, ier) += local_S_matrix * data->Sprime(data->X(i));
  }
  A = A - problem->dt*Sjac;
}

/**
 *reinitializes A to identity and b to zero, use after each iteration.
 */
void Solver::reinit(){
  A = Eigen::MatrixXd::Identity(data->getQVect().size(), data->getQVect().size());
  b = Eigen::VectorXd::Zero(data->getQVect().size(), 1);
}

#endif
