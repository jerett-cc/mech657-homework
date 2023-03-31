#ifndef SpatialLinearSystemCons_H
#define SpatialLinearSystemCons_H


#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"
#include "../Eigen/Sparse"
#include "QuasiEuler.h"
#include "Geometry.h"
#include <iomanip>
#include <cmath>

// needs to calculate the matrix, calculate the RHS, and solve the system for deltaQ

/**
 * this class should take the data structure, the problem (which has the pressure sensor and
 * dissipation calculation) and produce a matrix A, the rhs b, and can solve it with the iterative method
 * and the diagonal form
 */
class Solver{

  public:
    std::vector<Eigen::Vector3d> delta_Q;//solution vector
    double step_error = 1.0;
    void setupSystem();//does this need to include a parameter Data Class?
    void solveSystem();//updates delta_Q with the result of the system solution.



    Solver(ProblemData *data, QuasiEuler * problem)
      : data(data), problem(problem)
      , stencil_rows(data->q.size()*3)
      , stencil_cols(3*(data->q.size()+4))
    {
      delta_Q = data->q;//we will overwrite the values in this. when we compute delta_Q
      //this could introduce bugs. be aware?
//      std::cout << "Look here "  <<  data->getQVect().size() << std::endl;
      A = Eigen::MatrixXd::Identity(data->getQVect().size(), data->getQVect().size());
      b = Eigen::VectorXd::Zero(data->getQVect().size(), 1);
      assert(A.cols() == data->q.size() * 3);
      //todo remove this assertion, it will not be true in more than one dim
    }
    double error();


//  private:
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    int stencil_cols;
    int stencil_rows;
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
    void XXXcalcDx();
    void calcSx();
};

double Solver::error(){
  double tmp = 0;
  return -300.0;
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
  //calcDx();
  XXXcalcDx();
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
  std::cout << "dxA" << std::endl << TOSHOW << std::endl;
  A =  A + problem->dt*TOSHOW;
}

/**
 * this function calculates the linearization to the dissipation term Dx
 * TODO there is probably a bug here.
 */
void Solver::calculateAndAddL(){
  Eigen::MatrixXd VIEWSTENCIL =
    Eigen::MatrixXd::Zero(data->getQVect().size(), data->getQVect().size());
  Eigen::MatrixXd stencil_high_order(stencil_rows, stencil_cols);
  Eigen::MatrixXd stencil_low_order(stencil_rows, stencil_cols);

  Eigen::MatrixXd result(A.rows(), A.cols());
  result = Eigen::MatrixXd::Zero(A.rows(), A.cols());
  assert(result.cols()==A.cols()
      && result.rows() == A.rows()
      && "matrices for dissipation and system not same dimensions");

//  std::cout<< result << std::endl;
  //solve for the high order stencil
  std::cout << "Solving high order stencil." << std::endl;
  for (int i = 0; i < data->q.size(); ++i)//loop through each node
  {
    std::vector<double> gamma(2);
    for (int j = 0; j < 2; ++j)
    {
      if (i-1+j <0)
      {
        std::cout << "i-1+j <0" << std::endl;
        double u_0j = data->Velocity(i-1+j);
        double u_1j = data->Velocity(i+j);
        double a_0j = data->soundSpeed(i-1+j);
        double a_1j = data->soundSpeed(i+j);

        double sigma1j = problem->sensor_contributions(i+j,1) * (std::abs(u_1j) + a_1j);
        double sigma0j = sigma1j;
        gamma[j] = 0.5 * (sigma1j + sigma0j);

        std::cout << "GAMMA(-1) is \n" << gamma[0] << std::endl;
        std::cout << "done." << std::endl;
      }
      else if (i+j > problem->sensor_contributions.rows()-1)
      {
        std::cout << "i + j >rows-1" << std::endl;
        double u_0j = data->Velocity(i-1+j);
        double u_1j = data->Velocity(i+j);
        double a_0j = data->soundSpeed(i-1+j);
        double a_1j = data->soundSpeed(i+j);
        double sigma0j = problem->sensor_contributions(i-1+j,1) * (std::abs(u_0j) + a_0j);
        double sigma1j = sigma0j;
        gamma[j] = 0.5 * (sigma1j + sigma0j);
        std::cout << "done." << std::endl;
      } else {
        std::cout << "internal case" << std::endl;
        double u_0j = data->Velocity(i-1+j);
        double u_1j = data->Velocity(i+j);
        double a_0j = data->soundSpeed(i-1+j);
        double a_1j = data->soundSpeed(i+j);
        double sigma1j = problem->sensor_contributions(i+j,1) * (std::abs(u_1j) + a_1j);
        double sigma0j = problem->sensor_contributions(i-1+j,1) * (std::abs(u_0j) + a_0j);
        gamma[j] = 0.5 * (sigma1j + sigma0j);
        std::cout << "done." << std::endl;
      }
    }
//    double gammaj = problem->sensor_contributions(i,1)*;
    stencil_high_order.block<3,3>(i*3, i*3) = gamma[0]*(1*identity);
    stencil_high_order.block<3,3>(i*3, (i+1)*3) = -(gamma[1] + 3*gamma[0])*(-4*identity);
    stencil_high_order.block<3,3>(i*3, (i+2)*3) = (3*gamma[0] + 3*gamma[1])*(6*identity);
    stencil_high_order.block<3,3>(i*3, (i+3)*3) = -(3*gamma[1] + gamma[0])*(-4*identity);
    stencil_high_order.block<3,3>(i*3, (i+4)*3) = gamma[1]*(1*identity);
  }
  std::cout << "Done. Stencil high order is: " << std::endl;
  std::cout << stencil_high_order << std::endl;
  //solve for low order stencil
  std::cout << "Calculating low order stencil." << std::endl;
  for (int i = 0; i < data->q.size(); ++i)//loop through each node
  {
    std::vector<double> gamma(2);
    for (int j = 0; j < 2; ++j)
    {
      if (i-1+j <0)
      {
        double u_1j = data->Velocity(i+j);
        double a_1j = data->soundSpeed(i+j);

        double sigma1j = problem->sensor_contributions(i+j,0) * (std::abs(u_1j) + a_1j);
        double sigma0j = sigma1j;
        gamma[j] = 0.5 * (sigma1j + sigma0j);
      }
      else if (i+j > problem->sensor_contributions.rows()-1)
      {
        double u_0j = data->Velocity(i-1+j);
        double a_0j = data->soundSpeed(i-1+j);

        double sigma0j = problem->sensor_contributions(i-1+j,0) * (std::abs(u_0j) + a_0j);
        double sigma1j = sigma0j;

        gamma[j] = 0.5 * (sigma1j + sigma0j);
      } else {
        double u_0j = data->Velocity(i-1+j);
        double u_1j = data->Velocity(i+j);
        double a_0j = data->soundSpeed(i-1+j);
        double a_1j = data->soundSpeed(i+j);

        double sigma1j = problem->sensor_contributions(i+j,0) * (std::abs(u_1j) + a_1j);
        double sigma0j = problem->sensor_contributions(i-1+j,0) * (std::abs(u_0j) + a_0j);

        gamma[j] = 0.5 * (sigma1j + sigma0j);
      }
    }
    stencil_low_order.block<3,3>(i*3, (i+1)*3) = gamma[0]*(-1*identity);
    stencil_low_order.block<3,3>(i*3, (i+2)*3) = -(gamma[0] + gamma[1])*(2*identity);
    stencil_low_order.block<3,3>(i*3, (i+3)*3) = gamma[1]*(-1*identity);
  }
  std::cout << "done" << std::endl;
  stencil_high_order += stencil_low_order;
//   std::cout << stenci  l_high_order << "\n stencil ^^" <<std::endl;

  //extract correct matrix
  for (int i=6; i < stencil_cols-6; ++i)
  {
    result.col(i - 6) = stencil_high_order.col(i);
  }
  result = problem->dt * result;
  std::cout << "matrix dissipation from e4 and e2 contribution\n" << result  <<std::endl;
  A = A + result;
//  std::cout << "A is now " << std::endl << A << std::endl;
}

/**
 * this function should calculate dxE for each node and place that in the vector b
 */
void Solver::calcDe(){
  Eigen::VectorXd DxE = Eigen::VectorXd::Zero(data->getQVect().size(), 1);
  int j = 0;
  std::cout << "Printing E ____" <<std::endl;
  for(int i = 0; i<data->q.size(); ++i)
  {
    Eigen::Vector3d Ei = data->E(i+1) - data->E(i-1);
    std::cout << std::setprecision(16) <<  data->E(i) << std::endl;
    double e1 = Ei(0);
    double e2 = Ei(1);
    double e3 = Ei(2);
    DxE(j) = e1;
    //DxE(j) = 1;
    j++;
    DxE(j) = e2;
    //DxE(j) = 1;
    j++;
    DxE(j) = e3;
    //DxE(j) = 1;
    j++;
  }
  std::cout << "DxE with /2dx\n" << std::setprecision(16) << DxE/(2.0*data->dx) << std::endl;
  b = b - problem->dt*DxE/(2.0*data->dx);
  //std::cout << "h*DxE\n" << std::setprecision(16) << problem->dt*DxE/(2.0*data->dx) << std::endl;
}

/**
 * calculates the dissipation contribution to the RHS, and adds it to the member function b
 */
void Solver::calcDx(){
  std::vector<Eigen::Vector3d> tmp_h(data->q.size());//we will overwrite this could introduce bugs
  std::vector<Eigen::Vector3d> tmp_l(data->q.size());//these are indermediate
  Eigen::VectorXd c = Eigen::VectorXd::Zero(data->getQVect().size());
  for(int i = 0; i<data->q.size(); ++i)//loop through all nodes, this is before the last backward differencing.
  {
//    std::cout << "index dummy " << i << std::endl;
    tmp_h[i] = problem->highOrderDifferencing(i);
    tmp_l[i] = problem->lowOrderDifferencing(i);
//    std::cout << "-------------" << std::endl
//              << tmp_h[i] << std::endl
//              << "-------------" << std::endl
//              << tmp_l[i] << std::endl
//              << "-------------" << std::endl;
  }
  //std::cout << "B before \n" << b << std::endl;
  for(int i = 0; i < data->q.size(); ++i )//loop through all nodes
  {
    int ier = 3*i;
    Eigen::Vector3d back_tmph;
    Eigen::Vector3d back_tmpl;
    if (i-1 <0)
    {
      back_tmph = 0*tmp_h[i];
      back_tmpl = 0*tmp_l[i];
    }
    else
    {
      back_tmph = tmp_h[i] - tmp_h[i-1];
      back_tmpl = tmp_l[i] - tmp_l[i-1];
    }
    c(ier    ) = -problem->dt * (-back_tmph(0) + back_tmpl(0));
    c(ier + 1) = -problem->dt * (-back_tmph(1) + back_tmpl(1));
    c(ier + 2) = -problem->dt * (-back_tmph(2) + back_tmpl(2));

  }
  std::cout << "Dx is \n" << c << std::endl;
  std::cout << "____________________________\n" << "Q0:\n" << data->Q(0) << std::endl;
  std::cout<< "Q-1\n" << data->Q(-1) << "\n and Q-2 \n" << data->Q(-2) << std::endl;
  std::cout<< "Q1\n" << data->Q(1) << "\n and Q2 \n" << data->Q(2) << std::endl;
}

//FIXME: add a description here
void Solver::XXXcalcDx(){
  Eigen::VectorXd Dx = Eigen::VectorXd::Zero(data->getQVect().size());
  std::vector<Eigen::Vector3d> tmp_h(data->q.size());
  std::vector<Eigen::Vector3d> tmp_l(data->q.size());

  for (int i = 0; i < data->q.size(); ++i)//loop through all nodes
  {
    //start at 1, end at size
    std::cout << "i: " << i << std::endl;
    std::cout << "lambda2(i): \n" << std::setprecision(16) << problem->calc_lambda2_half(i) << std::endl;
    std::cout << "calc low order: \n" << std::setprecision(16)<< problem->lowOrderDifferencing(i) << std::endl;
    std::cout << "calc high order: \n" << std::setprecision(16)<< problem->highOrderDifferencing(i) << std::endl;

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
    Dx(ier  ) =  1/data->dx * (tmp_l[i](0) - tmp_h[i](0));//FIXME: compare this with andy's code
    Dx(ier+1) =  1/data->dx * (tmp_l[i](1) - tmp_h[i](1));
    Dx(ier+2) =  1/data->dx * (tmp_l[i](2) - tmp_h[i](2));

//    Dx(ier  ) =  1/(10.0/60.0) * (tmp_l[i](0) - tmp_h[i](0));//FIXME: delete this
//    Dx(ier+1) =  1/(10.0/60.0) * (tmp_l[i](1) - tmp_h[i](1));
//    Dx(ier+2) =  1/(10.0/60.0) * (tmp_l[i](2) - tmp_h[i](2));
  }

  std::cout << "____________Dissipation (divided by dx) is: _______________" << std::endl;
  std::cout << Dx << std::endl;
  b = b + problem->dt * Dx;
}

void Solver::calcSx(){
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(data->getQVect().size());

  for(unsigned int i = 0; i < data->q.size(); ++i)
  {
    tmp(3*i) = 0;
    tmp(3*i+1) = data->Pressure(i) * data->Sprime(data->X(i));
    tmp(3*i+2) = 0;
    //tmp(3*i) = 1;
    //tmp(3*i+1) = 1;
    //tmp(3*i+2) = 1;
  }
  //std::cout << "pressure at index 0 is " << data->Pressure(0) << std::endl;
  std::cout << "S contributions to RHS\n" << std::setprecision(16) << tmp << std::endl;
  b = b + problem->dt * tmp;
}

/**
 *
 *Working!
 *
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

  std::cout << "Sjac is \n" << Sjac << std::endl;
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
