#ifndef SpatialLinearSystemCons_H
#define SpatialLinearSystemCons_H


#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"
#include "../Eigen/Sparse"
#include "QuasiEuler.h"
#include "Geometry.h"

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
      assert(A.cols() == data->q.size()* 3);
      //todo remove this assertion, it will not be true in more than one dim
    }



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
    void reinit();

    void setup_b();
    void calcDe();
    void calcDx();

};

void Solver::setupSystem(){
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
}

/**
 *this function calculates An and L and adds them to A. */
void Solver::setup_A(){
  calculateAndAddSpatialMatrix();
  calculateAndAddL();
};
/**
 * this function calculates all the local flux jacobians and does an second order approx
 * to the total flux jacobian in space. It stitches each together in the system
 * matrix A.
 */
void Solver::calculateAndAddSpatialMatrix(){
  int local_matrix_size = data->q[0].size();
  Eigen::MatrixXd Ai(local_matrix_size, local_matrix_size);
  Eigen::MatrixXd Ai_n(local_matrix_size, local_matrix_size);

  for (int i = 0; i < data->q.size()-1; ++i)//loop through all nodes but last one
  {
//    std::cout<< i << std::endl;
    Ai = problem->calculateLocalFluxJacobian(i);//calculate flux at node i
    Ai_n = problem->calculateLocalFluxJacobian(i+1);//calculate flux at node i+1

    int sub_index = local_matrix_size * (i);

    assert(sub_index+local_matrix_size + 3<=A.cols() && "indexing outside of bounds in the matrix");
//        std::cout << "(" << sub_index << " , " <<sub_index+local_matrix_size << ")" << std::endl;
    A.block<3,3>(sub_index, sub_index+local_matrix_size) += 0.5*problem->dt/data->dx *Ai_n;
    //TODO divide by dx here??
    A.block<3,3>(sub_index + local_matrix_size, sub_index) += -0.5*problem->dt/data->dx*Ai;
  }
//    std::cout << "System matrix excluding first and last " << std::endl << A << std::endl;
}

/**
 * thid function calculates the linearization to the dissipation term Dx
 * todo there is probably a bug here.
 */
void Solver::calculateAndAddL(){
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
  std::cout << "Done." << std::endl;
//  std::cout << stencil_high_order << std::endl;

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

  std::cout << result << "\n stencil from e4 and e2 contribution ^^" <<std::endl;

  A += result;

//  std::cout << "A is now " << std::endl << A << std::endl;
}
/**
 * this function should calculate dxE for each node and place that in the vector b
 */
void Solver::calcDe(){
//  std::cout << "calc de called here ++++++++++++++++" << std::endl;
  int j = 0;
  for(int i = 0; i<data->q.size(); ++i)
  {
    Eigen::Vector3d Ei= data->E(i+1) - data->E(i-1);
//    std::cout << "index is " << i <<" next E is \n " << data->E(i+1) << std::endl;
//    std::cout << "index is " << i <<" previous E is \n " << data->E(i-1) << std::endl;
    std::cout << "index is " << i <<" Ei is \n " << Ei << std::endl;

    double e1 = Ei(0);
    double e2 = Ei(1);
    double e3 = Ei(2);
//    b(j) = e1/(2);
//    b(j+1) = e2/(2);
//    b(j+2) = e3/(2);
//    std::cout << "j is " << j << std::endl;
    b(j) = e1/(2*data->dx);
    b(j+1) = e2/(2*data->dx);
    b(j+2) = e3/(2*data->dx);
    j+=3;
  }
  std::cout << "dex is " << std::endl << -1*problem->dt * b << std::endl;
  b = -1*problem->dt * b;
}

/**
 * calculates the dissipation contribution to the RHS, and adds it to the member function b
 */
void Solver::calcDx(){
  std::vector<Eigen::Vector3d> tmp_h(data->q.size());//we will overwrite this could introduce bugs
  std::vector<Eigen::Vector3d> tmp_l(data->q.size());//these are indermediate

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
  std::cout << "B before \n" << b << std::endl;
  for(int i = 0; i < data->q.size(); ++i )//loop through all nodes
  {
    for(int j = 0; j<3; ++j)//loop through all components
    {
      b(i+j) -= 1/data->dx * (tmp_l[i](j) - tmp_l[i-1](j)) - 1/data->dx * (tmp_l[i](j) - tmp_l[i-1](j));
    }
  }


  std::cout << "E-Dx is \n" << b << std::endl;

}
/**
 *reinitializes A to identity and b to zero, use after each iteration. */
void Solver::reinit(){
  A =  Eigen::MatrixXd::Identity(data->getQVect().size(), data->getQVect().size());
  b = Eigen::VectorXd::Zero(data->getQVect().size(), 1);

}















#endif
