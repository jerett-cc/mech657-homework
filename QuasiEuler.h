#ifndef QUASIEULER_H
#define QUASIEULER_H

/*
 * This class includes all the functions to calculate local fluxes,
 * RHS contributions and all functions needed to calculate
 * physical quantities of interest from the solved quantities,
 * relations include sound speed, temperature, and pressure.
 *
 * this class also senses shocks with a pressure sensor found in
 * "Fundamentals Algorithms in Computational Fluid Dynamics"
 * by Pulliam and Zingg 2014 pp. 97
 */

#include <cmath>
#include <ostream>
#include <cassert>

#include "../Eigen/Dense"
#include "Geometry.h"


/**
 * this class should store the sensor contributions and also be able to calculate the
 * dissipation at an index idx of the data.
 *
 * dissipation will be used for the rhs,
 *
 * the sensor contributions will be used for ^^ and also for the contribution to L?
 */
class QuasiEuler{

  public:
    double dt;
    ProblemData *data;
    Eigen::MatrixXd sensor_contributions;

    QuasiEuler(ProblemData *Data, double dt)
      : data(Data), dt(dt)
    {
      sensor_contributions = Eigen::MatrixXd::Zero(data->q.size(), 2);
    }

    void calculateSensorContributions();

    Eigen::Matrix3d calculateLocalFluxJacobian(const int idx);

    Eigen::Vector3d lowOrderDifferencing(const int idx);//todo
    Eigen::Vector3d highOrderDifferencing(const int idx);//todo


};

/**
 * calculates he coefficient of the second derivative term used for the rhs at each
 * interior node point only,
 * stores kappa2 in column 0 and kappa4 in column 1 of the matrix.
 */
void QuasiEuler::calculateSensorContributions(){

  double kappa2 = 0.5;
  double kappa4 = 0.02;


  assert(data->q.size()==sensor_contributions.rows());

  for (int i = 0; i < data->q.size(); ++i)//loop through all points on mesh
  {
//    std::cout << "iteration "<< i << "--------------------------------------------" << std::endl;
    double topi = data->Pressure(i+1) - 2*data->Pressure(i) + data->Pressure(i-1);
//    std::cout << "Top i " << std::endl
//        << topi << std::endl;
    double bottomi = data->Pressure(i+1) + 2*data->Pressure(i) + data->Pressure(i-1);
//    std::cout << "bot i " << std::endl
//        << bottomi << std::endl;

    double topi_next = data->Pressure(i+2) - 2*data->Pressure(i+1) + data->Pressure(i);
//    std::cout << "Top i+1 " << std::endl
//        << topi_next << std::endl;
    double bottomi_next = data->Pressure(i+2) + 2*data->Pressure(i+1) + data->Pressure(i);
//    std::cout << "bot i+1 " << std::endl
//        << bottomi_next << std::endl;

    double topi_prev = data->Pressure(i) - 2*data->Pressure(i-1) + data->Pressure(i-2);
//    std::cout << "Top i-1 " << std::endl
//        << topi_prev << std::endl;
    double bottomi_prev = data->Pressure(i) + 2*data->Pressure(i-1) + data->Pressure(i-2);
//    std::cout << "bot i-1 " << std::endl
//        << bottomi_prev << std::endl;


    double GAMMA_i = std::abs(topi/bottomi);
    double GAMMA_i_next = std::abs(topi_next/bottomi_next);
    double GAMMA_i_prev = std::abs(topi_prev/bottomi_prev);
    double epsilon2 =  kappa2*std::max(GAMMA_i, std::max(GAMMA_i_next, GAMMA_i_prev));
    sensor_contributions(i,0) = epsilon2;
    sensor_contributions(i,1) = std::max(0., kappa4 - epsilon2);

  }

}

/**
 * this function should return a 3x3 local flux jacobian to an interior point at index
 * idx. this should only be used for the construction of the system matrix
 *
 */
Eigen::Matrix3d QuasiEuler::calculateLocalFluxJacobian(const int idx){
  Eigen::Matrix3d local_flux;

//  std::cout << local_flux.cols() << "  " << local_flux.rows() << std::endl;

  double gamma = data->parameter.gamma;
  double q1 = data->q[idx](0);
  double q2 = data->q[idx](1);
  double q3 = data->q[idx](2);
  local_flux(0,0) = 0;
  local_flux(0,1) = 1;
  local_flux(0,2) = 0;

  local_flux(1,0) = (gamma-3)/2*std::pow(data->q[idx](1)/data->q[idx](0), 2);
  local_flux(1,1) = (3-gamma)*data->q[idx](1)/data->q[idx](0);
  local_flux(1,2) = gamma-1;

  local_flux(2,0) = (gamma-1)*std::pow(q2/q1, 3) - gamma*(q3/q1)*(q2/q1);
  local_flux(2,1) = gamma*(data->q[idx](2)/data->q[idx](0))
                    - (3*(gamma-1))/(2)*std::pow(data->q[idx](1)/data->q[idx](0),2);
  local_flux(2,2) = gamma*(q2/q1);

  return local_flux;
}

Eigen::Vector3d QuasiEuler::lowOrderDifferencing(const int idx){
  double u_1 = data->Velocity(idx+1);
  double u_0 = data->Velocity(idx);
  double a_1 = data->soundSpeed(idx+1);
  double a_0 = data->soundSpeed(idx);
  double sigma_1 = sensor_contributions(idx + 1)*(std::abs(u_1) + a_1);
  double sigma_0 = sensor_contributions(idx)*(std::abs(u_0) + a_0);

  return (sigma_1 + sigma_0)/2 * (data->Q(idx+1) - data->Q(idx));
}

Eigen::Vector3d QuasiEuler::highOrderDifferencing(const int idx){
  double u_1 = data->Velocity(idx+1);
  double u_0 = data->Velocity(idx);
  double a_1 = data->soundSpeed(idx+1);
  double a_0 = data->soundSpeed(idx);
  double sigma_1 = sensor_contributions(idx + 1)*(std::abs(u_1) + a_1);
  double sigma_0 = sensor_contributions(idx)*(std::abs(u_0) + a_0);

  return (sigma_1 + sigma_0)/2 * (data->Q(idx+2) - 3* data->Q(idx+1) + 3*data->Q(idx) - data->Q(idx-1));
}



























#endif
