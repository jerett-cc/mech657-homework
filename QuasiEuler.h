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
    double dt,cfl;
    ProblemData *data;
    Eigen::MatrixXd sensor_contributions;
    QuasiEuler(ProblemData *Data, double cfl, double max_vel)
      : data(Data), cfl(cfl)
    {
      dt = cfl * data->dx / (max_vel + 300.);//FIXME: how to change this for adaptive time stepping?
      sensor_contributions = Eigen::MatrixXd::Zero(data->q.size(), 2);
    }
    void calculateSensorContributions();
    double sigma(const int idx);
    double calc_lambda2_half(const int);
    double calc_lambda4_half(const int);
    Eigen::Matrix3d calculateLocalFluxJacobian(const int idx);
    Eigen::Vector3d lowOrderDifferencing(const int idx);//todo
    Eigen::Vector3d highOrderDifferencing(const int idx);//todo
    Eigen::Vector3d boundaryDifferencing(const int idx);
};

/**
 * calculates he coefficient of the second derivative term used for the rhs at each
 * interior node point only,
 * stores kappa2 in column 0 and kappa4 in column 1 of the matrix.
 */
void QuasiEuler::calculateSensorContributions(){
  double kappa2 = 0.5;
  //kappa2 = 0.0;//FIXME
  double kappa4 = 0.02;
  //kappa4 = 0.;
  //kappa2 = 0.;
  assert(data->q.size()==sensor_contributions.rows());
  for (int i = 0; i < data->q.size(); ++i)
  {
    double topi = data->Pressure(i+1) - 2*data->Pressure(i) + data->Pressure(i-1);
    double bottomi = data->Pressure(i+1) + 2*data->Pressure(i) + data->Pressure(i-1);

    double topi_next = data->Pressure(i+2) - 2*data->Pressure(i+1) + data->Pressure(i);
    double bottomi_next = data->Pressure(i+2) + 2*data->Pressure(i+1) + data->Pressure(i);

    double topi_prev = data->Pressure(i) - 2*data->Pressure(i-1) + data->Pressure(i-2);
    double bottomi_prev = data->Pressure(i) + 2*data->Pressure(i-1) + data->Pressure(i-2);

    double GAMMA_i = std::abs(topi/bottomi);
    double GAMMA_i_next = std::abs(topi_next/bottomi_next);
    double GAMMA_i_prev = std::abs(topi_prev/bottomi_prev);

    if(i-1<0)
    {
      GAMMA_i_prev = GAMMA_i;
    }
    else if (i+1==data->q.size())
    {
      GAMMA_i_next = GAMMA_i;
    }
    double epsilon2 =  kappa2*std::max(GAMMA_i, std::max(GAMMA_i_next, GAMMA_i_prev));
    sensor_contributions(i,0) = epsilon2;
    sensor_contributions(i,1) = std::max(0., kappa4 - epsilon2);
  }
}

/**
 * calculates the largest eigenvalue of the flux jacobian matrix times the low order sensor
 * used to compute RHS dissipation and the matrix L
 */
double
QuasiEuler::calc_lambda2_half(const int a_i)
{
  double lambda_2;
  if (a_i < data->q.size()-1)
  {
    lambda_2 = 0.5 *
        (sensor_contributions(a_i,0) * sigma(a_i) + sensor_contributions(a_i+1, 0) * sigma(a_i+1));
  }
  else
  {
    lambda_2 = sensor_contributions(a_i,0) * sigma(a_i);
  }
  return lambda_2;
}

/**
 * calculates the largest eigenvalue of the flux jacobian matrix times the high order sensor
 * used to compute RHS dissipation and the matrix L
 */
double
QuasiEuler::calc_lambda4_half(const int a_i)
{
  double lambda_4;
  if (a_i < data->q.size())
  {
    lambda_4 = 0.5 *
        (sensor_contributions(a_i,1) * sigma(a_i) + sensor_contributions(a_i+1, 1) * sigma(a_i+1));
  }
  else
  {
    lambda_4 = sensor_contributions(a_i,1) * sigma(a_i);
  }
  return lambda_4;
}

/**
 * this function should return a 3x3 local flux jacobian to an interior point at index
 * idx. this should only be used for the construction of the system matrix
 *
 */
Eigen::Matrix3d
QuasiEuler::calculateLocalFluxJacobian(const int idx)
{
  Eigen::Matrix3d local_flux;
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

double QuasiEuler::sigma(const int idx){
  return std::abs(data->Velocity(idx)) + data->soundSpeed(idx);
}

Eigen::Vector3d QuasiEuler::lowOrderDifferencing(const int idx){
  return (data->Q(idx+1) - data->Q(idx));
}

Eigen::Vector3d QuasiEuler::highOrderDifferencing(const int idx){
  return (data->Q(idx+2) - 3* data->Q(idx+1) + 3*data->Q(idx) - data->Q(idx-1));
}

Eigen::Vector3d QuasiEuler::boundaryDifferencing(const int idx){
  return (-data->Q(idx+1) + 2*data->Q(idx) - data->Q(idx-1));
}

#endif
