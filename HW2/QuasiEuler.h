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

class QuasiEuler{

	public:

		const int local_matrix_size;

		//public functions
		double S(const double x);
			//TODO: write conversion functions and test them, could copy from hw1?
		Eigen::MatrixXd calculateLocalInviscidFluxJacobian(const StructuredGrid &data,
														   const int i) const;//do this for node i
		void calculateLocalDissipation();

		//helper functions
		void pressureSensor(StructuredGrid &data);
		void updateE(StructuredGrid &data);
		void updatePressure(StructuredGrid & data);
		void updateMach(StructuredGrid & data);
		void updateDensity(StructuredGrid & data);
		void updateTemp(StructuredGrid & data);
		void updatePhysicalQuantities(StructuredGrid &data);
		double calculateSoundSpeed(const StructuredGrid &data, const int index) const;
		void setInitialCondition(StructuredGrid &data);
		double nonlinearFunctionToSolveP1(double M,const StructuredGrid &data);
		double nonlinearFunctionToSolveP1Deriv(double M,const StructuredGrid &data);

		//public parameters (the problem parameters)
		double gamma, inlet_pressure, total_temperature, s_star;
		double initial_pressure_left, initial_pressure_right, time;
		double R;
		double dt = 0.1; //TODO: initialize this to something.


		QuasiEuler(StructuredGrid &mesh, double gamm, double inletP, double T, double sStar, double R_constant, double initial_time)
			: local_matrix_size(mesh.problem_dimension + 2), gamma(gamm),
			  inlet_pressure(inletP), total_temperature(T), s_star(sStar), time(initial_time),
			  R(R_constant)
		{
		  initial_pressure_left = -1;//TODO remove these...
		  initial_pressure_right = -1;
		};
		QuasiEuler operator=(const QuasiEuler& E){return *this;};

};

void QuasiEuler::pressureSensor(StructuredGrid &data){

	double kappa2 = 1./2.;
	double kappa4 = 1./50.;

	assert(data.sensor_contributions.cols() == 2 && "we only have two parameters for shock sensing, use only two columns");
	assert(data.mesh.size()==data.sensor_contributions.rows());

	for (int i = data.buffer_size; i < data.mesh.size()-data.buffer_size; ++i)//loop through all points on mesh
	{
		double topi = data.pressure(i+1) - 2*data.pressure(i) + data.pressure(i-1);
		double bottomi = data.pressure(i+1) + 2*data.pressure(i) + data.pressure(i-1);

		double topi_next = data.pressure(i+2) - 2*data.pressure(i+1) + data.pressure(i);
		double bottomi_next = data.pressure(i+2) + 2*data.pressure(i+1) + data.pressure(i);

		double topi_prev = data.pressure(i) - 2*data.pressure(i-1) + data.pressure(i-2);
		double bottomi_prev = data.pressure(i) + 2*data.pressure(i-1) + data.pressure(i-2);


		double GAMMA_i = std::abs(topi/bottomi);
		double GAMMA_i_next = std::abs(topi_next/bottomi_next);
		double GAMMA_i_prev = std::abs(topi_prev/bottomi_prev);

		data.sensor_contributions(i,0) = kappa2*std::max(GAMMA_i, std::max(GAMMA_i_next, GAMMA_i_prev));
		data.sensor_contributions(i,1) = kappa4*std::max(0., kappa4 - data.sensor_contributions(i,0));

	}
}

void QuasiEuler::setInitialCondition(StructuredGrid &data){
  double density = inlet_pressure/ (R* total_temperature);

  double mach = 0.5;

  while(std::fabs(nonlinearFunctionToSolveP1(mach, data))>1e-13)
  {
    mach = mach -nonlinearFunctionToSolveP1(mach, data)/nonlinearFunctionToSolveP1Deriv(mach, data);
  }

  double velocity = std::sqrt(gamma*inlet_pressure/density) * mach;
  double momentum = velocity * density;
  double energy = density * total_temperature * R / (gamma-1);
  double energy_momentum = density * energy;

  std::cout << "density = " << density << std::endl
            << "momentum = " << momentum << std::endl
            << "energy momentum = " << energy_momentum << std::endl
            << "mach = " << mach << std::endl;
  assert(data.Q.size() - data.buffer_size == data.stop_iteration_index && "you are attempting to index outside of the range of Q");
  for (int i = data.buffer_size; i <data.stop_iteration_index; ++i)
  {
    data.Q[i](0) = density;
    data.Q[i](1) = momentum;
    data.Q[i](2) = energy_momentum;
  }

}

Eigen::MatrixXd QuasiEuler::calculateLocalInviscidFluxJacobian(const StructuredGrid & data,
		const int i) const{

	Eigen::MatrixXd local_flux(local_matrix_size, local_matrix_size);
//TODO: make this so that I can construct the matrix independent of dimension, right now
	//it is hard coded.

	local_flux(0,0) = 0;
	local_flux(0,1) = 1;
	local_flux(0,2) = 0;
	local_flux(1,0) = (gamma-3)/2*std::pow(data.Q[i](1)/data.Q[i](0), 2);
	local_flux(1,1) = (3-gamma)*data.Q[i](1)/data.Q[i](0);
	local_flux(1,2) = gamma-1;
	local_flux(2,0) = (gamma-1)*std::pow(data.Q[i](1)/data.Q[i](0), 3)
							- gamma*(data.Q[i](2)/data.Q[i](0))*(data.Q[i](1)/data.Q[i](0));
	local_flux(2,1) = gamma*(data.Q[i](2)/data.Q[i](0))
							- (3*(gamma-1))/(2)*std::pow(data.Q[i](1)/data.Q[i](0),2);
	local_flux(2,2) = gamma*(data.Q[i](1)/data.Q[i](0));

//	std::cout<<"local flux at " << i << " is\n " << local_flux << std::endl;

	return local_flux;
}

void QuasiEuler::updatePhysicalQuantities(StructuredGrid &data){
  data.interpolateBoundary();
  updateDensity(data);
  updatePressure(data);
  updateMach(data);
  updateE(data);
}

void QuasiEuler::updatePressure(StructuredGrid & data){
	assert(data.Q_has_been_updated && "you need to update Q before updating pressure");
	for (int i = data.buffer_size; i < data.stop_iteration_index; ++i)
	{
		double q1 = data.Q[i](0);
		double q2 = data.Q[i](1);
		double q3 = data.Q[i](2);
		data.pressure(i) = (gamma-1)*(q3/q1 - 1/(2*q1)*std::pow(q2,2));
	}
	data.Pressure_has_been_updated = 1;
}

void QuasiEuler::updateMach(StructuredGrid & data){
	assert(data.Q_has_been_updated
	       && "you need to update Q before updating other variables");
	assert(data.Pressure_has_been_updated
	       && "you need to update pressure before updating the mach");
	assert(data.Density_has_been_updated
	       && "you need to update density before updating the mach");
	for (int i = data.buffer_size; i < data.stop_iteration_index; ++i)
	{
	  double q1 = data.Q[i](0);
	  std::cout << "q1 is" << q1 << std::endl;
	  double q2 = data.Q[i](1);
	  std::cout << "q2 is" << q2 << std::endl;
	  double q3 = data.Q[i](2);
	  std::cout << "q3 is" << q3 << std::endl;
		data.mach(i) = (q2/q1)/(QuasiEuler::calculateSoundSpeed(data, i));
	}
}

void QuasiEuler::updateDensity(StructuredGrid & data){
	assert(data.Q_has_been_updated && "you need to update Q before updating density");
	for (int i = data.buffer_size; i < data.stop_iteration_index; ++i)
	{
	  double q1 = data.Q[i](0);
	  double q2 = data.Q[i](1);
	  double q3 = data.Q[i](2);
		data.density(i) = q1;
	}
	data.Density_has_been_updated = 1;
}

void QuasiEuler::updateTemp(StructuredGrid & data){
  assert(0 && "not implemented yet.");
  assert(data.Q_has_been_updated
         && "you need to update Q before updating other variables");
  assert(data.Pressure_has_been_updated
         && "you need to update pressure before updating the temperature");
  assert(data.Density_has_been_updated
         && "you need to update density before updating the temperature");

  for (int i = data.buffer_size; i < data.stop_iteration_index; ++i)
  {
    double q1 = data.Q[i](0);
    double q2 = data.Q[i](1);
    double q3 = data.Q[i](2);
//    data.temperature(i) = data.pressure(i)/(data);
  }
}

void QuasiEuler::updateE(StructuredGrid &data){
  assert(data.Pressure_has_been_updated
         && "need pressure to calculate E");
  for (int i = 0; i < data.E.size(); ++i)
    {
      double q1 = data.Q[i](0);
      double q2 = data.Q[i](1);
      double q3 = data.Q[i](2);
      data.E[i](0) = q2;
      data.E[i](1) = q2 + data.pressure(i);
      data.E[i](1) = q3 + data.density(i)*data.pressure(i) - data.pressure(i);
    }
  data.E_has_been_updated = 1;
}

double QuasiEuler::calculateSoundSpeed(const StructuredGrid &data, const int index) const{
  std::cout << "density at index " << index << " is " << data.density(index) << std::endl;
  std::cout << "pressure at index " << index << " is " << data.pressure(index) << std::endl;
  return std::sqrt(gamma*data.pressure(index)/data.density(index));
//  return gamma*data.R* data.
}

double QuasiEuler::S(const double x){
  if(x <= 5)
  {
    return 1.0 + 1.5*std::pow(1-x/5,2);
  }
  else if(x >5)
  {
    return 1.0 + 0.5*std::pow(1-x/5,2);
  }
  else return 0.0;
}

double QuasiEuler::nonlinearFunctionToSolveP1(double M,const StructuredGrid &data){

  double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
  double exponent = (gamma + 1)/(2*(gamma-1));
  assert(fabs(exponent - 3)<1e-13);
  return pow(inside, exponent) - (S(data.mesh(data.buffer_size))*M)/(s_star);
}

double QuasiEuler::nonlinearFunctionToSolveP1Deriv(double M,const StructuredGrid &data){
  double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
  double exponent = (gamma + 1)/(2*(gamma-1));
  return M*pow(inside, exponent-1) - (S(data.mesh(data.buffer_size)))/(s_star);
}

















#endif
