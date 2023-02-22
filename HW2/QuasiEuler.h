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
		void S(const StructuredGrid &mesh);//TODO: return S(x) at mesh points
			//TODO: write conversion functions and test them, could copy from hw1?
		Eigen::MatrixXd calculateLocalInviscidFluxJacobian(const StructuredGrid &data,
					const int i) const;//do this for node i
		void calculateLocalDissipation();

		//helper functions
		void pressureSensor(StructuredGrid &data);
		void updatePressure(StructuredGrid & data);
		void updateMach(StructuredGrid & data);
		void updateDensity(StructuredGrid & data);

		//public parameters (the problem parameters)
		double gamma, inlet_pressure, total_temperature, s_star;
		double initial_pressure_left, initial_pressure_right, time;
		double dt = 0.1; //TODO: initialize this to something.


		QuasiEuler(StructuredGrid &mesh, double gamm)
			: local_matrix_size(mesh.problem_dimension + 2), gamma(gamm)
		{
		};
		QuasiEuler operator=(const QuasiEuler& E){return *this;};

};

void QuasiEuler::pressureSensor(StructuredGrid &data){
	//assert(0 && "pressure sensor needs data outside of range, fix this before using.");
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

//void QuasiEuler::updatePressure(StructuredGrid & data){
//	assert(data.Q_has_been_updated && "you need to update Q before updating other variables");
//	for (int i = data.buffer_size; i < data.stop_iteration_index; ++i)
//	{
//		double q1 = data.Q(i);
//		double q2 = data.Q(i+1);
//		double q3 = data.Q(i+2);
//		data.pressure(i) = (gamma-1)*(q3/q1 - 1/(2*q1)*std::pow(q2,2));
//	}
//}
//
//void QuasiEuler::updateMach(StructuredGrid & data){
//	assert(data.Q_has_been_updated && "you need to update Q before updating other variables");
//	for (int i = data.buffer_size; i < data.stop_iteration_index; ++i)
//	{
//		double q1 = data.Q(i);
//		double q2 = data.Q(i+1);
//		double q3 = data.Q(i+2);
//		data.pressure(i) = (gamma-1)*(q3/q1 - 1/(2*q1)*std::pow(q2,2));
//	}
//}
//
//void QuasiEuler::updateDensity(StructuredGrid & data){
//	assert(data.Q_has_been_updated && "you need to update Q before updating other variables");
//	for (int i = data.buffer_size; i < data.stop_iteration_index; ++i)
//	{
//		double q1 = data.Q(i);
//		double q2 = data.Q(i+1);
//		double q3 = data.Q(i+2);
//		data.density(i) = q1;
//	}
//}



















#endif
