/*
 * This class includes all the functions to calculate local fluxes, RHS contributions and all functions needed to calculate
 * physical quantities of interest from the solved quantities, relations include sound speed, temperature, and pressure.
 *
 * this class also senses shocks with a pressure sensor found in "Fundamentals Algorithms in Computational Fluid Dynamics"
 * by Pulliam and Zingg 2014 pp. 97
 */
#include <cmath>

#include "../Eigen/Dense"
#include "Geometry.h"

class QuasiEuler{

	public:

		//public functions
		void S(const StructuredGrid &mesh);//return S(x) at mesh points
			//TODO: write conversion functions and test them, could copy from hw1?
		void calculateLocalInviscidFluxJacobian();//do this for node i
		void calculateLocalDissipation();

		//helper functions
		void pressureSensor(const StructuredGrid & mesh);


		//public parameters (the problem parameters)
		double gamma, inlet_pressure, total_temperature, s_star;
		double initial_pressure_left, initial_pressure_right, time;

};


void QuasiEuler::pressureSensor(const StructuredGrid & mesh){
	//do something
}
