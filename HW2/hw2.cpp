
#include "Geometry.h"
#include <cmath>
#include <iostream>


//implements a second order finite difference algorithm
int main(){

	//Parameters for all problems
	double Left = 0.;
	double Right = 10.;
	double gamma = 1.4;
	double R = 287.;
	int num_cells = 10;

	//problem 1 additional parameters
	double S_star_1 = 0.8;
	double total_temperature_1 = 300;
	double total_inlet_pressure_1 = 1e5;

	//problem 2 additional parameters
	double shock_location = 7.;
	double S_star_2 = 1.;

	//problem 3 additional parameters
	double initial_pressure_right_3 = 1e4;
	double initial_pressure_left_3 = 1e5;
	double initial_density_right_3 = 0.125;
	double initial_density_left_3 = 1;
	double end_time_3 = 0.0061;


	//outer time loop will go here. for the first two setups, we will solve a steady state
	//for the last, we will solve a time dependent problem up from t= [0,6.1ms]

	{
		StructuredGrid grid_p1(Left, Right, num_cells);
		//QuasiEuler euler_problem;

		//
	}

return 0;
}
