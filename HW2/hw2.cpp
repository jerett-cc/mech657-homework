
#include "Geometry.h"
#include "QuasiEuler.h"
#include "SpatialLinearSystemConstruction.h"
#include <cmath>
#include <iostream>
#include "../Eigen/Dense"

/*TODO: either pad the vectors, or create a conditional in loops that require data
 * out of bounds for the mesh size. 
 */


//implements a second order finite difference algorithm
int main(){

	//Parameters for all problems
	int dimension = 1;
	double Left = 0.;
	double Right = 10.;
	double gamma = 1.4;
	double R = 287.;
	int num_cells = 2;

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
		StructuredGrid grid_p1(Left, Right, num_cells, dimension);

		//set up test Q
		Eigen::VectorXd test((dimension+2)*(num_cells+1));
		for (int i=0; i<test.size(); ++i)
		{
			test(i) = i;
		}
		//setup a random Q
		for (int i=0; i<grid_p1.Q.size(); ++i)
		{
//		std::cout << "size q is " << grid_p1.Q.size() <<std::endl;
//		std::cout << "i is " << i <<std::endl;
		Eigen::Vector3d a = Eigen::Vector3d::Random();
		if (std::fabs(a(0)-0)<1e-14)
		{
			a(0) = 1.;
		}
		else if (std::fabs(a(1)-0)<1e-14)
		{
			a(1) = 1.;
		}
		grid_p1.Q[i] += a;
		}
//		std::cout << "size q is now " << grid_p1.Q.size() <<std::endl;
		QuasiEuler euler_problem(grid_p1, gamma);
		SystemConstructionAndSolution solver(grid_p1);
//		//testing, clean up later.
		Eigen::MatrixXd temp;
		std::cout << "testing local inviscid flux calculation "<< std::endl;
		temp = euler_problem.calculateLocalInviscidFluxJacobian(grid_p1, 2);
		temp = euler_problem.calculateLocalInviscidFluxJacobian(grid_p1, 1);

		for (int i=0; i<grid_p1.pressure.size(); ++i)
				{
					if (i<=grid_p1.pressure.size()/2)
					{
						grid_p1.pressure(i) = 1.;
					}
					else {
						grid_p1.pressure(i) = 0.;
					}
				}
		std::cout << "Testing pressure sensor " << std::endl;
		euler_problem.pressureSensor(grid_p1);
		std::cout <<  grid_p1.sensor_contributions << std::endl;

		Eigen::MatrixXd A(euler_problem.local_matrix_size, euler_problem.local_matrix_size);

		std::cout << "Testing spatial matrix construction "<< std::endl;

		solver.calculateSpatialMatrix(grid_p1, euler_problem);

		std::cout << "testing dissipation matrix " << std::endl;

		solver.calculateDissipationContributionToMatrix(grid_p1);

//		std::cout<< solver.dense_system_matrix << std::endl;







		//
	}

return 0;
}
