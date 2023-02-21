
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
	int num_cells = 4;

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
		std::cout << "size q is " << grid_p1.Q.size() <<std::endl;
		grid_p1.Q = Eigen::VectorXd::Random(grid_p1.Q.size(), 1);
		std::cout << "size q is now " << grid_p1.Q.size() <<std::endl;
		QuasiEuler euler_problem(grid_p1, gamma);
		SystemConstructionAndSolution solver(grid_p1);
		//testing, clean up later.
		Eigen::MatrixXd temp;
		temp = euler_problem.calculateLocalInviscidFluxJacobian(grid_p1, 2);
		temp = euler_problem.calculateLocalInviscidFluxJacobian(grid_p1, 1);

		Eigen::MatrixXd e_contributions(grid_p1.mesh.size(), 2);

		for (int i=0; i<grid_p1.pressure.size(); ++i)
				{
					if (i<=grid_p1.pressure.size()/2) {grid_p1.pressure(i) = 0.1;}
					else {
						grid_p1.pressure(i) = 0.;
					}
				}
		euler_problem.pressureSensor(grid_p1, e_contributions);

		std::cout << "Testing pressure sensor " << std::endl <<  e_contributions << std::endl;
		std::cout << "Testing pressure sensor @@@ " << std::endl <<  grid_p1.sensor_contributions << std::endl;

		Eigen::MatrixXd A(euler_problem.local_matrix_size, euler_problem.local_matrix_size);



//		std::cout << "Testing system matrix construction "<< std::endl;
//
//		solver.calculateSpatialMatrix(grid_p1, euler_problem);


		//
	}

return 0;
}
