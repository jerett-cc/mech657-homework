
#include "Geometry.h"
#include "QuasiEuler.h"
#include "SpatialLinearSystemConstruction.h"
#include <cmath>
#include <iostream>
#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"

/*TODO:
 *    +++1. write conversion functions
 * 		2. add printing functionality for results
 * 		+++3. add extrapolation functions for filling values ad the padded boundaries (linear??)
 * 		+++4. calculate the RHS dissipation, using averaging the matrix values A
 * 		5. change artificial dissipation numbers for the problems
 * 			i. subsonic kappa2 = 0, kappa4 = 0.02
 * 			ii. transonic and shock kappa2 = 0.5, kappa4 = 0.02
 * 			iii. i might not need to be implemented...
 * 			iv. calculate a timestep from equation 4.138 based on a courant number (of about 2?)
 * 		6. add a calculation of numerical error? need to incorporate the other code to calculate exact solutions
 * 			i. use equation 4.176 from book
 * 		+++7. figure out initial conditions (state at the inflow?
 * 		8. calculate a courant number using this limits the time step for the shock problem
 * 			i. u = 300m/s, a = 315m/s
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
		QuasiEuler euler_problem(grid_p1, gamma, total_inlet_pressure_1, total_temperature_1, S_star_1, R, 0);

		//set up test Q
		std::cout << "setting initial condition " << std::endl;
		euler_problem.setInitialCondition(grid_p1);
		std::cout << "q in middle is\n" << grid_p1.Q[4] << std::endl;

//		Eigen::VectorXd test((dimension+2)*(num_cells+1));
//		for (int i=0; i<test.size(); ++i)
//		{
//			test(i) = i;
//		}
		//setup a random Q
//		for (int i=0; i<grid_p1.Q.size(); ++i)
//		{
//  //		std::cout << "size q is " << grid_p1.Q.size() <<std::endl;
//  //		std::cout << "i is " << i <<std::endl;
//      Eigen::Vector3d a = Eigen::Vector3d::Random();
//      a(0) = total_inlet_pressure_1 / (R*total_temperature_1);
//      a(1) = 0.1;
//      a(2) = 0.1;
//      grid_p1.Q[i] += a;
//		}
//		std::cout << "size q is now " << grid_p1.Q.size() <<std::endl;

		SystemConstructionAndSolution solver(grid_p1);
//		//testing, clean up later.
		Eigen::MatrixXd temp;
		std::cout << "testing local inviscid flux calculation "<< std::endl;
		temp = euler_problem.calculateLocalInviscidFluxJacobian(grid_p1, 2);
		temp = euler_problem.calculateLocalInviscidFluxJacobian(grid_p1, 1);

		std::cout << "Testing pressure sensor " << std::endl;
		euler_problem.pressureSensor(grid_p1);
		std::cout <<  grid_p1.sensor_contributions << std::endl;

		Eigen::MatrixXd A(euler_problem.local_matrix_size, euler_problem.local_matrix_size);

		std::cout << "Testing spatial matrix construction "<< std::endl;

		solver.calculateSpatialMatrix(grid_p1, euler_problem);

		std::cout << "testing dissipation matrix " << std::endl;

		solver.calculateDissipationMatrix(grid_p1, euler_problem);

		std::cout<< solver.dense_system_matrix << std::endl;

		std::cout << "testing update quantities " << std::endl;
		euler_problem.updatePhysicalQuantities(grid_p1);
//		grid_p1.clearUpdates();
//		std::cout << "mach is: \n" << grid_p1.mach << std::endl;
//		std::cout << "density is: \n" << grid_p1.density << std::endl;
//		std::cout << "pressure is: \n" << grid_p1.pressure << std::endl;

		std::cout << "testing rhs construction " << std::endl;
		solver.calculateResidualVector(grid_p1, euler_problem);

		std::cout << "residual vector:\n" << solver.system_rhs << std::endl;


		Eigen::BiCGSTAB<Eigen::MatrixXd> solv;

		solv.compute(solver.dense_system_matrix);
		Eigen::VectorXd solution = solv.solve(solver.system_rhs);

		std::cout<< "solution for first time step\n" << solution << std::endl;






		//
	}

return 0;
}
