
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
	int num_nodes = 6;

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


	// step 1: get initial conditions and initialize Q to this taking into account S(x)

	// step 2: calculate A = I + delta_t*dxA - delta_t*L
	//  this requires A and L

	// step 3: calculate RHS = -delta_t*E + delta_t*Dx + BC -->does this need to be multiplied by anything?
	//  this requires a way to calculate E from Q, and a way to caculate Dx --> this needs ghost cells??
	//  we also need to describe the boundary conditions and add them to the right piece of the RHS

	// step 4: solve the system
	parameters param1(total_inlet_pressure_1, total_temperature_1, gamma, R, S_star_1);
	parameters param2(total_inlet_pressure_1, total_temperature_1, gamma, R, S_star_2);

	ProblemData data1(num_nodes, 1, Left, Right, param1);
//	ProblemData data2(num_nodes, 1, Left, Right, param2);

	QuasiEuler problem1(&data1, 0.5);
	Solver solver(&data1, &problem1);

	data1.setInitialCondition(gamma, total_inlet_pressure_1, total_temperature_1, R, S_star_1);
//  std::cout << "-----------------------------" << std::endl;
//	std::cout << data1[-2] << std::endl;
//  std::cout << "-----------------------------" << std::endl;
//  std::cout << data1[0] << std::endl;
//  std::cout << "-----------------------------" << std::endl;
//  std::cout << data1[3] << std::endl;
//  std::cout << "-----------------------------" << std::endl;
//  std::cout << "testing get_Q " << std::endl;
//  std::cout << data1.getQVect() << std::endl;
//  std::cout << "-----------------------------" << std::endl;
//  std::cout << "testing pressure at L endpoint match initial pressure?" << std::endl;
//  std::cout << data1.Pressure(-1) << std::endl;

  data1.printQuantities("initial_problem");
//
//  problem1.calculateSensorContributions();
//  std::cout << "-----------------------------" << std::endl;
//  std::cout << "testing pressure sensor, expect zero except maybe at end" << std::endl;
//  std::cout << problem1.sensor_contributions << std::endl;
//
//  std::cout << "-----------------------------" << std::endl;
//  std::cout << "testing local flux jacobian" << std::endl;
//  std::cout << problem1.calculateLocalFluxJacobian(0) << std::endl;

//  std::cout << "-----------------------------" << std::endl;
//  std::cout << "testing setup A " << std::endl;
//  solver.setup_A();

//  std::cout << "-----------------------------" << std::endl;
//  std::cout << "testing setup b " << std::endl;
//  solver.calcDe();
//  solver.calcDx();

   solver.reinit();


  //solve the system
  std::cout << "-----------------------------" << std::endl;
  std::cout << "testing solve system " << std::endl;
  problem1.calculateSensorContributions();
  solver.setupSystem();
  solver.solveSystem();
  solver.reinit();

  std::cout << "-----------------------------" << std::endl;
  std::cout << "new q is" << std::endl;
  std::cout << data1.getQVect() << std::endl;

  data1.printQuantities("one step");


  std::cout << "-----------------------------" << std::endl;
  std::cout << "testing solve system step 2 " << std::endl;
  problem1.calculateSensorContributions();
  solver.setupSystem();
  solver.solveSystem();
  solver.reinit();
//
//
//  std::cout << "------------Test 4 iterations----------------" << std::endl;
//  std::cout << "new q at step 2 is" << std::endl;
//  std::cout << data1.getQVect() << std::endl;
//
//  data1.printQuantities("two step");
//
//   for (int i = 0; i<4;++i)
//   {
//     problem1.calculateSensorContributions();
//     solver.setupSystem();
//     solver.solveSystem();
//     std::cout << "new q is" << std::endl;
//     std::cout << data1.getQVect() << std::endl;
//     std::cout << "________Iteration over__________________________"<< std::endl;
//     solver.reinit();
//   }


return 0;
}
