
#include "Geometry.h"
#include "QuasiEuler.h"
#include "SpatialLinearSystemConstruction.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
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

const double tol = 1e-9;

template <typename vecT>
void printError(vecT vec, std::string f_name)
{
  std::ofstream a_file;
  a_file.open(f_name + ".out", std::ios::out | std::ios::trunc);

  for (unsigned int i = 0; i< vec.size(); ++i)
  {
    a_file << vec[i] << std::endl;
  }
  a_file.close();

}

void solveSteadyProblem(ProblemData &data,
                        QuasiEuler &problem,
                        Solver &solver,
                        parameters &param,
                        const unsigned int num_iterations,
                        const std::string fname)
{
    std::vector<double> err;
    data.setInitialCondition(param.gamma, param.inlet_pressure, param.total_temp, param.R, param.s_star);
    unsigned int iteration=0;
    double error = 1;
    while (error > tol && iteration < num_iterations)
   {
     solver.reinit();
     problem.calculateSensorContributions();
     //std:: cout << "contributions ____" << std::endl;
     //std::cout << problem.sensor_contributions << std::endl;
     //std::cout << "____________________________________" << std::endl;
     solver.setupSystem();
     solver.solveSystem();
     data.printQuantities(fname);
     //std::cout << "new q is" << std::endl;
     //std::cout << data.getQVect() << std::endl;
     //std::cout << "________Iteration over__________________________"<< std::endl;
     ++iteration;
     error = solver.L2Error();

     //std::cout << "Delta Q" << '\n';
     //data.updateBV();
     err.push_back(error);
   }
    for(unsigned int i = 0; i < data.q.size(); ++i)
    {
      std::cout << solver.delta_Q[i] << '\n';
    }
    std::cout << "End Error: " << std::setprecision(16) <<solver.L2Error() << "\n"
              << "Num iterations: " << iteration <<  std::endl;
    printError<std::vector<double>>(err, fname+"_error"+std::to_string((int)problem.cfl));
};

void solveTimeDependentProblem(ProblemData &data,
                               QuasiEuler &problem,
                               Solver &solver,
                               parameters &param,
                               std::string fname)
{
    data.setInitialCondition(param.gamma, param.inlet_pressure, param.total_temp, param.R, param.s_star);
    unsigned int iteration=0;
    double error = 1;
    double t = 0;
    while (t < solver.end_time)
   {
     solver.reinit();
     problem.calculateSensorContributions();
     std:: cout << "contributions ____" << std::endl;
     std::cout << problem.sensor_contributions << std::endl;
     std::cout << "____________________________________" << std::endl;
     solver.setupSystem();
     solver.solveSystem();
     data.printQuantities(fname);
     std::cout << "new q is" << std::endl;
     std::cout << data.getQVect() << std::endl;
     std::cout << "________Iteration over__________________________"<< std::endl;
     t += problem.dt;
     ++iteration;
     error = solver.L2Error();
   }
    std::cout << "Num iterations: " << iteration <<  std::endl;
};

//implements a second order finite difference algorithm
int main(){

	//Parameters for all problems
	int dimension = 1;
	double Left = 0.;
	double Right = 10.;
	double gamma = 1.4;
	double R = 287.;
	int num_nodes = 99;
	unsigned int max_num_iterations = 2000;
	double max_vel = 315.;

	//problem 1 additional parameters
	double S_star_1 = 0.8;
	double total_temperature_1 = 300;
	double total_inlet_pressure_1 = 1e5;
    double cfl1 = 1000;

	//problem 2 additional parameters
	double shock_location = 7.;
	double S_star_2 = 1.;
    double cfl2 = 1;

	//problem 3 additional parameters
	double initial_pressure_right_3 = 1e4;
	double initial_pressure_left_3 = 1e5;
	double initial_density_right_3 = 0.125;
	double initial_density_left_3 = 1;
	double end_time_3 = 0.0061;
	double cfl3 = 0.5;

	parameters param1(total_inlet_pressure_1, total_temperature_1, gamma, R, S_star_1);
	parameters param2(total_inlet_pressure_1, total_temperature_1, gamma, R, S_star_2);
	parameters param3(initial_pressure_right_3, 1, gamma, R, S_star_1);

  //FIXME: what is the 61 here??
	ProblemData data1(num_nodes,61, Left, Right, param1);
	ProblemData data2(num_nodes, 1, Left, Right, param2);
	ProblemData data3(num_nodes, 1, Left, Right, param3);

	QuasiEuler problem2(&data2, cfl1, max_vel);
	Solver solver2(&data2, &problem2);

	QuasiEuler problem1(&data1, cfl1, max_vel);
	Solver solver1(&data1, &problem1);

	QuasiEuler problem3(&data3, cfl3, max_vel);
	Solver solver3(&data3, &problem3, end_time_3);

  data1.setInitialCondition(gamma, total_inlet_pressure_1, total_temperature_1, R, S_star_1);
  data2.setInitialCondition(gamma, total_inlet_pressure_1, total_temperature_1, R, S_star_2);

  //data1.printQuantities("initial_problem");
  //data2.printQuantities("initial_problem2");

  problem1.calculateSensorContributions();
  std::cout << "-----------------------------" << std::endl;
  std::cout << "testing pressure sensor" << std::endl;
  std::cout << problem1.sensor_contributions << std::endl;

  std::cout << "-----------------------------" << std::endl;
  std::cout << "testing setupSystem" << std::endl;
  solver1.reinit();
  solver1.setupSystem();

  std::cout << "Testing 100 iterations" << std::endl;
//solveSteadyProblem(data2, problem2, solver2, param2, max_num_iterations, "./problem2figs/problem2");
solveTimeDependentProblem(data3, problem3, solver3, param3, "./problem3figs/problem3");
//solveSteadyProblem(data1, problem1, solver1, param1, max_num_iterations, "./problem1figs/problem1");
return 0;
}
