#ifndef GEOMETRY_H
#define GEOMETRY_H


#include<cassert>
#include<cmath>
#include<vector>

#include<iostream>

#include "../Eigen/Dense"
//TODO: include a way to describe the face of a cell once we start doing cell-cenered Finite Volume code

class StructuredGrid{

	public:

	    int problem_dimension, buffered_length;
		int stop_iteration_index;
		const int buffer_size = 2;
		double dx;
		bool Q_has_been_updated = 1;
		bool E_has_been_updated = 0;
		bool Pressure_has_been_updated = 0;
		bool Density_has_been_updated = 0;
		bool Temperature_has_been_updated = 0;
		bool boundary_has_been_updated = 0;

		Eigen::VectorXd mesh, temperature, pressure, density, mach;
		std::vector<Eigen::Vector3d> Q, E;
		Eigen::MatrixXd sensor_contributions;
//TODO: write this in a more dimension independent way, num cells to node points,
// 		also, this class is only good for finite difference, make it work for finite volume?

//TODO: fix readability of constructor
		StructuredGrid(double start, double stop, int num_cells, int dim) :
			L(start), R(stop), numCell(num_cells), mesh(num_cells+1), problem_dimension(dim), num_node(num_cells+1),
			Q(num_cells+1+2*buffer_size), E(num_cells+1+2*buffer_size), pressure(num_cells+1+2*buffer_size), density(num_cells+1+2*buffer_size),
			temperature(num_cells+1+2*buffer_size), mach(num_cells+1+2*buffer_size), sensor_contributions(mesh.size(), 2)
			{
				dx = (R - L)/num_cells;
				buffered_length = (num_cells+1+2*buffer_size)*(problem_dimension + 2);//do I need still?
//				std::cout << "buff len " << buffered_length << std::endl;
				stop_iteration_index = buffer_size + num_cells + 1;//TODO: rename this across project to something like node iterator stop_index
				num_components = dim +2;
					for (int j = 0; j < mesh.size(); ++j)
					{
						mesh(j) = L + dx*j;
					}
				assert(std::fabs(mesh(numCell)-R)<1e-15 &&
						"Mesh should end at the specified stop value.");
			};

		int get_size()const ;
		int get_system_size() const;
		int estimateNumberNonzeroElements();
		void clearUpdates();
		void interpolateBoundaryE();
		void interpolateBoundaryQ();
		void interpolateBoundaryQuantities();

		double L, R;
		int numCell;
		int num_node;
		int num_components;
};

int StructuredGrid::get_size() const{
	return mesh.size()*(problem_dimension+2);
}

int StructuredGrid::get_system_size() const{
  return (mesh.size()-buffer_size-2) *(problem_dimension + 2);
}

//TODO: decide if this function is defunct
int StructuredGrid::estimateNumberNonzeroElements(){
	return 25*(numCell+1);
}

//TODO setup the Q update flag
void StructuredGrid::clearUpdates(){
//  Q_has_been_updated = 0;
  E_has_been_updated = 0;
  Pressure_has_been_updated = 0;
  Density_has_been_updated = 0;
  Temperature_has_been_updated = 0;
}

//interpolates the boundary values to use for our stencil
void StructuredGrid::interpolateBoundaryE(){
  for (int i = buffer_size-1; i>-1; --i)
    {
//    std::cout <<"boundary index: " <<  i<< std::endl;
    Eigen::Vector3d left_slope_E = (E[i + 2] - E[i + 1])/dx;
//    std::cout << "Left slope \n" << left_slope<< std::endl;
    Eigen::Vector3d right_slope_E = (E[E.size()-1-i-1] - E[E.size()-1-i-2])/dx;
//    std::cout << "R slope \n" << right_slope<< std::endl;
    E[i] = E[i+1] - left_slope_E*dx;
//    std::cout << "left value \n" << Q[i]<< std::endl;
    E[E.size()-1-i] = E[E.size()-1-i-1] + right_slope_E*dx;
//    std::cout << "right value \n" << Q[Q.size()-1-i]<< std::endl;
    }

}

void StructuredGrid::interpolateBoundaryQ(){
  for (int i = buffer_size-1; i>-1; --i)
  {
//    std::cout <<"boundary index: " <<  i<< std::endl;
    Eigen::Vector3d left_slope = (Q[i + 2] - Q[i + 1])/dx;
//    std::cout << "Left slope \n" << left_slope<< std::endl;
    Eigen::Vector3d right_slope = (Q[Q.size()-1-i-1] - Q[Q.size()-1-i-2])/dx;
//    std::cout << "R slope \n" << right_slope<< std::endl;
    Q[i] = Q[i+1] - left_slope*dx;
//    std::cout << "left value \n" << Q[i]<< std::endl;
    Q[Q.size()-1-i] = Q[Q.size()-1-i-1] + right_slope*dx;
//    std::cout << "right value \n" << Q[Q.size()-1-i]<< std::endl;

  }
//  std::cout<< "----------------------------" << std::endl;
//  boundary_has_been_updated = 1;

//  std::cout << Q[0] << Q[1] << std::endl;
}

void StructuredGrid::interpolateBoundaryQuantities(){
  //interpolate pressure
  for (int i = buffer_size-1; i>-1; --i)
    {
    double left_slope_P = (pressure(i + 2) - pressure(i + 1))/dx;
    double right_slope_P = (pressure(pressure.size()-1-i-1) - pressure(pressure.size()-1-i-2))/dx;
    pressure(i) = pressure(i+1) - left_slope_P*dx;
    pressure(pressure.size()-1-i) = pressure(pressure.size()-1-i-1) + right_slope_P*dx;
    }

  //interpolate density
  for (int i = buffer_size-1; i>-1; --i)
  {
    double left_slope_D = (density(i + 2) - density(i + 1))/dx;
    double right_slope_D = (density(density.size()-1-i-1) - density(density.size()-1-i-2))/dx;
    density(i) = density(i+1) - left_slope_D*dx;
    density(density.size()-1-i) = density(density.size()-1-i-1) + right_slope_D*dx;
  }
  //interpolate mach
  for (int i = buffer_size-1; i>-1; --i)
  {
    double left_slope_M = (mach(i + 2) - mach(i + 1))/dx;
    double right_slope_M = (mach(mach.size()-1-i-1) - mach(mach.size()-1-i-2))/dx;
    mach(i) = mach(i+1) - left_slope_M*dx;
    mach(mach.size()-1-i) = mach(mach.size()-1-i-1) + right_slope_M*dx;
  }
    //todo add temp?
}








































#endif
