#ifndef GEOMETRY_H
#define GEOMETRY_H


#include<cassert>
#include<cmath>

#include<iostream>

#include "../Eigen/Dense"
//TODO: include a way to describe the face of a cell once we start doing cell-cenered Finite Volume code

class StructuredGrid{

	public:

	    int problem_dimension, buffered_length;
		int stop_iteration_index;
		const int buffer_size = 2;
		double dx;
		bool Q_has_been_updated = 0;

		Eigen::VectorXd mesh;
		Eigen::VectorXd temperature, pressure, density, mach, Q;
//TODO: write this in a more dimension independent way, num cells to node points,
// 		also, this class is only good for finite difference, make it work for finite volume?

//TODO: fix readability of constructor
		StructuredGrid(double start, double stop, int num_cells, int dim) :
			L(start), R(stop), numCell(num_cells), mesh(num_cells+1), problem_dimension(dim),
			Q((num_cells+1+2*buffer_size)*(problem_dimension + 2)), pressure((num_cells+1+2*buffer_size)*(problem_dimension + 2)),
			density((num_cells+1+2*buffer_size)*(problem_dimension + 2)), temperature((num_cells+1+2*buffer_size)*(problem_dimension + 2)),
			mach((num_cells+1+2*buffer_size)*(problem_dimension + 2))
			{
				dx = (R - L)/num_cells;
				buffered_length = (num_cells+1+2*buffer_size)*(problem_dimension + 2);
//				std::cout << "buff len " << buffered_length << std::endl;
				stop_iteration_index = buffer_size + num_cells + 1;//TODO: rename this across project to something like node iterator stop_index

					for (int j = 0; j < mesh.size(); ++j)
					{
						mesh(j) = L + dx*j;
					}
				assert(std::fabs(mesh(numCell)-R)<1e-15 &&
						"Mesh should end at the specified stop value.");
			};

		int get_size();

		double L, R;
		int numCell;
};

int StructuredGrid::get_size(){
	return mesh.size();
}

#endif
