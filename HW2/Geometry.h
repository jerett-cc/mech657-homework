#ifndef GEOMETRY_H
#define GEOMETRY_H


#include<cassert>
#include<cmath>

#include<iostream>

#include "../Eigen/Dense"
//TODO: include a way to describe the face of a cell once we start doing cell-cenered Finite Volume code

class StructuredGrid{

	public:

		const unsigned int problem_dimension, buffered_length;
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
			buffered_length((num_cells+1+buffer_size)*(problem_dimension + 2)), Q((num_cells+1+buffer_size)*(problem_dimension + 2)),
			pressure((num_cells+1+buffer_size)*(problem_dimension + 2)), density((num_cells+1+buffer_size)*(problem_dimension + 2)),
			temperature((num_cells+1+buffer_size)*(problem_dimension + 2)), mach((num_cells+1+buffer_size)*(problem_dimension + 2))
			{
				dx = (R - L)/num_cells;
				stop_iteration_index = buffered_length - buffer_size;
//				Eigen::VectorXd tempMesh(num_cells + 1);
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
