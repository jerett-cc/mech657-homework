#ifndef GEOMETRY_H
#define GEOMETRY_H


#include<cassert>
#include<cmath>

#include<iostream>

#include "../Eigen/Dense"
//TODO: include a way to describe the face of a cell once we start doing cell-cenered Finite Volume code

class StructuredGrid{

	public:

		const unsigned int problem_dimension;

		Eigen::VectorXd mesh, temperature, pressure, density, mach, Q;
//TODO: write this in a more dimension independent way, num cells to node points,
// 		also, this class is only good for finite difference, make it work for finite volume?
		StructuredGrid(double start, double stop, int num_cells, int dim) :
			L(start), R(stop), numCell(num_cells), mesh(num_cells+1), problem_dimension(dim)
			, Q((num_cells+1)*(problem_dimension+2)), pressure((num_cells+1)*(problem_dimension+2))
			{
				dx = (R - L)/num_cells;
//				Eigen::VectorXd tempMesh(num_cells + 1);
					for (int j = 0; j < mesh.size(); ++j)
					{
						mesh(j) = L + dx*j;
					}
				assert(std::fabs(mesh(num_cells)-R)<1e-15 &&
						"Mesh should end at the specified stop value.");
			};

		int get_size();


	private:
		double L, R, dx;
		int numCell;
};

int StructuredGrid::get_size(){
	return mesh.size();
}

#endif
