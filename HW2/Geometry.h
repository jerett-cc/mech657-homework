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
		bool Q_has_been_updated = 0;

		Eigen::VectorXd mesh, temperature, pressure, density, mach;
		std::vector<Eigen::Vector3d> Q;
		Eigen::MatrixXd sensor_contributions;
//TODO: write this in a more dimension independent way, num cells to node points,
// 		also, this class is only good for finite difference, make it work for finite volume?

//TODO: fix readability of constructor
		StructuredGrid(double start, double stop, int num_cells, int dim) :
			L(start), R(stop), numCell(num_cells), mesh(num_cells+1), problem_dimension(dim), num_node(num_cells+1),
			Q(num_cells+1), pressure(num_cells+1), density(num_cells+1), temperature(num_cells+1),
			mach(num_cells+1), sensor_contributions(mesh.size(), 2)
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
		int estimateNumberNonzeroElements();

		double L, R;
		int numCell;
		int num_node;
		int num_components;
};

int StructuredGrid::get_size() const{
	return mesh.size()*(problem_dimension+2);
}

int StructuredGrid::estimateNumberNonzeroElements(){
	return 25*(numCell+1);
}


#endif
