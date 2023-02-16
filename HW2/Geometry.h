
#include<cassert>
#include<cmath>

#include<iostream>

#include "../Eigen/Dense"
//TODO: include a way to describe the face of a cell once we start doing cell-cenered Finite Volume code

class StructuredGrid{

	public:

		Eigen::VectorXd mesh, temperature, pressure, density, mach, Q;

		StructuredGrid(double start, double stop, int num_cells) :
			L(start), R(stop), numCell(num_cells), mesh(num_cells+1)
			{
				dx = (R - L)/num_cells;
//				Eigen::VectorXd tempMesh(num_cells + 1);
					for (int j = 0; j < mesh.size(); ++j)
					{
						mesh(j) = L + dx*j;
					}
				assert(std::fabs(mesh(num_cells)-R)<1e-15 && "Mesh should end at the specified stop value.");
				std::cout << mesh << std::endl;
			};

		int get_size();


	private:
		double L, R, dx;
		int numCell;
};

int StructuredGrid::get_size(){
	return mesh.size();
}
