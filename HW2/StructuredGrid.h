
#include<cassert>
#include<cmath>

#include<iostream>

#include "../Eigen/Dense"


class StructuredGrid{

	public:

		Eigen::VectorXd mesh, temperature, pressure, density, mach, Q;

		StructuredGrid(double start, double stop, int num_cells) :
			L(start), R(stop), numCell(num_cells)
			{
				dx = (R - L)/num_cells;
				Eigen::VectorXd tempMesh(num_cells + 1);
					for (int j = 0; j < tempMesh.size(); ++j)
					{
						tempMesh(j) = L + dx*j;
					}
				mesh = tempMesh;
				assert(std::fabs(mesh(num_cells)-R)<1e-15 && "Mesh should end at the specified stop value.");
			};

	private:
		double L, R, dx;
		int numCell;
};
