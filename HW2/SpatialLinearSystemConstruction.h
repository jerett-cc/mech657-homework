
#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"
#include "../Eigen/Sparse"
#include "QuasiEuler.h"
#include "Geometry.h"


class SystemConstructionAndSolution{

	public:

		//constructor
		SystemConstructionAndSolution(StructuredGrid & mesh)
			:system_matrix(3*mesh.get_size(), 3*mesh.get_size()), system_rhs(3*mesh.get_size())
		{

		}
		//public functions
		void constructSystemMatrix();//do this block row by block row?
		void constructRHS();
		void calculateResidualVector();//this is R(Q), do this point by point and

			//TODO: add functions for the linear solve step

		//public variables
		Eigen::SparseMatrix<double, 1> system_matrix;//row major system matrix
		Eigen::VectorXd system_rhs;//rhs of system


};


