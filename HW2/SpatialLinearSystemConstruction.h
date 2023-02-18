#ifndef SpatialLinearSystemCons_H
#define SpatialLinearSystemCons_H


#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"
#include "../Eigen/Sparse"
#include "QuasiEuler.h"
#include "Geometry.h"

typedef Eigen::Triplet<double> T;

class SystemConstructionAndSolution{

	public:

		//constructor
		SystemConstructionAndSolution(StructuredGrid & mesh)
			:system_matrix(3*mesh.get_size(), 3*mesh.get_size()), system_rhs(3*mesh.get_size())
		{
		}
		//public functions
		void constructSystemMatrix(const StructuredGrid data, const QuasiEuler &problem);//do this block row by block row?
		void constructRHS();
		void calculateResidualVector();//this is R(Q), do this point by point and

			//TODO: add functions for the linear solve step

		//public variables
		Eigen::SparseMatrix<double, 1> system_matrix;//row major system matrix
		Eigen::VectorXd system_rhs;//rhs of system
		T identity;//TODO: make a triplet list of identity elements, either here or in
};


//getting an error here I do not understand
void SystemConstructionAndSolution::constructSystemMatrix(const StructuredGrid data, const QuasiEuler &problem){
	//iterate over all nodes
	std::vector<T> coefficients;
	coefficients.reserve(25*(data.numCell+1));//reserve an estimate of number of nonzero entries

	for (int i = data.buffer_size; i<data.stop_iteration_index; ++i)
	{
		Eigen::MatrixXd A(problem.local_matrix_size, problem.local_matrix_size);
//		A = problem.calculateLocalInviscidFluxJacobian(data, );//error here??
		/*
		 * calculate local A_i,A_i-1, A+1, then add to a triplet list by looping through the proper indices
		 * coefficients.pushBack(T(i,j,entry_ij)
		 *
		 * then when the list of triplets is done, do
		 * system_matrix.setFromTriplets(coeffifients.begin(), coefficients.end())
		 *
		 *
		 *need to include a calculation of the linearization L of the artificial dissipation
		 */
	}
}

#endif
