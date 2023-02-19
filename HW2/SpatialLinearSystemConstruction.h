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
		void constructSystemMatrix(StructuredGrid &data, QuasiEuler &problem);//do this block row by block row?
		void constructRHS();
		void calculateResidualVector();//this is R(Q), do this point by point and

			//TODO: add functions for the linear solve step

		//public variables
		Eigen::SparseMatrix<double, 1> system_matrix;//row major system matrix
		Eigen::VectorXd system_rhs;//rhs of system
		T identity;//TODO: make a triplet list of identity elements, either here or in
};


//getting an error here I do not understand
void SystemConstructionAndSolution::constructSystemMatrix(StructuredGrid &data, QuasiEuler &problem){
	//TODO: handle the first actual node, and the last, separately outside of the loop?
	std::vector<T> coefficients;
	coefficients.reserve(25*(data.numCell+1));//reserve an estimate of number of nonzero entries
	Eigen::MatrixXd Ai(problem.local_matrix_size, problem.local_matrix_size);
	Eigen::MatrixXd Ai_p(problem.local_matrix_size, problem.local_matrix_size);
	Eigen::MatrixXd Ai_n(problem.local_matrix_size, problem.local_matrix_size);

	for (int i = data.buffer_size; i<data.stop_iteration_index; ++i)
	{

//		assert(i>=0 && i<problem.local_matrix_size && "constructSystemMatrix problem");

		Ai = problem.calculateLocalInviscidFluxJacobian(data, i);//TODO: pass by reference, const?
		Ai_p = problem.calculateLocalInviscidFluxJacobian(data, i-1);
		Ai_n = problem.calculateLocalInviscidFluxJacobian(data, i+1);
		std::cout << "____________________Previous______________________________" << std::endl;
		std::cout << Ai_p << std::endl;
		std::cout << "____________________Current______________________________" <<std::endl;
		std::cout << Ai << std::endl;
		std::cout << "____________________Next______________________________" <<std::endl;
		std::cout << Ai_n << std::endl;
		std::cout << "__________________________________________________" <<std::endl;
//		std::cout << data.stop_iteration_index << " stop " << std::endl;

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
