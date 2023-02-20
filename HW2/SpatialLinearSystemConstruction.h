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
			:system_matrix(mesh.get_size(), mesh.get_size()), system_rhs(mesh.get_size())
		{
		}
		//public functions
		void constructSystemMatrix(StructuredGrid &data, QuasiEuler &problem);//do this block row by block row?
		void constructRHS();
		void calculateResidualVector();//this is R(Q), do this point by point and
		void calculateSpatialMatrix(StructuredGrid &data, QuasiEuler &problem);
		void calculateDissipationContributionToMatrix();

			//TODO: add functions for the linear solve step

		//public variables
		Eigen::SparseMatrix<double, 1> system_matrix;//row major system matrix
		Eigen::VectorXd system_rhs;//rhs of system
		T identity;//TODO: make a triplet list of identity elements, either here or in
};


//getting an error here I do not understand
void SystemConstructionAndSolution::calculateSpatialMatrix(StructuredGrid &data, QuasiEuler &problem){
	//TODO: handle the boundary nodes separately outside of the loop?

	std::vector<T> coefficients;
	coefficients.reserve(data.estimateNumberNonzeroElements());//reserve an estimate of number of nonzero entries
	Eigen::MatrixXd Ai(problem.local_matrix_size, problem.local_matrix_size);
	Eigen::MatrixXd Ai_p(problem.local_matrix_size, problem.local_matrix_size);
	Eigen::MatrixXd Ai_n(problem.local_matrix_size, problem.local_matrix_size);

	for (int i = data.buffer_size+1; i<data.stop_iteration_index-1; ++i)//iterate over all internal nodes
	{

//		assert(i>=0 && i<problem.local_matrix_size && "constructSystemMatrix problem");
//		std::cout << "iteration___________________________________________________________"<< std::endl;
		Ai = problem.calculateLocalInviscidFluxJacobian(data, i);
		Ai_p = problem.calculateLocalInviscidFluxJacobian(data, i-1);
		Ai_n = problem.calculateLocalInviscidFluxJacobian(data, i+1);

		int position_actual = problem.local_matrix_size * (i - data.buffer_size);

//		std::cout << position_actual << " block top " << std::endl;

		for (int j=0; j<problem.local_matrix_size; ++j)//iterate over all components at each node
		{
			coefficients.push_back(T(position_actual+j,position_actual+j,1));//create identity
			for (int k = 0; k < problem.local_matrix_size; ++k)
			{
//				std::cout << "i: " << position_actual << " j: " << position_actual + j << " k: " << k+position_actual + problem.local_matrix_size << std::endl;
				coefficients.push_back(T(position_actual + j, k+position_actual + problem.local_matrix_size, 1./2.*problem.dt*Ai_n(j,k)));
				coefficients.push_back(T(position_actual + j, k+position_actual - problem.local_matrix_size, -1./2.*problem.dt*Ai_p(j,k)));
			}
		}
		/*
		 *need to include a calculation of the linearization L of the artificial dissipation
		 *need to include
		 */
	}
	system_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
	std::cout << "System matrix excluding first and last " << std::endl << Eigen::MatrixXd(system_matrix) << std::endl;

}

void SystemConstructionAndSolution::calculateDissipationContributionToMatrix(){
	std::vector<T> coefficients;
	coefficients.reserve(system_matrix.cols()*5);//reserve an estimate of number of nonzero entries
	//TODO: need to include the coefficients to the pressure sensor here. do both individually, and add them to the main marix
	for (int i = 0; i <system_matrix.rows(); ++i)
	{
		for (int j = 0; j<system_matrix.cols(); ++j)
		{
			coefficients.push_back(T(i,j, 1));//TODO: remove this placeholder, make pentadiagonal logic.
		}
	}

		system_matrix.setFromTriplets(coefficients.begin(), coefficients.end());
		std::cout << "System matrix excluding first and last " << std::endl << Eigen::MatrixXd(system_matrix) << std::endl;
}

#endif
