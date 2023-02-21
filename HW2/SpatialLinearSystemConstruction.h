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
			:system_matrix(mesh.get_size(), mesh.get_size()), system_rhs(mesh.get_size()),
			 sensor_factors(mesh.num_node, 2)
		{
			dense_system_matrix = Eigen::MatrixXd::Identity(mesh.get_size(), mesh.get_size());
			std::cout << dense_system_matrix<<std::endl;
		}
		//public functions
		void constructSystemMatrix(StructuredGrid &data, QuasiEuler &problem);//do this block row by block row?
		void constructRHS();
		void calculateResidualVector();//this is R(Q), do this point by point and
		void calculateSpatialMatrix(StructuredGrid &data, const QuasiEuler &problem);
		void calculateDissipationContributionToMatrix();

			//TODO: add functions for the linear solve step

		//public variables
		Eigen::SparseMatrix<double, 1> system_matrix;//row major system matrix
		Eigen::MatrixXd dense_system_matrix;
		Eigen::VectorXd system_rhs;//rhs of system
		Eigen::MatrixXd sensor_factors;
		T identity;//TODO: make a triplet list of identity elements, either here or in
};


void SystemConstructionAndSolution::calculateSpatialMatrix(StructuredGrid &data, const QuasiEuler &problem){

	Eigen::MatrixXd Ai(problem.local_matrix_size, problem.local_matrix_size);
	Eigen::MatrixXd Ai_n(problem.local_matrix_size, problem.local_matrix_size);

	for (int i = data.buffer_size; i <data.stop_iteration_index-1; ++i)//loop through all internal nodes
	{
		Ai = problem.calculateLocalInviscidFluxJacobian(data, i);//calculate flux at node i
		Ai_n = problem.calculateLocalInviscidFluxJacobian(data, i+1);//calculate flux at node i+1
		int sub_index = problem.local_matrix_size * (i - data.buffer_size);

		std::cout << "node " << sub_index << " , " << sub_index+problem.local_matrix_size << std::endl;
		dense_system_matrix.block<3,3>(sub_index, sub_index+problem.local_matrix_size) += 0.5*problem.dt*Ai_n;
		dense_system_matrix.block<3,3>(sub_index + problem.local_matrix_size, sub_index) += -0.5*problem.dt*Ai;
	}

	std::cout << "System matrix excluding first and last " << std::endl << dense_system_matrix << std::endl;

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
