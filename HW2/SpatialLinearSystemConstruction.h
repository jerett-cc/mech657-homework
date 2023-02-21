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
			 sensor_factors(mesh.num_node, 2), system_size(mesh.get_size())
		{
			dense_system_matrix = Eigen::MatrixXd::Identity(mesh.get_size(), mesh.get_size());
			identity = Eigen::MatrixXd::Identity(mesh.num_components, mesh.num_components);

		}
		//public functions
		void constructSystemMatrix(StructuredGrid &data, QuasiEuler &problem);//do this block row by block row?
		void constructRHS();
		void calculateResidualVector();//this is R(Q), do this point by point and
		void calculateSpatialMatrix(StructuredGrid &data, const QuasiEuler &problem);
		void calculateDissipationContributionToMatrix(const StructuredGrid &data);

			//TODO: add functions for the linear solve step

		//public variables
		Eigen::SparseMatrix<double, 1> system_matrix;//row major system matrix
		Eigen::MatrixXd dense_system_matrix;
		Eigen::VectorXd system_rhs;//rhs of system
		Eigen::MatrixXd sensor_factors;
		Eigen::MatrixXd identity;

		const int system_size;

};


void SystemConstructionAndSolution::calculateSpatialMatrix(StructuredGrid &data, const QuasiEuler &problem){

	Eigen::MatrixXd Ai(problem.local_matrix_size, problem.local_matrix_size);
	Eigen::MatrixXd Ai_n(problem.local_matrix_size, problem.local_matrix_size);

	for (int i = data.buffer_size; i <data.stop_iteration_index-1; ++i)//loop through all internal nodes
	{
		Ai = problem.calculateLocalInviscidFluxJacobian(data, i);//calculate flux at node i
		Ai_n = problem.calculateLocalInviscidFluxJacobian(data, i+1);//calculate flux at node i+1

		int sub_index = problem.local_matrix_size * (i - data.buffer_size);

		dense_system_matrix.block<3,3>(sub_index, sub_index+problem.local_matrix_size) += 0.5*problem.dt*Ai_n;
		dense_system_matrix.block<3,3>(sub_index + problem.local_matrix_size, sub_index) += -0.5*problem.dt*Ai;
	}

//	std::cout << "System matrix excluding first and last " << std::endl << dense_system_matrix << std::endl;

}

void SystemConstructionAndSolution::calculateDissipationContributionToMatrix(const StructuredGrid &data){

	Eigen::MatrixXd stencil_high_order(data.num_components*data.num_node, data.Q.size());
	Eigen::MatrixXd stencil_low_order(data.num_components*data.num_node, data.Q.size());
	Eigen::MatrixXd result(data.get_size(), data.get_size());

	result = Eigen::MatrixXd::Zero(data.get_size(), data.get_size());

	for (int i = 0; i < data.num_node; ++i)//loop through each node
	{
		stencil_high_order.block<3,3>(i*data.num_components, i*data.num_components) = data.sensor_contributions(i,1)*(1*identity);
		stencil_high_order.block<3,3>(i*data.num_components, (i+1)*data.num_components) = data.sensor_contributions(i,1)*(-4*identity);
		stencil_high_order.block<3,3>(i*data.num_components, (i+2)*data.num_components) = data.sensor_contributions(i,1)*(6*identity);
		stencil_high_order.block<3,3>(i*data.num_components, (i+3)*data.num_components) = data.sensor_contributions(i,1)*(-4*identity);
		stencil_high_order.block<3,3>(i*data.num_components, (i+4)*data.num_components) = data.sensor_contributions(i,1)*(1*identity);
	}

	for (int i = 0; i < data.num_node; ++i)//loop through each node
		{
			stencil_low_order.block<3,3>(i*data.num_components, (i+1)*data.num_components) = data.sensor_contributions(i,0)*(-1*identity);
			stencil_low_order.block<3,3>(i*data.num_components, (i+2)*data.num_components) = data.sensor_contributions(i,0)*(2*identity);
			stencil_low_order.block<3,3>(i*data.num_components, (i+3)*data.num_components) = data.sensor_contributions(i,0)*(-1*identity);
		}

	stencil_high_order += stencil_low_order;

//	std::cout << stencil_high_order << "\n stencil ^^" <<std::endl;

	//extract correct matrix
	for (int i=data.num_components*data.buffer_size; i < data.stop_iteration_index*data.num_components; ++i)
	{
		result.col(i - data.num_components*data.buffer_size) = stencil_high_order.col(i);
	}
	result = 1/data.dx * result;
//	std::cout << result << "\n stencil from e4 contribution ^^" <<std::endl;

	dense_system_matrix += result;

}


#endif
