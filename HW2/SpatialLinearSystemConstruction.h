
#include "../Eigen/Dense"
#include "../Eigen/IterativeLinearSolvers"

class SystemConstruction{

	public:
		//public functions
		void constructSpatialMatrix();//do this row by row?
		void calculateLocalFluxJacobian();//do this for node i
		void calculateResidualVector();//this is R(Q)
			//TODO: add functions for the linear solve step

		//public parameters
};
