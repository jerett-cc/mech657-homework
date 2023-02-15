

#include <cmath>

#include "../Eigen/Dense"
#include "Geometry.h"

class QuasiEuler{

	public:

		//public functions
		void S(StructuredGrid &mesh);//return S(x) at mesh points
			//TODO: write conversion functions and test them, could copy from hw1?
		void pressureSensor();


		//public parameters (the problem parameters)
		double gamma, inlet_pressure, total_temperature, s_star;
		double initial_pressure_left, initial_pressure_right, time;

};
