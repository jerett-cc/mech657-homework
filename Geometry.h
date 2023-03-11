#ifndef GEOMETRY_H
#define GEOMETRY_H


#include<cassert>
#include<cmath>
#include<vector>
#include<string>
#include<fstream>
#include<stdio.h>

#include<iostream>

#include "../Eigen/Dense"
//TODO: include a way to describe the face of a cell once we start doing cell-cenered Finite Volume code

struct parameters{
    const double inlet_pressure, total_temp, gamma, R, s_star;

    parameters(double Inlet_Pressure, double Total_Temp, double Gamma, double gas_constant_R,
               double S_Star):
                 inlet_pressure(Inlet_Pressure), total_temp(Total_Temp), gamma(Gamma),
                 R(gas_constant_R), s_star(S_Star){}
    //copy constructor is default
    parameters(const parameters&) = default;
};

/**
 * this class stores the solution vector and boundary conditions for the homework,
 * as well as methods to calculate the basic variables like pressure, density, mach, temp.
 */
class ProblemData{
  public:
    std::vector<Eigen::Vector3d> q;
    Eigen::VectorXd x;
    Eigen::Vector3d boundary_El;
    Eigen::Vector3d boundary_Ql;
    Eigen::Vector3d boundary_Er;//todo need this?
    Eigen::Vector3d boundary_Qr;//todo need this?

    double dx, boundary_velocity, boundary_density, boundary_pressure;
    const parameters parameter;

    Eigen::Vector3d E(const int idx);
    double Pressure(const int idx);
    double Density(const int idx);
    double Temperature(const int idx);//todo make temp at each node point in interior
    double Mach(const int idx);
    double Velocity(const int idx);
    double initial_pressure;


    const Eigen::Vector3d& operator[](const int idx) const;
    Eigen::Vector3d operator[](const int idx);
    void operator+=(const std::vector<Eigen::Vector3d> &a_vec);
    void setInitialCondition(const double gamma,
                             const double inlet_pressure,
                             const double total_temperature,
                             const double R,
                             const double s_star);
    Eigen::VectorXd getQVect();
    double soundSpeed(const int position);
    double Energy(const int idx);
    Eigen::Vector3d Q(const int idx);
    void printQuantities(std::string f_name);

    ProblemData(int number_nodes, int number_ghost, double interval_l, double interval_r, parameters param)
    :Left(interval_l), Right(interval_r), highest_index(number_nodes-1),
     x(number_nodes), q(number_nodes), problem_size(number_nodes*3),
     parameter(param), R_index(number_nodes), L_index(-1)
    {
//      std::cout << number_nodes << std::endl;
      dx = (Right-Left)/(number_nodes+1);
//      std::cout << "dx is " <<  dx << std::endl;
      for (int i=0; i < x.size(); ++i)
      {
        x(i) = Left + (i+1)*dx;
//        std::cout << X(i) << std::endl;
      }
//      std::cout << highest_index << std::endl;

    }

//TODO: need to add an update boundary values function
  //private:
    double Left, Right;
    int highest_index, problem_size, R_index, L_index;
    double S(const double x);
    double X(const int idx);
    double nonlinearFunctionToSolveP1(double M,double s_star, double gamma, double x);
    double nonlinearFunctionToSolveP1Deriv(double M,double s_star, double gamma, double x);
    bool outOfBounds(int idx);
};

/**
 * this function returns the energy of the system
 *
 *              e = epsilon + v^2/2
 *
 * which can be determined from Q by dividing the first and last entries in vector Q.
 */
double ProblemData::Energy(const int idx){
      return Q(idx)(2)/S(X(idx));
}

/**
 * returns 1 if index is out of bounds for Q, and 0 if not.
 */
bool ProblemData::outOfBounds(int idx){
  bool out_of_bounds = (idx <0 || idx > highest_index) ? 1:0;
  return out_of_bounds;
}

/**
 * Returns constant access to Q at a specified index, supporting indexing out of bounds by simply using
 * a constant extrapolation.
 */
const Eigen::Vector3d& ProblemData::operator[](const int idx) const{
  bool out_of_bounds = (idx <0 || idx > highest_index) ? 1:0;

  if (out_of_bounds)
  {
    if (idx<0)
    {
      return boundary_Ql;
    }
    else
    {
      return boundary_Qr;
    }
  }
  else
  {
    return q[idx];
  }
}

/**
 * Return a modifiable Q at a specified index, supporting indexing out of bounds by
 * a constant extrapolation.
 */
Eigen::Vector3d ProblemData::operator[](const int idx){

  if (outOfBounds(idx))
  {

    if (idx<0)
    {
      return boundary_Ql;
    }
    else
    {
      return boundary_Qr;
    }
  }
  else
  {
    return q[idx];
  }
}

/**
 * Add a vector to Q (alias data), call would be data += step_solution
 * todo test this.
 */
void ProblemData::operator+=(const std::vector<Eigen::Vector3d> &a_vec){
  std::cout << "made it into +=\n";
  assert(a_vec.size()==q.size() && "trying to add a vector of different size to q");
  for(int i = 0; i<q.size(); ++i)
  {
    q[i] = q[i] + a_vec[i];
  }
}

/**
 * return Q as a single vector, rather than a std::vector of 3 component vectors
 * this would be the points on the interior, FIXME test this
 */
Eigen::VectorXd ProblemData::getQVect(){
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(problem_size);
  int j = 0;
  for(int i = 0; i<q.size(); ++i)
  {
    tmp(j) = q[i](0);
    tmp(j+1) = q[i](1);
    tmp(j+2) = q[i](2);
    j+=3;
  }
  return tmp;
}
/**
 * this should return E at the specified index, if idx is out of bounds, it returns a constant
 * vector that is the left most or right most allowed value of E FIXME TEST THIS
 */
Eigen::Vector3d ProblemData::E(const int idx){//todo need to verify this works
  Eigen::Vector3d E = Eigen::Vector3d::Zero();

  if (idx < 0)
  {
    return boundary_El;
  }
  else if (idx > highest_index)
  {
    return boundary_Er;
  }
  else
  {
    E(0) = Density(idx) * Velocity(idx) * S(X(idx));
    E(1) = (Density(idx) * Velocity(idx) * Velocity(idx) +  Pressure(idx)) * S(X(idx));
    E(2) = Velocity(idx) * (Density(idx) * Energy(idx) + Pressure(idx)) * S(X(idx));

    return E;
  }
}

/**
 * return double pressure at specified index
 */
double ProblemData::Pressure(const int idx){
  return (parameter.gamma -1) * (Energy(idx) - 1/(2*Density(idx)) * (std::pow(Density(idx) * Velocity(idx), 2)));
}

/**
 * return double velocity at specified index
 */
double ProblemData::Velocity(const int idx){
  double q1 = Q(idx)(0);
  double q2 = Q(idx)(1);
  return q2/q1;
}

/**
 * return double mach at specified index
 */
double ProblemData::Mach(const int idx){
  double v = Velocity(idx);
  return v/soundSpeed(idx);
}

/**
 * return double density at specified index
 */
double ProblemData::Density(const int idx){
  return Q(idx)(0)/S(X(idx));
}

/**
 * return double temperature at specified index
 */
double ProblemData::Temperature(const int idx){
  return Pressure(idx)/(parameter.R * Density(idx));
}

/**
 * return Q at specified index
 */
Eigen::Vector3d ProblemData::Q(const int idx){
  return operator[](idx);
}

/**
 * return X at specified index
 */
double ProblemData::X(const int idx){
//  std::cout << outOfBounds(idx) <<std::endl;
  if (outOfBounds(idx))
  {
    if (idx<0)
    {
      return Left;
    }
    else
    {
      return Right;
    }
  }
  else
  {
    return x(idx);
  }
}

/**
 * returns sound speed in specified index
 *  needs pressure and density.
 */
double ProblemData::soundSpeed(const int idx){
  return std::sqrt(parameter.gamma * Pressure(idx)/ Density(idx));
}

/**
 * set initial conditions based on typical problem parameters
 * TODO what to do is we are on problem 3? probably need to calculate
 */
void ProblemData::setInitialCondition(const double gamma,
                                        const double inlet_pressure,
                                        const double total_temperature,
                                        const double R,
                                      const double s_star){
  //calculate mach numbers on left and right of intervals.
  double machl = 0.1;
  double machr = machl;
  while(std::fabs(nonlinearFunctionToSolveP1(machl, s_star, gamma, Left))>1e-13)
  {
    machl = machl -nonlinearFunctionToSolveP1(machl, s_star, gamma, Left)/
       nonlinearFunctionToSolveP1Deriv(machl, s_star, gamma, Left);
    machr = machr -nonlinearFunctionToSolveP1(machr, s_star, gamma, Right)/
       nonlinearFunctionToSolveP1Deriv(machr, s_star, gamma, Right);
  }


  std::cout << "Mach left is  " << machl << std::endl;
  std::cout << "Mach right is " << machr << std::endl;
  std::cout << "____________Left________________" << std::endl;
  //calculate pressure and velocity on left of interval, and put into boundary Q_l
  double insidel = 1 + (gamma-1)/(2)*std::pow(machl,2);
  double pressurel = inlet_pressure * std::pow(insidel, -gamma/(gamma-1));
  double temperaturel = total_temperature / insidel;
  double densityl = pressurel / (R* temperaturel);
  double velocityl = std::sqrt(gamma*pressurel/densityl) * machl;
  double epsilonl = parameter.R/(parameter.gamma - 1) * pressurel / (densityl * parameter.R);
  double el = densityl * (epsilonl + 0.5 * std::pow(velocityl,2));

  boundary_velocity = velocityl;
  boundary_density = densityl;

  boundary_Ql(0) = densityl * S(Left);
  boundary_Ql(1) = densityl * velocityl * S(Left);

  boundary_El(0) = boundary_Ql(1);
  boundary_El(1) = std::pow(boundary_Ql(1),2)/boundary_Ql(0) + pressurel*S(Left);

  std::cout << "pressure L    = " << pressurel << std::endl;
  std::cout << "density L     = " << densityl << std::endl;
  std::cout << "sound speed L = " << std::sqrt(parameter.gamma * pressurel / densityl) << std::endl;
  std::cout << "velocity L    = " << velocityl << std::endl;
  std::cout << "Energy L      = " << el << std::endl;
  std::cout << "boundary Q(1) = " << boundary_Ql(1) << std::endl;
  std::cout << "boundary E(0) = " << boundary_El(0) << std::endl;
  std::cout << "boundary E(1) = " << boundary_El(1) << std::endl;


  //calculate pressure and velocity on right of the interval, and put into boundary_Qr
  //this one we can only prescribe the pressure, the other variables should be interpolated from the interior
  std::cout << "____________Right________________" << std::endl;
  double insider = 1 + (gamma-1)/(2)*std::pow(machr,2);

  double pressurer = inlet_pressure * std::pow(insider,-gamma/(gamma-1));
  double temperaturer = total_temperature /insider;
  double densityr = pressurer / (R* temperaturer);
  double velocityr = std::sqrt(gamma*pressurer/densityr) * machr;
  double epsilonr = parameter.R/(parameter.gamma - 1) * pressurer / (densityl * parameter.R);
  double er = densityl * (epsilonr + 0.5 * std::pow(velocityl,2));

  boundary_pressure = pressurer;

  boundary_Qr(2) = densityr * (pressurer + velocityr*velocityr /2) * S(Right);

  boundary_Er(2) = velocityr * (densityr * (pressurer + velocityr*velocityr /2) + pressurer) * S(Right);

  std::cout << "pressure R    = " << pressurer << std::endl;
  std::cout << "density R     = " << densityr << std::endl;
  std::cout << "sound speed R = " << std::sqrt(parameter.gamma * pressurer / densityr) << std::endl;
  std::cout << "velocity R    = " << velocityr << std::endl;
  std::cout << "Energy L      = " << er/densityl << std::endl;
  std::cout << "boundary E(2) = " << boundary_Er(2) << std::endl;
  std::cout << "____________---________________" << std::endl;

  //initialize the other components of the boundary vectors with the proper values.

  boundary_Ql(2) = densityl * (pressurer + velocityl*velocityl/2) *S(Left);
  boundary_El(2) = velocityl * (densityl * (pressurer + velocityl*velocityl/2) + pressurer) * S(Left);

  boundary_Qr(0) = densityl * S(Right);
  boundary_Qr(1) = densityl * velocityl * S(Right);
  boundary_Er(0) = boundary_Qr(1);
  boundary_Er(1) = densityl * velocityl * velocityl + pressurel;

  std::cout << "_____Initial vectors___________" << std::endl;
  std::cout << boundary_Ql << std::endl << boundary_Qr << std::endl;
  std::cout << boundary_El << std::endl << boundary_Er << std::endl;
  std::cout << "_______________________________" << std::endl;

  /*FIXME why are the values not zero?? Probably need to initialize with the primitive variables densityl, velocityl, and pressurer*/
  //initialize the entire domain to the same values as the boundary, with the proper
  for (int i = 0; i < q.size(); ++i)
    {
	  q[i](0) = densityl * S(X(i));
	  q[i](1) = densityl * velocityl * S(X(i));
	  q[i](2) = er * S(X(i));
//      q[i](0) = boundary_Ql(0)*S(x(i))/S(Left) ;
//      q[i](1) = boundary_Ql(1)*S(x(i))/S(Left) ;
//      q[i](2) = boundary_Ql(2)*S(x(i))/S(Left) ;
      std::cout << q[i] << std::endl;
    }


  //double a = std::pow(gamma*pressure/dens, 0.5);
  //double a2 = std::sqrt(gamma*R* temperature);
  //double momentum = velocity * dens;
  //double epsilon = temperature * R / (gamma-1);
  //double energy = (temperature * R / (gamma-1) + std::pow(velocity, 2)/2);
  //double energy_momentum = dens * energy;
//  std::cout << "density = " << dens << std::endl
//            << "sound speed is "<< a << std::endl
//            << "sound speed2 is "<< a2 << std::endl
//            << "pressure is "<< pressure << std::endl
//            << "temperature is "<< temperature << std::endl
//            << "velocity is = " << velocity << std::endl
//            << "momentum = " << momentum << std::endl
//            << "energy moment= " << energy_momentum/dens << std::endl
//            << "mach = " << mach << std::endl
//            << "R = " << R << std::endl
//            << "epsilon = " << epsilon << std::endl
//            << "epsilon +vel^2= " << epsilon + 0.5 * velocity* velocity
//              << "S is " << S(Left) << std::endl
//              << "energy is " << energy_momentum << std::endl
//              << std::endl;
  //set the left boundary.
  //boundary_Ql(0) = dens*S(Left);
  //boundary_Ql(1) = momentum*S(Left);
  //boundary_Ql(2) = energy_momentum*S(Left);

  //boundary_El(0) = boundary_Ql(1);
  //boundary_El(1) = std::pow(boundary_Ql(1),2)/boundary_Ql(0) + pressure*S(Left);
  //boundary_El(2) = boundary_Ql(1)*boundary_Ql(2)/boundary_Ql(0)
    //                 + pressure*S(Left)* boundary_Ql(1)/boundary_Ql(0);


  //for (int i = 0; i < q.size(); ++i)
  //{
    //q[i](0) = dens * S(x(i));
    //q[i](1) = momentum * S(x(i));
    //q[i](2) = energy_momentum * S(x(i));
//    std::cout << "S at index " << i << " is " << S(x(i)) << std::endl;
  //}
  //std::cout << "Left boundary condition is " <<std::endl
    //        <<  boundary_Ql << std::endl
      //      << "Left boundary condition on E is " <<std::endl
        //    <<  boundary_El << std::endl;

  //std::cout<< "Initial Q (interior) is " << std::endl;
  //for (int i = 0; i< q.size(); ++i)
  //{
    //std::cout << q[i] << std::endl;
  //}
  //std::cout<< "Initial E (interior) is " << std::endl << "__________________________________________________" << std::endl;
  //for (int i = 0; i< q.size(); ++i)
    //{
      //std::cout << E(i) << std::endl;
    //}

}

double ProblemData::S(const double x){
  if(x <= 5)
    {
      return 1.0 + 1.5*std::pow(1-x/5,2);
    }
    else if(x >5)
    {
      return 1.0 + 0.5*std::pow(1-x/5,2);
    }
    else return 0.0;
}

double ProblemData::nonlinearFunctionToSolveP1(double M, double s_star, double gamma, double x){

  double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
  double exponent = (gamma + 1)/(2*(gamma-1));
  assert(fabs(exponent - 3)<1e-13);
  return pow(inside, exponent) - (S(x)*M)/(s_star);
}

double ProblemData::nonlinearFunctionToSolveP1Deriv(double M, double s_star, double gamma, double x){
  double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
  double exponent = (gamma + 1)/(2*(gamma-1));
  return M*pow(inside, exponent-1) - (S(x))/(s_star);
}
/**
 * calculates and prints quantities to csv. To be used after the calculation is done.
 */
void ProblemData::printQuantities(std::string f_name){

  std::ofstream a_file;
  a_file.open(f_name+ "_velocity.csv", std::ios::out | std::ios::trunc);

  for (int i = L_index; i< R_index+1; ++i)
  {
    a_file <<  " " << Left + (i+1)*dx  << ", " << Velocity(i) << std::endl;
  }
  a_file.close();

//  std::ofstream a_file1;
//  a_file1.open(f_name+ "_temp.csv", std::ios::out | std::ios::trunc);
//
//  for (int i = L_index; i< R_index+1; ++i)
//  {
//    a_file1 <<  " " << X(i)  << ", " << Temperature(i) << std::endl;
//  }
//  a_file1.close();

  std::ofstream a_file2;
  a_file2.open(f_name+ "_density.csv", std::ios::out | std::ios::trunc);

  for (int i = L_index; i< R_index+1; ++i)
  {
    a_file2 <<  " " << Left + (i+1)*dx  << ", " << Density(i) << std::endl;
  }
  a_file2.close();

  std::ofstream a_file3;
  a_file3.open(f_name+ "_pressure.csv", std::ios::out | std::ios::trunc);

  for (int i = L_index; i< R_index+1; ++i)
  {
    a_file3 <<  " " << Left + (i+1)*dx  << ", " << Pressure(i) <<std::endl;
  }
  a_file3.close();

}


























#endif
