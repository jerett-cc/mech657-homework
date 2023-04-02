#ifndef GEOMETRY_H
#define GEOMETRY_H


#include<cassert>
#include<cmath>
#include<vector>
#include<string>
#include<fstream>
#include<stdio.h>
#include<regex>

#include<iostream>
#include<iomanip>

#include "../Eigen/Dense"
#include "test_with_exact.h"

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
    double get_q1(const int idx);
    double get_q2(const int idx);
    double get_q3(const int idx);

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

    double convert_pressure_to_energy(const double pressure, const double density, const double velocity, const double gamma){
//      std::cout
//                << pressure << std::endl
//                << density << std::endl
//                << velocity << std::endl
//                << gamma << std::endl;
      return pressure/(density*(gamma-1.0)) + 0.5 * velocity*velocity;
    }

    void updateBV();

//TODO: need to add an update boundary values function
  //private:
    double Left, Right;
    int highest_index, problem_size, R_index, L_index;
    double S(const double x);
    double Sprime(const double x);
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
  for(unsigned int i = 0; i<q.size(); ++i)
  {
    q[i] = q[i] + a_vec[i];
  }
}

/**
 * return Q as a single vector, rather than a std::vector of 3 component vectors
 * this would be the points on the interior
 */
Eigen::VectorXd ProblemData::getQVect(){
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(problem_size);
  int j = 0;
  for(unsigned int i = 0; i<q.size(); ++i)
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
    E(2) = Velocity(idx) * (Energy(idx) + Pressure(idx)) * S(X(idx));

    return E;
  }
}

/**
 * return double pressure at specified index
 */
double ProblemData::Pressure(const int idx){
  //return (parameter.gamma -1.0) * (Energy(idx) - 0.5/Density(idx) * (std::pow(Density(idx) * Velocity(idx), 2)));
return (parameter.gamma -1.0) * (Energy(idx) - 0.5/Density(idx) * Density(idx) * Velocity(idx) * Density(idx) * Velocity(idx));
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
 * TODO what to do if we are on problem 3? probably need to calculate
 */
void ProblemData::setInitialCondition(const double gamma,
                                        const double inlet_pressure,
                                        const double total_temperature,
                                        const double R,
                                      const double s_star){
  Eigen::Vector3d boundary;//FIXME do I need this vector anymore??
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
  //std::cout << "Mach left is  " << machl << std::endl;
  //std::cout << "Mach right is " << machr << std::endl;
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

  boundary(0) = densityl;
  boundary(1) = densityl * velocityl;

  std::cout << "pressure L    = " << pressurel << std::endl;
  std::cout << "density L     = " << densityl << std::endl;
  std::cout << "sound speed L = " << std::sqrt(parameter.gamma * pressurel / densityl) << std::endl;
  std::cout << "velocity L    = " << velocityl << std::endl;
  std::cout << "Energy L      = " << el/densityl << std::endl;
//  std::cout << "boundary Q(1) = " << boundary_Ql(1) << std::endl;
//  std::cout << "boundary E(0) = " << boundary_El(0) << std::endl;
//  std::cout << "boundary E(1) = " << boundary_El(1) << std::endl;
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

  //std::cout << "pressure R    = " << pressurer << std::endl;
  //std::cout << "density R     = " << densityr << std::endl;
  //std::cout << "sound speed R = " << std::sqrt(parameter.gamma * pressurer / densityr) << std::endl;
  //std::cout << "velocity R    = " << velocityr << std::endl;
  //std::cout << "Energy R      = " << er/densityl << std::endl;
  //std::cout << "boundary E(2) = " << boundary_Er(2) << std::endl;
  //std::cout << "____________---________________" << std::endl;

//  output the vectors on the boundary, should only differ by S(Right)/S(Left)
//  std::cout << "_______________________________" << std::endl;

  /*FIXME why are the values not zero?? Probably need to initialize with the primitive variables densityl, velocityl, and pressurer*/
  //initialize the entire domain to the same values as the boundary, with the proper

#if 1
  std::cout << "______________Initializing with constant value from BC____________" << std::endl;
  for(unsigned int i = 0; i<q.size(); i++)
  {
    q[i](0) = densityl*S(X(i));
    q[i](1) = densityl*velocityl*S(X(i));
    q[i](2) = densityl*convert_pressure_to_energy(pressurer, densityl, velocityl, parameter.gamma)*S(X(i));

//    q[i](0) = 1;
//    q[i](1) = 1;
//    q[i](2) = 1;
    std::cout << q[i] << std::endl;
  }

  boundary_Ql(0) = densityl*S(Left);//prescribed bc
  boundary_Ql(1) = densityl*velocityl*S(Left);
  boundary_Ql(2) = densityl*convert_pressure_to_energy(pressurel, densityl, velocityl, parameter.gamma)*S(Left);

  boundary_Qr(0) = densityr*S(Right);//prescribed bc
  boundary_Qr(1) = densityr*velocityr* S(Right);
  boundary_Qr(2) = densityr*convert_pressure_to_energy(pressurer, densityr, velocityr, parameter.gamma)*S(Right);

  boundary_El(0) = densityl * velocityl * S(Left);
  boundary_El(1) = (densityl * velocityl * velocityl + pressurel)* S(Left);
  double energyl = densityl * (pressurel/(densityl * (parameter.gamma - 1)) + std::pow(velocityl,2)/2);
  boundary_El(2) = velocityl * (energyl + pressurel) * S(Left);

  boundary_Er(0) = densityr * velocityr * S(Right);
  boundary_Er(1) = (densityr * velocityr * velocityr + pressurer)* S(Right);
  double energyr = densityr * (pressurer/(densityr * (parameter.gamma - 1)) + std::pow(velocityr,2)/2);
  boundary_Er(2) = velocityr * (energyr + pressurer) * S(Right);

#endif

#if 0
  std::cout << "______________Initializing with exact solution____________" << std::endl;
  Eigen::VectorXd density_vec(20+2);
  Eigen::VectorXd velocity_vec(20+2);
  Eigen::VectorXd pressure_vec(20+2);
  Eigen::VectorXd energy_vec(20+2);
  get_data("density_number.txt", density_vec);
  get_data("velocity_number.txt", velocity_vec);
  get_data("pressure_number.txt", pressure_vec);
  for(unsigned int i = 0; i<q.size(); i++)
  {
    q[i](0) = density_vec(i+1)*S(X(i));
    q[i](1) = density_vec(i+1)*velocity_vec(i+1)*S(X(i));
    q[i](2) = density_vec(i+1)*convert_pressure_to_energy(pressure_vec(i+1), density_vec(i+1), velocity_vec(i+1), parameter.gamma)*S(X(i));
    //std::cout << q[i] << std::endl;
  }
  boundary_Ql(0) = density_vec(0)*S(X(-1));//prescribed bc
  boundary_Ql(1) = density_vec(0)*velocity_vec(0)*S(X(-1));
  boundary_Ql(2) = density_vec(0)*convert_pressure_to_energy(pressure_vec(0), density_vec(0), velocity_vec(0), parameter.gamma)*S(X(-1));

  boundary_Qr(0) = density_vec(21)*S(X(22));//prescribed bc
  boundary_Qr(1) = density_vec(21)*velocity_vec(21)* S(X(22));
  boundary_Qr(2) = density_vec(21)*convert_pressure_to_energy(pressure_vec(21), density_vec(21), velocity_vec(21), parameter.gamma)*S(X(22));

  boundary_El(0) = density_vec(0) * velocity_vec(0) * S(X(-1));
  boundary_El(1) = (density_vec(0) * velocity_vec(0) * velocity_vec(0) + pressure_vec(0))* S(X(-1));
  double energyl = density_vec(0) * (pressure_vec(0)/(density_vec(0) * (parameter.gamma - 1)) + std::pow(velocity_vec(0),2)/2);
  boundary_El(2) = velocity_vec(0) * (energyl + pressure_vec(0)) * S(X(-1));

  boundary_Er(0) = density_vec(21) * velocity_vec(21) * S(X(22));
  boundary_Er(1) = (density_vec(21) * velocity_vec(21) * velocity_vec(21) + pressure_vec(21))* S(X(22));
  double energyr = density_vec(21) * (pressure_vec(21)/(density_vec(21) * (parameter.gamma - 1)) + std::pow(velocity_vec(21),2)/2);
  boundary_Er(2) = velocity_vec(21) * (energyr + pressure_vec(21)) * S(X(22));

#endif
  std::cout << "_____Initial vectors Ql Qr___________" << std::endl;
  std::cout << boundary_Ql << std::endl << " + " << std::endl <<  boundary_Qr << std::endl;
  std::cout << "_____Initial vectors El Er___________" << std::endl;
  std::cout << boundary_El << std::endl << " + " << std::endl <<  boundary_Er << std::endl;
}

void
ProblemData::updateBV()
{
  boundary_Ql(0) = boundary_density*S(Left);//prescribed bc
  boundary_Ql(1) = boundary_density*boundary_velocity*S(Left);
  boundary_Ql(2) = boundary_density*convert_pressure_to_energy(Pressure(0), boundary_density, boundary_velocity, parameter.gamma)*S(Left);

  boundary_Qr(0) = boundary_density*S(Right);//prescribed bc
  boundary_Qr(1) = boundary_density*boundary_velocity* S(Right);
  boundary_Qr(2) = boundary_density*convert_pressure_to_energy(Pressure(highest_index), boundary_density, boundary_velocity, parameter.gamma)*S(Right);

//  boundary_El(0) = densityl * velocityl * S(X(-1));
//  boundary_El(1) = (densityl * velocityl * velocityl + pressurel)* S(X(-1));
//  double energyl = densityl * (pressurel/(densityl * (parameter.gamma - 1)) + std::pow(velocityl,2)/2);
//  boundary_El(2) = velocityl * (energyl + pressurel) * S(X(-1));
//
//  boundary_Er(0) = densityr * velocityr * S(X(22));
//  boundary_Er(1) = (densityr * velocityr * velocityr + pressurer)* S(X(22));
//  double energyr = densityr * (pressurer/(densityr * (parameter.gamma - 1)) + std::pow(velocityr,2)/2);
//  boundary_Er(2) = velocityr * (energyr + pressurer) * S(X(22));
}

double ProblemData::S(const double x){
  if(x <= 5.0)
    {
      return 1.0 + 1.5*std::pow((1.0-x/5.0),2);
    }
    else if(x >5)
    {
      return 1.0 + 0.5*std::pow((1.0-x/5.0),2);
    }
    else return 0.0;
}

double ProblemData::Sprime(const double x){
  if(x < 5)
  {
    return -3./5.*(1.-x/5.);
  }
  else
  {
    return -1./5. * (1.-x/5.);
  }
}

double ProblemData::get_q1(const int idx){
  return Q(idx)(0)/S(X(idx));
}

double ProblemData::get_q2(const int idx){
  return Q(idx)(1)/S(X(idx));
}

double ProblemData::get_q3(const int idx){
  return Q(idx)(2)/S(X(idx));
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
  a_file.open(f_name+ "_velocity.out", std::ios::out | std::ios::trunc);

  for (int i = L_index; i< R_index+1; ++i)
  {
    a_file << std::setprecision(16) << X(i) << ", " << Velocity(i) << std::endl;
  }
  a_file.close();

  std::ofstream a_file1;
  a_file1.open(f_name+ "_energy.out", std::ios::out | std::ios::trunc);

  for (int i = L_index; i< R_index+1; ++i)
  {
    a_file1 <<  std::setprecision(16) << X(i) << ", " << Energy(i)/Density(i) << std::endl;
  }
  a_file1.close();

  std::ofstream a_file2;
  a_file2.open(f_name+ "_density.out", std::ios::out | std::ios::trunc);

  for (int i = L_index; i< R_index+1; ++i)
  {
    a_file2 << std::setprecision(16) << X(i) << ", " << Density(i) << std::endl;
  }
  a_file2.close();

  std::ofstream a_file3;
  a_file3.open(f_name+ "_pressure.out", std::ios::out | std::ios::trunc);

  for (int i = L_index; i< R_index+1; ++i)
  {
    a_file3 << std::setprecision(16) << X(i) << ", " << Pressure(i) <<std::endl;
  }
  a_file3.close();


  std::ofstream a_file4;
  a_file3.open(f_name+ "_convert.out", std::ios::out | std::ios::trunc);

  for (int i = 0; i< R_index; ++i)
  {
    a_file3 << std::setprecision(16) << X(i) << ", " << convert_pressure_to_energy(Pressure(i), Density(i), Velocity(i), parameter.gamma) <<std::endl;
  }
  a_file3.close();


}


























#endif
