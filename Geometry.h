#ifndef GEOMETRY_H
#define GEOMETRY_H


#include<cassert>
#include<cmath>
#include<vector>

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

class ProblemData{
  public:
    std::vector<Eigen::Vector3d> q;
    Eigen::VectorXd x;
    Eigen::Vector3d boundary_El;//todo
    Eigen::Vector3d boundary_Ql;
    Eigen::Vector3d boundary_Er;//todo
    Eigen::Vector3d boundary_Qr;//todo need this?

    double dx;

    Eigen::Vector3d E(const int idx);
    double Pressure(const int idx);
    double Density(const int idx);//todo make density at each node point in interior
    double Temperature(const int idx);//todo make temp at each node point in interior
    double Mach(const int idx);//todo make mach at each node point in interior


    const Eigen::Vector3d& operator[](const int idx) const;
    Eigen::Vector3d& operator[](const int idx);
    void setInitialCondition(const double gamma,
                             const double inlet_pressure,
                             const double total_temperature,
                             const double R,
                             const double s_star);
    Eigen::VectorXd getQVect();

    ProblemData(int number_nodes, int number_ghost, double interval_l, double interval_r, parameters param)
    :Left(interval_l), Right(interval_r), highest_index(number_nodes-1),
     x(number_nodes), q(number_nodes), problem_size(number_nodes*3),
     parameter(param)
    {
//      std::cout << number_nodes << std::endl;
      dx = (Right-Left)/(number_nodes+1);
      std::cout << "dx is " <<  dx << std::endl;
      for (int i=0; i < x.size(); ++i)
      {
        x(i) = Left + (i+1)*dx;
//        std::cout << X(i) << std::endl;
      }
//      std::cout << highest_index << std::endl;

    }


  private:
    double Left, Right;
    int highest_index, problem_size;
    const parameters parameter;
    double soundSpeed(const int position);
    double S(const double x);
    Eigen::Vector3d Q(const int idx);//todo
    double X(const int idx);//todo
    double nonlinearFunctionToSolveP1(double M,double s_star, double gamma);
    double nonlinearFunctionToSolveP1Deriv(double M,double s_star, double gamma);
    bool outOfBounds(int idx);
};

/**
 * returns 1 if index is out of bounds for Q, and 0 if not.
 */
bool ProblemData::outOfBounds(int idx){
  bool out_of_bounds = (idx <0 || idx > highest_index) ? 1:0;
  return out_of_bounds;
}

/**
 * Returns Q at a specified index, supporting indexing out of bounds by simply using
 * a constant extrapolation.
 */
const Eigen::Vector3d& ProblemData::operator[](const int idx) const{
  bool out_of_bounds = (idx <0 || idx > highest_index) ? 1:0;
  //if we are accessing outside of our bounds, do this logic
  if (out_of_bounds)
  {
    // if we are trying to access a value left of our leftmost, return out left boundary value
    if (idx<0)
    {
      return boundary_Ql;
    }
    else // if we are trying to access a value right of out rightmost value, just
         // do a constant extrapolation of our R value.
    {
      return q[highest_index];
    }
  }
  else // if not out of bounds, return Q at this point
  {
    return q[idx];
  }
}
/**
 * Return a modifiable Q at a specified index, supporting indexing out of bounds by
 * a constant extrapolation.
 */
Eigen::Vector3d& ProblemData::operator[](const int idx){
  bool out_of_bounds = (idx <0 || idx > highest_index) ? 1:0;
  //if we are accessing outside of our bounds, do this logic
  if (outOfBounds(idx))
  {
    // if we are trying to access a value left of our leftmost, return out left boundary value
    if (idx<0)
    {
      return boundary_Ql;
    }
    else // if we are trying to access a value right of out rightmost value, just
         // do a constant extrapolation of our R value.
    {
      return q[highest_index];
    }
  }
  else // if not out of bounds, return Q at this point
  {
    return q[idx];
  }
}

/**
 * return Q as a single vector, rather than a std::vector of 3 component vectors
 * this would be the points on the interior
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
 * vector that is the left most or right most allowed value of E
 */
Eigen::Vector3d ProblemData::E(const int idx){
  Eigen::Vector3d E = Eigen::Vector3d::Zero();
    double q1 = q[idx](0);
    double q2 = q[idx](1);
    double q3 = q[idx](2);
    E(0) = q2;
    E(1) = q2 + Pressure(idx);
    E(2) = q3 + Density(idx)*Pressure(idx) - Pressure(idx);

    return E;
}
/**
 * return double pressure at specified index
 */
double ProblemData::Pressure(const int idx){
  std::cout << "++++" << std::endl
            << "x[-1] is " << std::endl
            << X(-1) << std::endl
            <<  "S(X[-1]) is " << std::endl
            << S(X(idx)) << std::endl;
  double q1 = Q(idx)(0)/(S(X(idx)));
  double q2 = Q(idx)(1)/(S(X(idx)));
  double q3 = Q(idx)(2)/(S(X(idx)));
  return (parameter.gamma-1) * (q3 - 1/ (2*q1)*std::pow(q2,2));
}

/**
 * return double mach at specified index
 */
double ProblemData::Mach(const int idx){
//  std::cout << "++++" << std::endl
//            << "q[-1] is " << std::endl
//            << Q(-1) << std::endl;
  double q1 = Q(idx)(0);
  double q2 = Q(idx)(1);
  double q3 = Q(idx)(2);
  return (q2/q1)/(soundSpeed(idx));
}

/**
 * return double density at specified index
 */
double ProblemData::Density(const int idx){
  double q1 = Q(idx)(0);
  double q2 = Q(idx)(1);
  double q3 = Q(idx)(2);
  return q1;
}

/**
 * return double density at specified index
 */
double ProblemData::Temperature(const int idx){
  assert(0 && " calculating temp not implemented yet");
  double q1 = Q(idx)(0);
  double q2 = Q(idx)(1);
  double q3 = Q(idx)(2);
  return q1;//todo calculate this.
}

/**
 * return Q at specified index, this is for internal class computations only
 */
Eigen::Vector3d ProblemData::Q(const int idx){
  return operator [](idx);
}

/**
 * return X at specified index, this is for internal class computations only
 */
double ProblemData::X(const int idx){
//  std::cout << outOfBounds(idx) <<std::endl;
  if (outOfBounds(idx))
  {
    return dx*(idx+1) + Left;
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
  double mach = 0.1;
  while(std::fabs(nonlinearFunctionToSolveP1(mach, s_star, gamma))>1e-13)
  {
    mach = mach -nonlinearFunctionToSolveP1(mach, s_star, gamma)/
        nonlinearFunctionToSolveP1Deriv(mach, s_star, gamma);
  }
  double inside = 1 + (gamma-1)/(2)*std::pow(mach,2);
  double pressure = inlet_pressure * std::pow(inside,-gamma/(gamma-1));
  double temperature = total_temperature /inside;
  double dens = pressure / (R* temperature);
  double velocity = std::sqrt(gamma*pressure/dens) * mach ;
  double a = std::pow(gamma*pressure/dens, 0.5);
  double a2 = std::sqrt(gamma*R* temperature);
  double momentum = velocity * dens;
  double epsilon = temperature * R / (gamma-1);
  double energy = (temperature * R / (gamma-1) + std::pow(velocity, 2)/2);
  double energy_momentum = dens * energy;
  std::cout //<< "density = " << dens << std::endl
//            << "sound speed is "<< a << std::endl
//            << "sound speed2 is "<< a2 << std::endl
            << "pressure is "<< pressure << std::endl
//            << "temperature is "<< temperature << std::endl
//            << "velocity is = " << velocity << std::endl
//            << "momentum = " << momentum << std::endl
//            << "energy moment= " << energy_momentum/dens << std::endl
//            << "mach = " << mach << std::endl
//            << "R = " << R << std::endl
//            << "epsilon = " << epsilon << std::endl
//            << "epsilon +vel^2= " << epsilon + 0.5 * velocity* velocity
              << std::endl;
  boundary_Ql(0) = dens*S(Left);
  boundary_Ql(1) = momentum*S(Left);
  boundary_Ql(2) = energy_momentum*S(Left);

  for (int i = 0; i<q.size(); ++i)
  {
    q[i](0) = dens * S(x(i));
    q[i](1) = momentum * S(x(i));
    q[i](2) = energy_momentum * S(x(i));
    std::cout << "S at index " << i << " is " << S(x(i)) << std::endl;
  }
  std::cout << "Left boundary condition is " <<std::endl
            <<  boundary_Ql << std::endl;

  std::cout<< "Initial Q (interior) is " << std::endl;
  for (int i = 0; i< q.size(); ++i)
  {
    std::cout << q[i] << std::endl;
  }
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

double ProblemData::nonlinearFunctionToSolveP1(double M, double s_star, double gamma){

  double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
  double exponent = (gamma + 1)/(2*(gamma-1));
  assert(fabs(exponent - 3)<1e-13);
  return pow(inside, exponent) - (S(Left)*M)/(s_star);
}

double ProblemData::nonlinearFunctionToSolveP1Deriv(double M, double s_star, double gamma){
  double inside = 2/(gamma +1) + (gamma - 1)/(gamma+1)*pow(M,2);
  double exponent = (gamma + 1)/(2*(gamma-1));
  return M*pow(inside, exponent-1) - (S(Left))/(s_star);
}




























#endif