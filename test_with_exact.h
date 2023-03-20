#ifndef _TEST_WITH_EXACT_H_
#define _TEST_WITH_EXACT_H_


#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <cmath>
#include <cassert>
#include <vector>
#include "../Eigen/Dense"

// bool compare_double(double x, double y, float epsilon = 1e-14)
// {
//   if(fabs(x - y) < epsilon)
//     return true; //they are same
//   return false; //they are not same
// }
// bool compare_vector(Eigen::VectorXd x, Eigen::VectorXd y, float epsilon = 1e-14)
// {
//   assert(x.size() == y.size());
//   for(int i = 0; i < x.size(); i++)
//     {
//       if(compare_double(x(i), y(i)), epsilon)
//         {
//         }
//       else
//         {
//           return false;
//         }
//     }
//       return true;
// }


void
get_data(const std::string a_file_name, Eigen::VectorXd& a_vector)
{

  std::ifstream file;
  file.open(a_file_name);
  std::string operating_line;
  double val;
  int i = 0;
  while(file)
    {
      std::regex pattern("(.*)(  )(.*)");
      std::getline(file, operating_line);
      std::string operating_value = std::regex_replace(operating_line,pattern, "$1");
      pattern =(" ");
      operating_value = std::regex_replace(operating_value,pattern, "");
      if(! operating_value.empty())
        {
          val = std::stod(operating_value);
          a_vector(i) = val;
          i++;
        }
    }


}





#endif  // _TEST_WITH_EXACT_H_
