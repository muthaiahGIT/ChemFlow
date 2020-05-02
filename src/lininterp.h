#ifndef CTF_LININTERP_H_
#define CTF_LININTERP_H_

#include <vector>
#include <Eigen/Dense>


double lininterp(const double x, const std::vector<double>& xOld,
               const std::vector<double>& yOld);

double lininterp(const double x, const Eigen::VectorXd& xOld,
                 const Eigen::VectorXd& yOld);

// #ifdef NoRepository
//     #include "lininterp.cpp"
// #endif

#endif  // CTF_LININTERP_H_
