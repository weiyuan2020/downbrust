/*
 * @Author:  L.F.Wang
 * @Date: 2021-08-01 21:12:19 
 * @Last Modified by: weiyuan
 * @Last Modified time: 2021-08-01 21:12:44
 */
#pragma once
// Path of "ThirdParties" should be included in the "include path" and "library path".
/////////////////////////////////////////////
// If FORTRAN is available, add "$FORTRAN_INSTALL_LOCATION$\compiler\lib\ia32" into the "library path", otherwise comment these.
// #include "ThirdParties/cblas/cblasInc.h"
// #include "ThirdParties/lapacke/lapackeInc.h"
// #define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE
/////////////////////////////////////////////
//#define EIGEN_NO_DEBUG
//#include <ThirdParties/Eigen/Core>
#include "Eigen/Dense"
//#include <ThirdParties/Eigen/LU>