#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "MathInc.h"

using Eigen::Matrix3d;
using Eigen::Vector3d;
using MATH::SQ;
using std::cout;
using std::ifstream;
using std::string;

namespace DOWNBRUST {
class downbrust {
 private:
  /* init data */
  Vector3d dbPosP;      // prime downbrust position (m) forwad-right-down axis
  Vector3d dbAtt;       // downbrust attitude angle (deg)
  double dbCoreRadius;  // radius of core of rotating vortex material 涡丝半径
  double dbRingVortex;  // reference radius of ring vortex 涡环半径
  double WZREF;         // Wz reference speed (m/s)

  /* cal input data */
  Vector3d arcftCGPos;  // aircraft C.G. position (m)

  /* middle data */
  Vector3d dbPosM;      // mirror downbrust position
  Vector3d cg2dbPos;    // arcftCGPos - dbPos
  Matrix3d matp, matm;  // matrix of prime and mirror ring vortex

  double circrv;  // circulation strength of primary ring vortex,

  /* output data */
  Vector3d arcftCGSpd;  // aircraft speed

  int calculateSpd();
  //   int STRMFN(double hinput, double rinput, double* strfoutput);
  //   int calSpeed(Vector3d relativePosinput, double* woutput);
  int STRMFN_partialPart(Vector3d relePosinput, double* phi_r_zoutput);
  int calPointSpd(Vector3d posP, Vector3d posM, Vector3d* actualSpd);

 public:
  downbrust(/* args */);
  ~downbrust();
  int Init(double* initArray);
  int Init(string inputfilename);
  int SetcgPos(const Vector3d cgposInput);
  int GetcgSpd(Vector3d* cgspeedoutput);
  int GetPlaneSpd(string inputfilename, string outputfilename);
};
}  // namespace DOWNBRUST