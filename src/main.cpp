#include <iostream>

#include "downbrust.h"

using namespace DOWNBRUST;

int main() {
  downbrust test;
  //   double Array[10] = {0.0};
  //   int Idx = 0;
  //   Array[Idx++] = 0.0;     // dbPosP[0]
  //   Array[Idx++] = 0.0;     // dbPosP[1]
  //   Array[Idx++] = -610.0;  // dbPosP[2]
  //   Array[Idx++] = 0.0;     // dbAtt[0]
  //   Array[Idx++] = 0.0;     // dbAtt[1]
  //   Array[Idx++] = 0.0;     // dbAtt[2]
  //   Array[Idx++] = 458.0;   // dbCoreRadius 涡丝半径
  //   Array[Idx++] = 915.0;   // dbRingVortex 涡环半径
  //   Array[Idx++] = 10.0;    // WZREF
  //   test.Init(Array);
  test.Init("../input/ConfigData.txt");
  double x = 900;
  double z = -600;
  Vector3d pos;
  pos << x, 0.0, z;
  test.SetcgPos(pos);
  Vector3d Spd;
  test.GetcgSpd(&Spd);
  cout << "\nposition = " << pos.transpose();
  cout << "\nSpeed    = " << Spd.transpose();
  test.GetPlaneSpd("../input/SetPlaneSpeed.txt", "../output/PlaneSpeedoutput.txt");
  return 0;
}