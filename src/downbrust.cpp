#include "downbrust.h"

namespace DOWNBRUST {
downbrust::downbrust(/* args */) {}
downbrust::~downbrust() {}
int downbrust::Init(double* initArray) {
  int Idx = 0;
  dbPosP[0] = initArray[Idx++];
  dbPosP[1] = initArray[Idx++];
  dbPosP[2] = initArray[Idx++];
  dbAtt[0] = initArray[Idx++];
  dbAtt[1] = initArray[Idx++];
  dbAtt[2] = initArray[Idx++];
  dbCoreRadius = initArray[Idx++];
  dbRingVortex = initArray[Idx++];
  WZREF = initArray[Idx++];

  // 下涡环中心位置
  dbPosM = dbPosP;
  dbPosM[2] *= -1;
  return 0;
}

int downbrust::Init(string inputfilename) {
  ifstream infile;
  infile.open(inputfilename);

  int Idx = 0;
  double InputArray[9];
  while (!infile.eof() | Idx < 9) {
    infile >> InputArray[Idx++];
  }
  Init(InputArray);

  infile.close();
  return 0;
}

int downbrust::SetcgPos(const Vector3d cgposInput) {
  arcftCGPos = cgposInput;
  return 0;
}

int downbrust::GetcgSpd(Vector3d* cgspeedoutput) {
  calculateSpd();
  *cgspeedoutput = arcftCGSpd;
  return 0;
}

int downbrust::GetPlaneSpd(string inputfilename, string outputfilename) {
  ifstream infile;
  infile.open(inputfilename);

  int ArrIdx = 0;
  double InputArray[9];
  while (!infile.eof() | ArrIdx < 9) {
    infile >> InputArray[ArrIdx++];
  }
  infile.close();

  ArrIdx = 0;
  double xbgn = InputArray[ArrIdx++];
  double xstp = InputArray[ArrIdx++];
  double xend = InputArray[ArrIdx++];
  double zbgn = InputArray[ArrIdx++];
  double zstp = InputArray[ArrIdx++];
  double zend = InputArray[ArrIdx++];
  double ypos = InputArray[ArrIdx++];

  int xDim = ceil(abs(xend - xbgn) / xstp);
  int zDim = ceil(abs(zend - zbgn) / zstp);
  double* xArr = new double[xDim];
  double* zArr = new double[zDim];

  cout<<'\n';
  for (int xIdx = 0; xIdx < xDim; xIdx++) {
    xArr[xIdx] = xbgn + xstp * xIdx;
    cout<<xArr[xIdx]<<'\t';
  }
  for (int zIdx = 0; zIdx < zDim; zIdx++) {
    zArr[zIdx] = zbgn + zstp * zIdx;
    cout<<zArr[zIdx]<<'\t';
  }

  std::ofstream out;
  out.open(outputfilename);

  Vector3d positionVec, veloVec;

  for (int xIdx = 0; xIdx < xDim; xIdx++) {
    for (int zIdx = 0; zIdx < zDim; zIdx++) {
      positionVec << xArr[xIdx], ypos, -zArr[zIdx];
      // positionVec << ypos, xArr[xIdx], -zArr[zIdx];
      SetcgPos(positionVec);
      GetcgSpd(&veloVec);
      out << std::setprecision(8);
      out << positionVec[0] << '\t' << positionVec[1] << '\t' << positionVec[2];
      out << '\t';
      out << veloVec[0] << '\t' << veloVec[1] << '\t' << veloVec[2];
      out << '\n';
    }
  }
  out.close();

  delete xArr;
  delete zArr;
  return 0;
}

int downbrust::calculateSpd() {
  // data init
  double RDVORT = dbRingVortex;
  double HGVORT = dbPosP[2];
  // vortex circulation strength
  circrv = WZREF * 2.0 * RDVORT /
           (1.0 - 1.0 / pow(1.0 + SQ(2.0 * HGVORT / RDVORT), 1.5));

  // wz of the vertical axis of vortex ring
  double HGCG = arcftCGPos[2];
  double kwzp = pow(1.0 + SQ((HGVORT - HGCG) / RDVORT), 1.5);
  double kwzm = pow(1.0 + SQ((-HGVORT - HGCG) / RDVORT), 1.5);
  double WZ = circrv / 2.0 / RDVORT * (1.0 / kwzp - 1.0 / kwzm);

  // REGION 1: on the vertical axis of vortex
  // relative x,y and r = sqrt(x^2+y^2);
  MATH::ElrangToDcm21(213)(dbAtt, &matp);  // prime vortex: theta -> phi
  //   MATH::ElrangToDcm21(213)(-dbAtt, &matm);  // mirror vortex: -theta->-phi
  matm = matp.transpose();

  cg2dbPos = matp * (arcftCGPos - dbPosP);
  double XCGDBC = cg2dbPos[0];
  double YCGDBC = cg2dbPos[1];
  double ZCGDBC = cg2dbPos[2];
  double RVAXCG = sqrt(SQ(XCGDBC) + SQ(YCGDBC));

  double WX, WY;
  double R1PSML = sqrt(SQ(ZCGDBC) + SQ(RVAXCG - RDVORT));
  if (RVAXCG <= 1e-3) {
    WX = 0.0;
    WY = 0.0;
    arcftCGSpd = Vector3d(WX, WY, -WZ);
  }

  // REGION 3: in the vortex material
  else if (R1PSML < dbCoreRadius) {
    double lambda =
        R1PSML / dbCoreRadius;  // relative percent in vortex material
    // 涡环截面上的相对位置
    Vector3d insidePointPos = Vector3d(dbRingVortex * XCGDBC / RVAXCG,
                                       dbRingVortex * YCGDBC / RVAXCG, 0.0);
    // 涡丝边界点xn相对上涡环的位置
    Vector3d outsidePointPosP =
        (cg2dbPos + (lambda - 1) * insidePointPos) / lambda;
    // 将该点转至地轴系
    Vector3d outsidePointPosGround =
        matp.transpose() * outsidePointPosP + dbPosP;
    // xn相对下涡环的位置
    Vector3d outsidePointPosM = matm * (outsidePointPosGround - dbPosM);
    // xn的速度
    calPointSpd(outsidePointPosP, outsidePointPosM, &arcftCGSpd);
    // 通过lambda缩放得到实际速度
    arcftCGSpd *= lambda;
  } else {
    // REGION2: common
    Vector3d relativePosP = matp * (arcftCGPos - dbPosP);
    Vector3d relativePosM = matm * (arcftCGPos - dbPosM);
    // relative position to mirror vortex
    calPointSpd(relativePosP, relativePosM, &arcftCGSpd);
  }
  return 0;
}

int downbrust::STRMFN_partialPart(Vector3d relativePosinput,
                                  double* phi_r_zoutput) {
  // 对流函数求偏导数，参考赵燕勤的matlab程序，但导数公式没看懂，使用定义和mathematica重新推导
  // 在matlab中检查计算结果。绘图无误
  double xinput = relativePosinput[0];
  double yinput = relativePosinput[1];
  double zinput = relativePosinput[2];
  double rinput = sqrt(SQ(xinput) + SQ(yinput));

  // r1 and r2 of primary and mirror vortex
  double r1 = sqrt(SQ(zinput) + SQ(rinput - dbRingVortex));
  double r2 = sqrt(SQ(zinput) + SQ(rinput + dbRingVortex));
  double r1_z = zinput / r1;
  double r2_z = zinput / r2;
  double r1_r = (rinput - dbRingVortex) / r1;
  double r2_r = (rinput + dbRingVortex) / r2;

  // use mathematica for Differential phi, and check output in matlab
  double phi_r1 =
      ((-r1 + r2) * (-((0.75 * SQ(r2) * sqrt((r1 * r2) / SQ(r1 + r2))) / r1) +
                     r2 * (-0.75 - 3. * sqrt((r1 * r2) / SQ(r1 + r2))) +
                     r1 * (-0.25 - 2.25 * sqrt((r1 * r2) / SQ(r1 + r2))))) /
      (SQ(r1 + r2) * SQ(0.25 + 1.5 * sqrt((r1 * r2) / SQ(r1 + r2))));
  double phi_r2 =
      ((-r1 + r2) * ((0.75 * SQ(r1) * sqrt((r1 * r2) / SQ(r1 + r2))) / r2 +
                     r2 * (0.25 + 2.25 * sqrt((r1 * r2) / SQ(r1 + r2))) +
                     r1 * (0.75 + 3. * sqrt((r1 * r2) / SQ(r1 + r2))))) /
      (SQ(r1 + r2) * SQ(0.25 + 1.5 * sqrt((r1 * r2) / SQ(r1 + r2))));

  phi_r_zoutput[0] = phi_r1 * r1_r + phi_r2 * r2_r;
  phi_r_zoutput[1] = phi_r1 * r1_z + phi_r2 * r2_z;
  return 0;
}

int downbrust::calPointSpd(Vector3d posP, Vector3d posM, Vector3d* actualSpd) {
  // 给定一个点的坐标，计算该点的速度
  double STRMFNConst = -circrv / 2 / MATH::PI * 0.788;
  double STRMFNVaribP[2], STRMFNVaribM[2];
  STRMFN_partialPart(posP, STRMFNVaribP);
  STRMFN_partialPart(posM, STRMFNVaribM);

  double XCGDBC, YCGDBC, RVAXCG;
  XCGDBC = posP[0];
  YCGDBC = posP[1];
  RVAXCG = sqrt(SQ(XCGDBC) + SQ(YCGDBC));

  double WPX = STRMFNConst * (STRMFNVaribP[1]) * XCGDBC / RVAXCG / RVAXCG;
  double WPY = STRMFNConst * (STRMFNVaribP[1]) * YCGDBC / RVAXCG / RVAXCG;
  double WPZ = STRMFNConst * (STRMFNVaribP[0]) / RVAXCG;

  XCGDBC = posM[0];
  YCGDBC = posM[1];
  RVAXCG = sqrt(SQ(XCGDBC) + SQ(YCGDBC));
  double WMX = STRMFNConst * (STRMFNVaribM[1]) * XCGDBC / RVAXCG / RVAXCG;
  double WMY = STRMFNConst * (STRMFNVaribM[1]) * YCGDBC / RVAXCG / RVAXCG;
  double WMZ = STRMFNConst * (STRMFNVaribM[0]) / RVAXCG;

  Vector3d WP, WM;
  WP << WPX, WPY, WPZ;
  WM << WMX, WMY, WMZ;
  *actualSpd = matp.transpose() * WP - matm.transpose() * WM;
  return 0;
}

// int downbrust::STRMFN(double hinput, double rinput, double* strfouptut) {
//   double RDVORT = dbRingVortex;
//   double HGVORT = dbPos[2];
//   double HGCG = hinput;
//   double RVAXCG = rinput;

//   // r1 and r2 of primary and mirror vortex
//   double R1PSML = sqrt(SQ(HGCG - HGVORT) + SQ(RVAXCG - RDVORT));
//   double R1MSML = sqrt(SQ(HGCG + HGVORT) + SQ(RVAXCG - RDVORT));
//   double R2PLRG = sqrt(SQ(HGCG - HGVORT) + SQ(RVAXCG + RDVORT));
//   double R2MLRG = sqrt(SQ(HGCG + HGVORT) + SQ(RVAXCG + RDVORT));

//   // k-mod of ellipse integral
//   double KMODPV = (R2PLRG - R1PSML) / (R2PLRG + R1PSML);
//   double ACMPEP = (0.788 * SQ(KMODPV)) / (0.25 + 0.75 * sqrt(1 -SQ(KMODPV)));
//
//   double KMODMV = (R2MLRG - R1MSML) / (R2MLRG + R1MSML); double
//   ACMPEM = (0.788 * SQ(KMODMV)) / (0.25 + 0.75 * sqrt(1 - SQ(KMODMV)));

//   // steam funtion
//   double STRMFN = -circrv / 2 / MATH::PI *
//                   ((R1PSML + R2PLRG) * ACMPEP - (R1MSML + R2MLRG) * ACMPEM);
//   *strfouptut = STRMFN;
//   return 0;
// }

// int downbrust::calSpeed(Vector3d relativePosinput, double* woutput) {
//   double xinput = relativePosinput[0];
//   double yinput = relativePosinput[0];
//   double hinput = relativePosinput[0];
//   double rinput = sqrt(SQ(xinput) + SQ(yinput));
//   double RVAXCG = rinput;
//   double step = 0.01;
//   double STRFCG, STRFCZ, STRFCR;
//   STRMFN(hinput, rinput, &STRFCG);
//   STRMFN(hinput - step, rinput, &STRFCZ);
//   STRMFN(hinput, rinput + step, &STRFCR);

//   double wroutput = (STRFCZ - STRFCG) / RVAXCG;  // wr
//   woutput[0] = wroutput * xinput / RVAXCG;
//   woutput[1] = wroutput * yinput / RVAXCG;
//   woutput[2] = (STRFCG - STRFCR) / RVAXCG;  // wz
//   return 0;
// }

}  // namespace DOWNBRUST