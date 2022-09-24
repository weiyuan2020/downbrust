/*
 * @Author: L.F.Wang;
 * @Date: 2021-08-01 20:50:37 
 * @Last Modified by: weiyuan
 * @Last Modified time: 2021-08-01 20:52:44
 * @Description:
 */

#pragma once
#include "MathDef.h"
#include "EigenInc.h"
#define Dcm21ToElrang(ORDER) Dcm21ToElrang_##ORDER
#define ElrangToDcm21(ORDER) ElrangToDcm21_##ORDER
#define ElrangToDcm12(ORDER) ElrangToDcm12_##ORDER
#define CalcElrateMtx(AX_Type, ORDER) CalcElrateMtx##AX_Type##_##ORDER
#define CalcElrateMtxDeri(AX_Type, ORDER) CalcElrateMtxDeri##AX_Type##_##ORDER
using Eigen::Matrix3d;
using Eigen::Vector3d;
namespace MATH
{

    Matrix3d skew(const Vector3d &vec3d);
    double atan_4quad(double y, double x);

    int Dcm21ToElrang_123(const Matrix3d &InDcm_Ax1ToAx2, Vector3d *OutElrang_Ax1ToAx2);
    Vector3d Dcm21ToElrang_123(const Matrix3d &InDcm_Ax1ToAx2);
    int ElrangToDcm21_123(const Vector3d &InElrang_Ax1ToAx2, Matrix3d *OutDcm_Ax1ToAx2);
    Matrix3d ElrangToDcm21_123(const Vector3d &InElrang_Ax1ToAx2);
    int ElrangToDcm12_123(const Vector3d &InElrang_Ax1ToAx2, Matrix3d *OutDcm_Ax2ToAx1);
    Matrix3d ElrangToDcm12_123(const Vector3d &InElrang_Ax1ToAx2);
    int CalcElrateMtxPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1, Matrix3d *OutKintcMtxAx1Inv);
    int CalcElrateMtxDeriPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1Deri);
    int CalcElrateMtxAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2, Matrix3d *OutKintcMtxAx2Inv);
    int CalcElrateMtxDeriAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2Deri);
    Matrix3d CalcElrateMtxPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxInvPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxDeriPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2);
    Matrix3d CalcElrateMtxAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxInvAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxDeriAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2);

    int Dcm21ToElrang_213(const Matrix3d &InDcm_Ax1ToAx2, Vector3d *OutElrang_Ax1ToAx2);
    Vector3d Dcm21ToElrang_213(const Matrix3d &InDcm_Ax1ToAx2);
    int ElrangToDcm21_213(const Vector3d &InElrang_Ax1ToAx2, Matrix3d *OutDcm_Ax1ToAx2);
    Matrix3d ElrangToDcm21_213(const Vector3d &InElrang_Ax1ToAx2);
    int CalcElrateMtxPreAx_213(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1, Matrix3d *OutKintcMtxAx1Inv);
    int CalcElrateMtxDeriPreAx_213(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1Deri);
    int CalcElrateMtxAftAx_213(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2, Matrix3d *OutKintcMtxAx2Inv);
    int CalcElrateMtxDeriAftAx_213(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2Deri);
    Matrix3d CalcElrateMtxPreAx_213(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxInvPreAx_213(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxDeriPreAx_213(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2);
    Matrix3d CalcElrateMtxAftAx_213(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxInvAftAx_213(
        const Vector3d &InElrang_Ax1ToAx2);
    Matrix3d CalcElrateMtxDeriAftAx_213(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2);

}
