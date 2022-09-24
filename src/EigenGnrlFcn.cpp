/*
 * @Author: L.F.Wang
 * @Date: 2021-08-01 20:52:40
 * @Last Modified by: weiyuan
 * @Last Modified time: 2021-08-27 21:44:25
 */
#include "EigenGnrlFcn.h"
namespace MATH {
    Matrix3d skew(const Vector3d &vec3d) {
        Eigen::Matrix3d mtx3d;
        mtx3d << 0.0, -vec3d(2), vec3d(1),
            vec3d(2), 0.0, -vec3d(0),
            -vec3d(1), vec3d(0), 0.0;
        return mtx3d;
    }

    double atan_4quad(double y, double x) {
        if (x >= 0.0) {
            return atan2(y, x);
        }
        else if (y >= 0.0) {
            return atan2(y, x) + PI;
        }
        else {
            return atan2(y, x) - PI;
        }
    }

    int Dcm21ToElrang_123(const Matrix3d &InDcm_Ax1ToAx2, Vector3d *OutElrang_Ax1ToAx2) {
        double phi, the, psi;

        //-90 <= theta <= +90
        the = -std::asin(InDcm_Ax1ToAx2(0, 2));
        double CosThe = std::cos(the);
        //-180 < phi <= +180
        double CosPhi = InDcm_Ax1ToAx2(2, 2) / CosThe;
        if (CosPhi >= 0.0)
            phi = std::asin(InDcm_Ax1ToAx2(1, 2) / CosThe);
        else {
            double SinPhi = std::asin(InDcm_Ax1ToAx2(1, 2) / CosThe);
            double SgnSinPhi = SinPhi >= 0.0 ? 1.0 : -1.0;
            phi = SinPhi + SgnSinPhi * PI / 2.0;
        }
        //-180 < psi <= +180
        double CosPsi = InDcm_Ax1ToAx2(0, 0) / CosThe;
        if (CosPsi >= 0.0)
            psi = std::asin(InDcm_Ax1ToAx2(0, 1) / CosThe);
        else {
            double SinPsi = std::asin(InDcm_Ax1ToAx2(0, 1) / CosThe);
            double SgnSinPsi = SinPsi >= 0.0 ? 1.0 : -1.0;
            psi = SinPsi + SgnSinPsi * PI / 2.0;
        }
        (*OutElrang_Ax1ToAx2)(0) = phi;
        (*OutElrang_Ax1ToAx2)(1) = the;
        (*OutElrang_Ax1ToAx2)(2) = psi;
        return 0;
    }

    Eigen::Vector3d Dcm21ToElrang_123(const Matrix3d &InDcm_Ax1ToAx2) {
        double phi, the, psi;

        //-90 <= theta <= +90
        the = -std::asin(InDcm_Ax1ToAx2(0, 2));
        double CosThe = std::cos(the);
        //-180 < phi <= +180
        double CosPhi = InDcm_Ax1ToAx2(2, 2) / CosThe;
        if (CosPhi >= 0.0)
            phi = std::asin(InDcm_Ax1ToAx2(1, 2) / CosThe);
        else {
            double SinPhi = std::asin(InDcm_Ax1ToAx2(1, 2) / CosThe);
            double SgnSinPhi = SinPhi >= 0.0 ? 1.0 : -1.0;
            phi = SinPhi + SgnSinPhi * PI / 2.0;
        }
        //-180 < psi <= +180
        double CosPsi = InDcm_Ax1ToAx2(0, 0) / CosThe;
        if (CosPsi >= 0.0)
            psi = std::asin(InDcm_Ax1ToAx2(0, 1) / CosThe);
        else {
            double SinPsi = std::asin(InDcm_Ax1ToAx2(0, 1) / CosThe);
            double SgnSinPsi = SinPsi >= 0.0 ? 1.0 : -1.0;
            psi = SinPsi + SgnSinPsi * PI / 2.0;
        }

        return Vector3d(phi, the, psi);
    }

    int ElrangToDcm21_123(const Vector3d &InElrang_Ax1ToAx2, Matrix3d *OutDcm_Ax1ToAx2) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        (*OutDcm_Ax1ToAx2)(0, 0) = cThe * cPsi;
        (*OutDcm_Ax1ToAx2)(0, 1) = cThe * sPsi;
        (*OutDcm_Ax1ToAx2)(0, 2) = -sThe;

        (*OutDcm_Ax1ToAx2)(1, 0) = sThe * sPhi * cPsi - cPhi * sPsi;
        (*OutDcm_Ax1ToAx2)(1, 1) = sThe * sPhi * sPsi + cPhi * cPsi;
        (*OutDcm_Ax1ToAx2)(1, 2) = cThe * sPhi;

        (*OutDcm_Ax1ToAx2)(2, 0) = sThe * cPhi * cPsi + sPhi * sPsi;
        (*OutDcm_Ax1ToAx2)(2, 1) = sThe * cPhi * sPsi - sPhi * cPsi;
        (*OutDcm_Ax1ToAx2)(2, 2) = cThe * cPhi;

        return 0;
    }

    Eigen::Matrix3d ElrangToDcm21_123(const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        Matrix3d OutDcm_Ax1ToAx2;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        OutDcm_Ax1ToAx2(0, 0) = cThe * cPsi;
        OutDcm_Ax1ToAx2(0, 1) = cThe * sPsi;
        OutDcm_Ax1ToAx2(0, 2) = -sThe;

        OutDcm_Ax1ToAx2(1, 0) = sThe * sPhi * cPsi - cPhi * sPsi;
        OutDcm_Ax1ToAx2(1, 1) = sThe * sPhi * sPsi + cPhi * cPsi;
        OutDcm_Ax1ToAx2(1, 2) = cThe * sPhi;

        OutDcm_Ax1ToAx2(2, 0) = sThe * cPhi * cPsi + sPhi * sPsi;
        OutDcm_Ax1ToAx2(2, 1) = sThe * cPhi * sPsi - sPhi * cPsi;
        OutDcm_Ax1ToAx2(2, 2) = cThe * cPhi;

        return OutDcm_Ax1ToAx2;
    }

    int ElrangToDcm12_123(const Vector3d &InElrang_Ax1ToAx2, Matrix3d *OutDcm_Ax2ToAx1) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        (*OutDcm_Ax2ToAx1)(0, 0) = cThe * cPsi;
        (*OutDcm_Ax2ToAx1)(1, 0) = cThe * sPsi;
        (*OutDcm_Ax2ToAx1)(2, 0) = -sThe;

        (*OutDcm_Ax2ToAx1)(1, 1) = sThe * sPhi * cPsi - cPhi * sPsi;
        (*OutDcm_Ax2ToAx1)(1, 1) = sThe * sPhi * sPsi + cPhi * cPsi;
        (*OutDcm_Ax2ToAx1)(2, 1) = cThe * sPhi;

        (*OutDcm_Ax2ToAx1)(0, 2) = sThe * cPhi * cPsi + sPhi * sPsi;
        (*OutDcm_Ax2ToAx1)(1, 2) = sThe * cPhi * sPsi - sPhi * cPsi;
        (*OutDcm_Ax2ToAx1)(2, 2) = cThe * cPhi;

        return 0;
    }

    Matrix3d ElrangToDcm12_123(const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        Matrix3d OutDcm_Ax2ToAx1;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        OutDcm_Ax2ToAx1(0, 0) = cThe * cPsi;
        OutDcm_Ax2ToAx1(1, 0) = cThe * sPsi;
        OutDcm_Ax2ToAx1(2, 0) = -sThe;

        OutDcm_Ax2ToAx1(0, 1) = sThe * sPhi * cPsi - cPhi * sPsi;
        OutDcm_Ax2ToAx1(1, 1) = sThe * sPhi * sPsi + cPhi * cPsi;
        OutDcm_Ax2ToAx1(2, 1) = cThe * sPhi;

        OutDcm_Ax2ToAx1(0, 2) = sThe * cPhi * cPsi + sPhi * sPsi;
        OutDcm_Ax2ToAx1(1, 2) = sThe * cPhi * sPsi - sPhi * cPsi;
        OutDcm_Ax2ToAx1(2, 2) = cThe * cPhi;

        return OutDcm_Ax2ToAx1;
    }

    int CalcElrateMtxPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1, Matrix3d *OutKintcMtxAx1Inv) {
        using std::cos;
        using std::sin;
        using std::tan;

        //     double sPhi = sin( InElrang_Ax1ToAx2( 0 ) );
        //     double cPhi = cos( InElrang_Ax1ToAx2( 0 ) );
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        if (fabs(cThe) < MINRAD) {
            return 1;
        }
        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax1
        (*OutKintcMtxAx1) << cThe * cPsi, -sPsi, 0.0,
            cThe * sPsi, cPsi, 0.0,
            -sThe, 0.0, 1.0;

        (*OutKintcMtxAx1Inv) << cPsi / cThe, sPsi / cThe, 0.0,
            -sPsi, cPsi, 0.0,
            sThe * cPsi / cThe, sThe * sPsi / cThe, 1.0;

        return 0;
    }

    int CalcElrateMtxDeriPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1Deri) {
        using std::cos;
        using std::sin;

        //     double sPhi = sin( InElrang_Ax1ToAx2( 0 ) );
        //     double cPhi = cos( InElrang_Ax1ToAx2( 0 ) );
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        (*OutKintcMtxAx1Deri) << -InElrangVelo_Ax1ToAx2(1) * sThe * cPsi -
            InElrangVelo_Ax1ToAx2(2) * cThe * sPsi,
            -InElrangVelo_Ax1ToAx2(2) * cPsi,
            0.0,
            -InElrangVelo_Ax1ToAx2(1) * sThe * sPsi +
            InElrangVelo_Ax1ToAx2(2) * cThe * cPsi,
            -InElrangVelo_Ax1ToAx2(2) * sPsi,
            0.0,
            -InElrangVelo_Ax1ToAx2(1) * cThe,
            0.0,
            0.0;
        return 0;
    }

    int CalcElrateMtxAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2, Matrix3d *OutKintcMtxAx2Inv) {
        using std::cos;
        using std::sin;
        using std::tan;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        if (fabs(cThe) < MINRAD) {
            return 1;
        }
        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax2
        (*OutKintcMtxAx2) << 1.0, 0.0, -sThe,
            0.0, cPhi, sPhi * cThe,
            0.0, -sPhi, cPhi * cThe;

        (*OutKintcMtxAx2Inv) << 1.0, sPhi * sThe / cThe, cPhi * sThe / cThe,
            0.0, cPhi, -sPhi,
            0.0, sPhi / cThe, cPhi / cThe;
        return 0;
    }

    int CalcElrateMtxDeriAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2Deri) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        (*OutKintcMtxAx2Deri) << 0.0,
            0.0,
            -InElrangVelo_Ax1ToAx2(1) * cThe,
            0.0,
            -InElrangVelo_Ax1ToAx2(0) * sPhi,
            InElrangVelo_Ax1ToAx2(0) * cPhi * cThe -
            InElrangVelo_Ax1ToAx2(1) * sPhi * sThe,
            0.0,
            -InElrangVelo_Ax1ToAx2(0) * cPhi,
            -InElrangVelo_Ax1ToAx2(0) * sPhi * cThe -
            InElrangVelo_Ax1ToAx2(1) * cPhi * sThe;
        return 0;
    }

    Matrix3d CalcElrateMtxPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx1;

        //     double sPhi = sin( InElrang_Ax1ToAx2( 0 ) );
        //     double cPhi = cos( InElrang_Ax1ToAx2( 0 ) );
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax1
        OutKintcMtxAx1 << cThe * cPsi, -sPsi, 0.0,
            cThe * sPsi, cPsi, 0.0,
            -sThe, 0.0, 1.0;

        return OutKintcMtxAx1;
    }

    Matrix3d CalcElrateMtxInvPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx1Inv;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        if (fabs(cThe) < MINRAD) {
            return Matrix3d::Zero();
        }
        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax1
        OutKintcMtxAx1Inv << cPsi / cThe, sPsi / cThe, 0.0,
            -sPsi, cPsi, 0.0,
            sThe / cThe * cPsi, sThe / cThe * sPsi, 1.0;

        return OutKintcMtxAx1Inv;
    }

    Matrix3d CalcElrateMtxDeriPreAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2) {
        using std::cos;
        using std::sin;

        Matrix3d OutKintcMtxAx1Deri;

        //     double sPhi = sin( InElrang_Ax1ToAx2( 0 ) );
        //     double cPhi = cos( InElrang_Ax1ToAx2( 0 ) );
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        OutKintcMtxAx1Deri << -InElrangVelo_Ax1ToAx2(1) * sThe * cPsi -
            InElrangVelo_Ax1ToAx2(2) * cThe * sPsi,
            -InElrangVelo_Ax1ToAx2(2) * cPsi,
            0.0,
            -InElrangVelo_Ax1ToAx2(1) * sThe * sPsi +
            InElrangVelo_Ax1ToAx2(2) * cThe * cPsi,
            -InElrangVelo_Ax1ToAx2(2) * sPsi,
            0.0,
            -InElrangVelo_Ax1ToAx2(1) * cThe,
            0.0,
            0.0;
        return OutKintcMtxAx1Deri;
    }

    Matrix3d CalcElrateMtxAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx2;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax2
        OutKintcMtxAx2 << 1.0, 0.0, -sThe,
            0.0, cPhi, sPhi * cThe,
            0.0, -sPhi, cPhi * cThe;

        return OutKintcMtxAx2;
    }

    Matrix3d CalcElrateMtxInvAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx2Inv;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        if (fabs(cThe) < MINRAD) {
            return Matrix3d::Zero();
        }
        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax2
        OutKintcMtxAx2Inv << 1.0, sPhi * sThe / cThe, cPhi * sThe / cThe,
            0.0, cPhi, -sPhi,
            0.0, sPhi / cThe, cPhi / cThe;
        return OutKintcMtxAx2Inv;
    }

    Matrix3d CalcElrateMtxDeriAftAx_123(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2) {
        using std::cos;
        using std::sin;

        Matrix3d OutKintcMtxAx2Deri;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        OutKintcMtxAx2Deri << 0.0,
            0.0,
            -InElrangVelo_Ax1ToAx2(1) * cThe,
            0.0,
            -InElrangVelo_Ax1ToAx2(0) * sPhi,
            InElrangVelo_Ax1ToAx2(0) * cPhi * cThe -
            InElrangVelo_Ax1ToAx2(1) * sPhi * sThe,
            0.0,
            -InElrangVelo_Ax1ToAx2(0) * cPhi,
            -InElrangVelo_Ax1ToAx2(0) * sPhi * cThe -
            InElrangVelo_Ax1ToAx2(1) * cPhi * sThe;
        return OutKintcMtxAx2Deri;
    }

    int Dcm21ToElrang_213(const Matrix3d &InDcm_Ax1ToAx2, Vector3d *OutElrang_Ax1ToAx2) {
        using std::asin;

        (*OutElrang_Ax1ToAx2)(0) = asin(InDcm_Ax1ToAx2(1, 2));
        (*OutElrang_Ax1ToAx2)(1) = atan_4quad(-InDcm_Ax1ToAx2(0, 2), InDcm_Ax1ToAx2(2, 2));
        (*OutElrang_Ax1ToAx2)(2) = atan_4quad(-InDcm_Ax1ToAx2(1, 0), InDcm_Ax1ToAx2(1, 1));
        return 0;
    }

    Eigen::Vector3d Dcm21ToElrang_213(const Matrix3d &InDcm_Ax1ToAx2) {
        using std::asin;

        return Vector3d(
            asin(InDcm_Ax1ToAx2(1, 2)),
            atan_4quad(-InDcm_Ax1ToAx2(0, 2), InDcm_Ax1ToAx2(2, 2)),
            atan_4quad(-InDcm_Ax1ToAx2(1, 0), InDcm_Ax1ToAx2(1, 1)));
    }

    int ElrangToDcm21_213(const Vector3d &InElrang_Ax1ToAx2, Matrix3d *OutDcm_Ax1ToAx2) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        (*OutDcm_Ax1ToAx2)(0, 0) = cThe * cPsi - sThe * sPhi * sPsi;
        (*OutDcm_Ax1ToAx2)(0, 1) = cThe * sPsi + sThe * sPhi * cPsi;
        (*OutDcm_Ax1ToAx2)(0, 2) = -sThe * cPhi;

        (*OutDcm_Ax1ToAx2)(1, 0) = -cPhi * sPsi;
        (*OutDcm_Ax1ToAx2)(1, 1) = cPhi * cPsi;
        (*OutDcm_Ax1ToAx2)(1, 2) = sPhi;

        (*OutDcm_Ax1ToAx2)(2, 0) = sThe * cPsi + cThe * sPhi * sPsi;
        (*OutDcm_Ax1ToAx2)(2, 1) = sThe * sPsi - cThe * sPhi * cPsi;
        (*OutDcm_Ax1ToAx2)(2, 2) = cThe * cPhi;

        return 0;
    }

    Eigen::Matrix3d ElrangToDcm21_213(const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        Matrix3d OutDcm_Ax1ToAx2;
        OutDcm_Ax1ToAx2(0, 0) = cThe * cPsi - sThe * sPhi * sPsi;
        OutDcm_Ax1ToAx2(0, 1) = cThe * sPsi + sThe * sPhi * cPsi;
        OutDcm_Ax1ToAx2(0, 2) = -sThe * cPhi;

        OutDcm_Ax1ToAx2(1, 0) = -cPhi * sPsi;
        OutDcm_Ax1ToAx2(1, 1) = cPhi * cPsi;
        OutDcm_Ax1ToAx2(1, 2) = sPhi;

        OutDcm_Ax1ToAx2(2, 0) = sThe * cPsi + cThe * sPhi * sPsi;
        OutDcm_Ax1ToAx2(2, 1) = sThe * sPsi - cThe * sPhi * cPsi;
        OutDcm_Ax1ToAx2(2, 2) = cThe * cPhi;

        return OutDcm_Ax1ToAx2;
    }

    int CalcElrateMtxPreAx_213(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1, Matrix3d *OutKintcMtxAx1Inv) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        //     double sThe = sin( InElrang_Ax1ToAx2( 1 ) );
        //     double cThe = cos( InElrang_Ax1ToAx2( 1 ) );
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        (*OutKintcMtxAx1)(0, 0) = cPsi;
        (*OutKintcMtxAx1)(0, 1) = -cPhi * sPsi;
        (*OutKintcMtxAx1)(0, 2) = 0;
        (*OutKintcMtxAx1)(1, 0) = sPsi;
        (*OutKintcMtxAx1)(1, 1) = cPhi * cPsi;
        (*OutKintcMtxAx1)(1, 2) = 0;
        (*OutKintcMtxAx1)(2, 0) = 0;
        (*OutKintcMtxAx1)(2, 1) = sPhi;
        (*OutKintcMtxAx1)(2, 2) = 1;

        if (fabs(cPhi) < MINRAD)
            return 1;

        (*OutKintcMtxAx1Inv)(0, 0) = cPsi;
        (*OutKintcMtxAx1Inv)(0, 1) = sPsi;
        (*OutKintcMtxAx1Inv)(0, 2) = 0;
        (*OutKintcMtxAx1Inv)(1, 0) = -sPsi / cPhi;
        (*OutKintcMtxAx1Inv)(1, 1) = cPsi / cPhi;
        (*OutKintcMtxAx1Inv)(1, 2) = 0;
        (*OutKintcMtxAx1Inv)(2, 0) = sPsi * sPhi / cPhi;
        (*OutKintcMtxAx1Inv)(2, 1) = -cPsi * sPhi / cPhi;
        (*OutKintcMtxAx1Inv)(2, 2) = 1;

        return 0;
    }

    int CalcElrateMtxDeriPreAx_213(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx1Deri) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        //     double sThe = sin( InElrang_Ax1ToAx2( 1 ) );
        //     double cThe = cos( InElrang_Ax1ToAx2( 1 ) );
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        double dphi = InElrangVelo_Ax1ToAx2(0);
        double dpsi = InElrangVelo_Ax1ToAx2(2);

        (*OutKintcMtxAx1Deri)(0, 0) = -sPsi * dpsi;
        (*OutKintcMtxAx1Deri)(0, 1) = sPhi * sPsi * dphi - cPhi * cPsi * dpsi;
        (*OutKintcMtxAx1Deri)(0, 2) = 0;
        (*OutKintcMtxAx1Deri)(1, 0) = cPsi * dpsi;
        (*OutKintcMtxAx1Deri)(1, 1) = -sPhi * cPsi * dphi - cPhi * sPsi * dpsi;
        (*OutKintcMtxAx1Deri)(1, 2) = 0;
        (*OutKintcMtxAx1Deri)(2, 0) = 0;
        (*OutKintcMtxAx1Deri)(2, 1) = cPhi * dphi;
        (*OutKintcMtxAx1Deri)(2, 2) = 0;

        return 0;
    }

    int CalcElrateMtxAftAx_213(
        const Vector3d &InElrang_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2, Matrix3d *OutKintcMtxAx2Inv) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        (*OutKintcMtxAx2)(0, 0) = cThe;
        (*OutKintcMtxAx2)(0, 1) = 0;
        (*OutKintcMtxAx2)(0, 2) = -sThe * cPhi;
        (*OutKintcMtxAx2)(1, 0) = 0;
        (*OutKintcMtxAx2)(1, 1) = 1;
        (*OutKintcMtxAx2)(1, 2) = sPhi;
        (*OutKintcMtxAx2)(2, 0) = sThe;
        (*OutKintcMtxAx2)(2, 1) = 0;
        (*OutKintcMtxAx2)(2, 2) = cThe * cPhi;

        if (fabs(cPhi) < MINRAD)
            return 1;

        (*OutKintcMtxAx2Inv)(0, 0) = cThe;
        (*OutKintcMtxAx2Inv)(0, 1) = 0;
        (*OutKintcMtxAx2Inv)(0, 2) = sThe;
        (*OutKintcMtxAx2Inv)(1, 0) = sPhi * sThe / cPhi;
        (*OutKintcMtxAx2Inv)(1, 1) = 1;
        (*OutKintcMtxAx2Inv)(1, 2) = -sPhi * cThe / cPhi;
        (*OutKintcMtxAx2Inv)(2, 0) = -sThe / cPhi;
        (*OutKintcMtxAx2Inv)(2, 1) = 0;
        (*OutKintcMtxAx2Inv)(2, 2) = cThe / cPhi;

        return 2;
    }

    int CalcElrateMtxDeriAftAx_213(
        const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2,
        Matrix3d *OutKintcMtxAx2Deri) {
        using std::cos;
        using std::sin;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        double dphi = InElrangVelo_Ax1ToAx2(0);
        double dtheta = InElrangVelo_Ax1ToAx2(1);

        (*OutKintcMtxAx2Deri)(0, 0) = -sThe * dtheta;
        (*OutKintcMtxAx2Deri)(0, 1) = 0;
        (*OutKintcMtxAx2Deri)(0, 2) = -cThe * cPhi * dtheta + sThe * sPhi * dphi;
        (*OutKintcMtxAx2Deri)(1, 0) = 0;
        (*OutKintcMtxAx2Deri)(1, 1) = 0;
        (*OutKintcMtxAx2Deri)(1, 2) = cPhi * dphi;
        (*OutKintcMtxAx2Deri)(2, 0) = cThe * dtheta;
        (*OutKintcMtxAx2Deri)(2, 1) = 0;
        (*OutKintcMtxAx2Deri)(2, 2) = -sThe * cPhi * dtheta - cThe * sPhi * dphi;

        return 0;
    }

    Matrix3d CalcElrateMtxPreAx_213(const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx1;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        //     double sThe = sin( InElrang_Ax1ToAx2( 1 ) );
        //     double cThe = cos( InElrang_Ax1ToAx2( 1 ) );
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax1
        OutKintcMtxAx1 << cPsi, -cPhi * sPsi, 0,
            sPsi, cPhi * cPsi, 0,
            0, sPhi, 1;

        return OutKintcMtxAx1;
    }

    Matrix3d CalcElrateMtxInvPreAx_213(const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx1Inv;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        //     double sThe = sin( InElrang_Ax1ToAx2( 1 ) );
        //     double cThe = cos( InElrang_Ax1ToAx2( 1 ) );
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        if (InElrang_Ax1ToAx2(1) >= PI / 2.0 - MINRAD || InElrang_Ax1ToAx2(1) <= -PI / 2.0 + MINRAD) {
            return Matrix3d::Zero();
        }
        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax1
        OutKintcMtxAx1Inv << cPsi, sPsi, 0,
            -sPsi / cPhi, cPsi / cPhi, 0,
            sPsi * sPhi / cPhi, -cPsi * sPhi / cPhi, 1;
        return OutKintcMtxAx1Inv;
    }

    Matrix3d CalcElrateMtxDeriPreAx_213(const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2) {
        using std::cos;
        using std::sin;

        Matrix3d OutKintcMtxAx1Deri;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        double dphi = InElrangVelo_Ax1ToAx2(0);
        double dpsi = InElrangVelo_Ax1ToAx2(2);

        OutKintcMtxAx1Deri << -sPsi * dpsi,
            sPhi * sPsi * dphi - cPhi * cPsi * dpsi,
            0,
            cPsi * dpsi,
            -sPhi * cPsi * dphi - cPhi * sPsi * dpsi,
            0,
            0,
            cPhi * dphi,
            0;
        return OutKintcMtxAx1Deri;
    }

    Matrix3d CalcElrateMtxAftAx_213(const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx2;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax2
        OutKintcMtxAx2 << cThe,
            0,
            -sThe * cPhi,
            0,
            1,
            sPhi,
            sThe,
            0,
            cThe * cPhi;

        return OutKintcMtxAx2;
    }

    Matrix3d CalcElrateMtxInvAftAx_213(const Vector3d &InElrang_Ax1ToAx2) {
        using std::cos;
        using std::sin;
        using std::tan;

        Matrix3d OutKintcMtxAx2Inv;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        //     double sPsi = sin( InElrang_Ax1ToAx2( 2 ) );
        //     double cPsi = cos( InElrang_Ax1ToAx2( 2 ) );

        if (fabs(cPhi) < MINRAD) {
            return Matrix3d::Zero();
        }
        //KintcMtx convert Euler angular velocities to body fixed velocities in Ax2
        OutKintcMtxAx2Inv << cThe, 0, sThe,
            sPhi * sThe / cPhi, 1, -sPhi * cThe / cPhi,
            -sThe / cPhi, 0, cThe / cPhi;
        return OutKintcMtxAx2Inv;
    }

    Matrix3d CalcElrateMtxDeriAftAx_213(const Vector3d &InElrang_Ax1ToAx2, const Vector3d &InElrangVelo_Ax1ToAx2) {
        using std::cos;
        using std::sin;

        Matrix3d OutKintcMtxAx2Deri;

        double sPhi = sin(InElrang_Ax1ToAx2(0));
        double cPhi = cos(InElrang_Ax1ToAx2(0));
        double sThe = sin(InElrang_Ax1ToAx2(1));
        double cThe = cos(InElrang_Ax1ToAx2(1));
        double sPsi = sin(InElrang_Ax1ToAx2(2));
        double cPsi = cos(InElrang_Ax1ToAx2(2));

        double dphi = InElrangVelo_Ax1ToAx2(0);
        double dtheta = InElrangVelo_Ax1ToAx2(1);

        OutKintcMtxAx2Deri << -sThe * dtheta,
            0,
            -cThe * cPhi * dtheta + sThe * sPhi * dphi,
            0,
            0,
            cPhi * dphi,
            cThe * dtheta,
            0,
            -sThe * cPhi * dtheta - cThe * sPhi * dphi;
        return OutKintcMtxAx2Deri;
    }
}