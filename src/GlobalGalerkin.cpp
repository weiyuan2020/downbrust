/*
 * @Author: L.F.Wang
 * @Date: 2021-08-01 21:04:39
 * @Last Modified by:   weiyuan
 * @Last Modified time: 2021-08-01 21:04:39
 */
 /*****************************************************
 Name: GlobalGalerkin.cpp
 Copyright:
 Author: L.F.Wang; Y.Wei
 Date: 16/04/21 09:32
 Description:
 ******************************************************/
#include "GlobalGalerkin.h"

namespace MATH {
    //CLASS_GlobalGalerkin::CLASS_GlobalGalerkin(
    //	int _HarmoOrder, int _VarNum
    //) :
    //	HarmoOrder(_HarmoOrder),
    //	VarNum(_VarNum)
    //{
    //	CoeffMtx = MatrixXd::Zero(VarNum, 2 * HarmoOrder + 1);
    //	FourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
    //	dFourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
    //	d2FourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
    //}
    CLASS_GlobalGalerkin::CLASS_GlobalGalerkin() {
    }

    CLASS_GlobalGalerkin::~CLASS_GlobalGalerkin() {
    }

    int CLASS_GlobalGalerkin::Init( // initial function
        int _HarmoOrder,			//
        int _VarNum					//
    ) {
        HarmoOrder = _HarmoOrder;
        VarNum = _VarNum;
        CoeffMtx = MatrixXd::Zero(VarNum, 2 * HarmoOrder + 1);
        FourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
        dFourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
        d2FourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
        // 	HarmoOrder = _HarmoOrder;
        // VarNum = _VarNum;
        // CoeffMtx = MatrixXd::Zero(VarNum, 2 * HarmoOrder + 1);
        // (1, cos(Psi), sin(kPsi),..., cos(kPsi), sin(kPsi))
        // FourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
        // dFourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);
        // d2FourierBaseVec = VectorXd::Zero(2 * HarmoOrder + 1);

        return 0;
    }
    int CLASS_GlobalGalerkin::SetCoeffMtxAndBaseFrq(
        double _BaseFrq,
        const double *VarRowArr
    ) {
        BaseFrq = _BaseFrq;
        //std::cout << "BaseFrq: " << std::endl << BaseFrq << std::endl;
        for (int i = 0; i < VarNum; ++i) {
            for (int j = 0; j < 2 * HarmoOrder + 1; ++j) {
                CoeffMtx(i, j) = VarRowArr[i * (2 * HarmoOrder + 1) + j];
            }
        }
        //std::cout << "CoeffMtx = " << CoeffMtx << '\n';
        return 0;
    }

    int CLASS_GlobalGalerkin::SetTime(double t) {
        //std::cout << "in global Garlerkin FourierBaseVec: " << std::endl;
        FourierBaseVec(0) = 1.0;
        dFourierBaseVec(0) = 0.0;
        d2FourierBaseVec(0) = 0.0;
        for (int i = 0; i < HarmoOrder; ++i) {
            FourierBaseVec[2 * i + 1] = std::sin((i + 1) * BaseFrq * t);
            FourierBaseVec[2 * i + 2] = std::cos((i + 1) * BaseFrq * t);
            dFourierBaseVec[2 * i + 1] = (i + 1) * BaseFrq * FourierBaseVec[2 * i + 2];
            dFourierBaseVec[2 * i + 2] = -(i + 1) * BaseFrq * FourierBaseVec[2 * i + 1];
            d2FourierBaseVec[2 * i + 1] = -(i + 1) * BaseFrq * (i + 1) * BaseFrq * FourierBaseVec[2 * i + 1];
            d2FourierBaseVec[2 * i + 2] = -(i + 1) * BaseFrq * (i + 1) * BaseFrq * FourierBaseVec[2 * i + 2];
        }
        return 0;
    }

    Eigen::VectorXd CLASS_GlobalGalerkin::GetFourierBaseVec() {
        return FourierBaseVec;
    }

    Eigen::VectorXd CLASS_GlobalGalerkin::GetVarVal() {
        return CoeffMtx * FourierBaseVec;
    }

    Eigen::VectorXd CLASS_GlobalGalerkin::GetVarVelo() {
        return CoeffMtx * dFourierBaseVec;
    }

    Eigen::VectorXd CLASS_GlobalGalerkin::GetVarAcc() {
        return CoeffMtx * d2FourierBaseVec;
    }

    int CLASS_GlobalGalerkin::GetFourierBaseVec(Eigen::VectorXd *FourierBaseVec_) {
        *FourierBaseVec_ = FourierBaseVec;
        return 0;
    }

    int CLASS_GlobalGalerkin::GetVarVal(Eigen::VectorXd *VarVal_) {
        *VarVal_ = CoeffMtx * FourierBaseVec;
        return 0;
    }

    int CLASS_GlobalGalerkin::GetVarVelo(Eigen::VectorXd *VarVelo_) {
        *VarVelo_ = CoeffMtx * dFourierBaseVec;
        return 0;
    }

    int CLASS_GlobalGalerkin::GetVarAcc(Eigen::VectorXd *VarAcc_) {
        *VarAcc_ = CoeffMtx * d2FourierBaseVec;
        return 0;
    }
}