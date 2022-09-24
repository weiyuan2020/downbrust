/*
 * @Author: L.F.Wang
 * @Date: 2021-08-01 20:57:15 
 * @Last Modified by: weiyuan
 * @Last Modified time: 2021-08-01 21:14:45
 * Description:
 */

#pragma once

#include <iostream>

#include "stdlib.h"
#include "EigenInc.h"

namespace MATH
{
	using namespace Eigen;
	class CLASS_GlobalGalerkin // 伽辽金类
	{
	private:
		//const int HarmoOrder;
		//const int VarNum;
		int HarmoOrder; // 伽辽金法阶数
		int VarNum;		// 输入数据数量 ex: (Flap, Lag) VarNum = 2 (Flap) VarNum = 1
		double BaseFrq; //

		MatrixXd CoeffMtx;
		VectorXd FourierBaseVec;
		VectorXd dFourierBaseVec;
		VectorXd d2FourierBaseVec;

		//VectorXd DiffCalcAndApprxVarAcc;

	public:
		CLASS_GlobalGalerkin();
		~CLASS_GlobalGalerkin();
		int Init(int _HarmoOrder, int _VarNum);
		int SetCoeffMtxAndBaseFrq(double _BaseFrq, const double *VarRowArr);
		int SetTime(double t);

		VectorXd GetFourierBaseVec();
		VectorXd GetVarVal();
		VectorXd GetVarVelo();
		VectorXd GetVarAcc();
		int GetFourierBaseVec(Eigen::VectorXd *FourierBaseVec_);
		int GetVarVal(Eigen::VectorXd *VarVal_);
		int GetVarVelo(Eigen::VectorXd *VarVelo_);
		int GetVarAcc(Eigen::VectorXd *VarAcc_);
	};

}
