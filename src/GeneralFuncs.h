/*
 * @Author:  L.F.Wang
 * @Date: 2021-08-01 20:56:47
 * @Last Modified by: weiyuan
 * @Last Modified time: 2021-08-01 20:57:08
 */

#pragma once
#include <iostream>
#include <math.h>
#include "MathDef.h"

namespace MATH {
	int toBinaryNum(int num, const int& len, char* bi);

	int JacobiMtxFDM(
		int Nx
		, int Nf
		, double *TrimVal
		, int(*obf)(int, int, double*, double*)
		, double **JacobiMtx
		, double *PerbFwd
		, double *PerbAft = NULL
	);

	int Gauss(
		int Dim
		, double **A
		, double *x
		, double MinGauss = 1.E-5
	);

	int GradDescent(
		int Nop
		, double *Value
		, double *Perb
		, int(*obf)(int, double*, double*)
		, double *FinalRemains
		, int MaxLoop = 40
		, int MaxRlx = 20
		, double MaxRemains = 1.E10
		, double MinRemains = 1.E-7
		, double MinPerb = 1.E-7
		, double MinGauss = 1.E-5
	);

	template <typename T1, typename T2>
	int TmpBinarySearch(int low, int high, T1* data, const T2& key)
	{
		int index = (low + high) / 2;
		if (key == data[index]) {
			return index;
		}
		else if (key < data[index]) {
			if (index == low + 1)
				return low;
			else
				return TmpBinarySearch<T1, T2>(low, index, data, key);
		}
		else {
			if (index == high - 1)
				return index;
			else
				return TmpBinarySearch<T1, T2>(index, high, data, key);
		}

	}

	template <typename ClassName>
	int TmpGradDescent(
		int Nop											  // var num
		, double *Value									  // init value
		, double *Perb									  // step size
		, ClassName *obj
		, int (ClassName::*obf)(int, double*, double*)
		, double *FinalRemains
		, int MaxLoop = 40
		, int MaxRlx = 20
		, double MaxRemains = 1.E10
		, double MinRemains = 1.E-7
		, double MinPerb = 1.E-7
		, double MinGauss = 1.E-5
	) {
		int i, j, looptimes, rlxtimes, ErrorNewt, ErrorObf;
		double remains1, remains2, rlx, SumValue, SumValueD;
		double *f, *f1, *f2, *dx, **dfdx;
		//memory space allocate
		f = (double*)malloc(sizeof(double)*Nop);
		f1 = (double*)malloc(sizeof(double)*Nop);
		f2 = (double*)malloc(sizeof(double)*Nop);
		dx = (double*)malloc(sizeof(double)*Nop);
		dfdx = (double**)malloc(sizeof(double*)*Nop);
		for (i = 0; i < Nop; ++i) {
			dfdx[i] = (double*)malloc(sizeof(double)*(Nop + 1));
		}
		//Initialize
		remains1 = MaxRemains;
		remains2 = 1.E10;
		looptimes = 0;
		ErrorObf = (obj->*obf)(Nop, Value, f);
		if (ErrorObf) {
			ErrorNewt = 10 + ErrorObf;
			goto EndNewtonRlx;
		}
		remains2 = 0.;
		for (i = 0; i < Nop; ++i)
			remains2 += SQ(f[i]);
		//loop start
		while (remains2 >= MinRemains) {
			rlx = 1.;
			rlxtimes = 0;
			looptimes += 1;
			if (looptimes > MaxLoop) {
				ErrorNewt = 2;//reach the max loop times
				goto EndNewtonRlx;
			}
			//rlx loop start
			while (remains2 > remains1*(1. - 0.25*rlx)) {
				if (looptimes == 1) {
					ErrorNewt = 3;//initial value is not good enough
					goto EndNewtonRlx;
				}
				rlxtimes += 1;
				if (rlxtimes > MaxRlx) {
					ErrorNewt = 4;//reach the max rlx times
					goto EndNewtonRlx;
				}
				remains1 = remains2;
				for (i = 0; i < Nop; ++i)
					Value[i] += rlx*dx[i];
				ErrorObf = (obj->*obf)(Nop, Value, f);
				if (ErrorObf) {
					ErrorNewt = 10 + ErrorObf;
					goto EndNewtonRlx;
				}
				remains2 = 0.;
				for (i = 0; i < Nop; ++i)
					remains2 += SQ(f[i]);
				if (remains2 < MinRemains) {
					ErrorNewt = 0;//achieve the desired accuracy in rlx loop, return successful
					goto EndNewtonRlx;
				}
			}
			for (i = 0; i < Nop; ++i) {
				Value[i] -= Perb[i];
				ErrorObf = (obj->*obf)(Nop, Value, f1);
				if (ErrorObf) {
					ErrorNewt = 10 + ErrorObf;
					goto EndNewtonRlx;
				}
				Value[i] += 2 * Perb[i];
				ErrorObf = (obj->*obf)(Nop, Value, f2);
				if (ErrorObf) {
					ErrorNewt = 10 + ErrorObf;
					goto EndNewtonRlx;
				}
				Value[i] -= Perb[i];
				for (j = 0; j < Nop; ++j)
					dfdx[j][i] = (f2[j] - f1[j]) / (2 * Perb[i]);
				dfdx[i][Nop] = f[i];
			}
			if (Gauss(Nop, dfdx, dx, MinGauss)) {
				ErrorNewt = 5;
				goto EndNewtonRlx;
			}
			SumValue = 0.;
			SumValueD = 0.;
			for (i = 0; i < Nop; ++i) {
				Value[i] -= dx[i];
				SumValue += fabs(Value[i]);
				SumValueD += fabs(dx[i]);
			}
			ErrorObf = (obj->*obf)(Nop, Value, f);
			if (ErrorObf) {
				ErrorNewt = 10 + ErrorObf;
				goto EndNewtonRlx;
			}
			remains2 = 0.;
			for (i = 0; i < Nop; ++i)
				remains2 += SQ(f[i]);
			if (SumValueD < (SumValue*MinPerb)) {
				ErrorNewt = 6;
				goto EndNewtonRlx;
			}
		}
		ErrorNewt = 0;
	EndNewtonRlx:
		*FinalRemains = remains2;
		//memory space free
		for (i = 0; i < Nop; ++i) {
			free(dfdx[i]);
		}
		free(dfdx);
		free(dx);
		free(f2);
		free(f1);
		free(f);
		return ErrorNewt;
	}
	//typedef int(*ObfFDM)(int,int,double*,double*);
	template <typename ClassName>
	int JacobiMtxFDM(
		int Nx
		, int Nf
		, double *TrimVal
		, ClassName *obj
		, int (ClassName::*obf)(int, int, double*, double*)
		, double **JacobiMtx
		, double *PerbFwd
		, double *PerbAft = NULL
	) {
		int i, j, Err = 0;
		double *fFwd, *fAft, *fMid;

		fFwd = (double*)malloc(sizeof(double)*Nf);
		fMid = (double*)malloc(sizeof(double)*Nf);
		fAft = (double*)malloc(sizeof(double)*Nf);

		if (Err = (obj->*obf)(Nx, Nf, TrimVal, fMid)) {
			goto EndJacobi;
		}

		for (i = 0; i < Nx; ++i) {
			TrimVal[i] += PerbFwd[i];
			if (Err = (obj->*obf)(Nx, Nf, TrimVal, fFwd)) {
				TrimVal[i] -= PerbFwd[i];
				goto EndJacobi;
			}
			TrimVal[i] -= PerbFwd[i];

			if (PerbAft == NULL) {
				for (j = 0; j < Nf; ++j) {
					JacobiMtx[j][i] = fFwd[j] / PerbFwd[i];
				}
			}
			else {
				TrimVal[i] += PerbAft[i];
				if (Err = (obj->*obf)(Nx, Nf, TrimVal, fAft)) {
					TrimVal[i] -= PerbAft[i];
					goto EndJacobi;
				}
				TrimVal[i] -= PerbAft[i];
				for (j = 0; j < Nf; ++j) {
					JacobiMtx[j][i] = (fFwd[j] - fAft[j]) / (PerbFwd[i] - PerbAft[i]);
				}
			}
		}


	EndJacobi:
		free(fAft);
		free(fMid);
		free(fFwd);

		return Err;
	}
}