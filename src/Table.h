/*
 * @Author: L.F.Wang
 * @Date: 2021-04-16 21:03:28
 * @Last Modified by:   weiyuan
 * @Last Modified time: 2021-08-01 21:03:28
 */

#pragma once
#include <iostream>
#include <fstream>
#include "GeneralFuncs.h"
namespace MATH
{

	class Table
	{
	public:

		Table();
		~Table();
		int Init(int _Ndim, const int *_Ntitle, const double **_title, const double *_data);
		int Init(int _Ndim, int *_Ntitle, double **_title, double *_data);
		int Init(const char *outFile);
		double LookUp(const double *pos)const;

	private:
		bool ifInit;
		int* dimWide;

		int Ndim;
		int *Ntitle;
		double **title;
		double *data;

		double GetData(int* loc)const;
		void RlsSpace();

	};
}