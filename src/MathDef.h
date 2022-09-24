/*
 * @Author: L.F.Wang;
 * @Date: 2021-08-01 20:59:07 
 * @Last Modified by: weiyuan
 * @Last Modified time: 2021-08-01 20:59:52
 * Description: math constant and inline functions
 */

#pragma once
namespace MATH
{
	inline double sgn(double x) 
	{ 
		return ((x)>0. ? 1. : -1.);
	}
	inline double SQ(double x)
	{
		return ((x)*(x));
	}
	const double PI = 3.141592653589793115997963468544185161590576171875;	
	const double r2d = 57.29577951308232286464772187173366546630859375;		
	const double GRAVITY = 9.8;												
	const double MINRAD = 1.0E-5;											
	const int CharArrLen_Name = 256;										
}
