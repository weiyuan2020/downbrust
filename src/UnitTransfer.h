/*
 * @Author: weiyuan 
 * @Date: 2021-09-02 14:14:19 
 * @Last Modified by:   weiyuan 
 * @Last Modified time: 2021-09-02 14:14:19 
 */
#pragma once

const double ft2m = 0.3048;							// ft -> m
const double in2m = 0.0254;							// in -> m
const double mslug2kg = 14.5939;					// slug -> kg
const double mlb2kg = 0.45359237;					// lbm -> kg
const double strans = 4.4482216152605;				// slug*ft -> kg*m
const double itrans = 1.3558179483314004;			// slug*ft2 -> kg*m2
const double lbins2_kgm2 = 0.1129848290276167;		// lb*in*s2 -> kg*m2
const double ppsi2pa = 6894.757;					// psi=lbf/inch2 -> pa = N/m2
const double patm2pa = 101325.0;					// atm -> pa
const double vkt2kmh = 1.852;						// knot -> km/h
const double vkt2fts = 1.688;						// knot -> ft/s
const double vkt2m = 0.5145024;						// knot -> m/s
const double d2r = 0.017453292;						// degree -> rad
const double hp2w = 735.499;						// hp(MHP)-> W
const double rpm_rps = 0.10471975512;				// rpm-> rad/s

const double K2R = 1.8;								// thermodynamic temperature -> Rankine Equals

													/*
													1 degree Fahrenheit ( f) = 32 + 1 degree Celsius x 1.8

													Kelvin (k) = 273.15 + degree Celsius

													leech ( Re) = degree Celsius  1.25

													Lambert ( R) = (degree Celsius + 273.15) x 1.8 degree

													Celsius ( C) = (degree Fahrenheit - 32)  1.8
													*/
const double flb2N = 4.4482216152605005;			// lbf -> N