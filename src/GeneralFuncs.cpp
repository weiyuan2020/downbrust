/*
 * @Author: L.F.Wang
 * @Date: 2021-04-26 09:46:25
 * @Last Modified by:   weiyuan
 * @Last Modified time: 2021-08-01 21:04:25
 */

#include "GeneralFuncs.h"
namespace MATH {
    int toBinaryNum(int num, const int& len, char* bi) {
        if (pow(2.0, len) < num)
            return 1;
        int i = 0;
        while (num != 0 && i < len) {
            bi[i++] = (char)num % 2;
            num /= 2;
        }
        while (i < len)
            bi[i++] = (char)0;
        return 0;
    }

    int Gauss(
        int Dim
        , double **A
        , double *x
        , double MinGauss/*=1.E-5*/
    ) {
        int i, j, k;
        double Val, temp;

        for (k = 0; k < Dim; ++k) {
            Val = A[k][k];
            i = k;
            while (fabs(Val) < MinGauss&&++i < Dim)
                Val = A[i][k];
            if (i == Dim) {
                std::cout << "ErrDim:" << k << std::endl;
                for (i = 0; i < Dim; ++i) {
                    for (j = 0; j < Dim + 1; ++j) {
                        std::cout << A[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                return 1;
            }
            if (i != k) {
                for (j = k; j < Dim + 1; ++j) {
                    temp = A[i][j];
                    A[i][j] = A[k][j];
                    A[k][j] = temp;
                }
            }
            for (j = k + 1; j < Dim + 1; ++j)
                A[k][j] /= Val;
            for (i = k + 1; i < Dim; ++i)
                for (j = k + 1; j < Dim + 1; ++j)
                    A[i][j] -= A[i][k] * A[k][j];
        }
        for (k = Dim - 1; k >= 0; --k) {
            for (i = Dim - 1; i > k; --i)
                A[k][Dim] -= A[i][Dim] * A[k][i];
        }
        for (i = 0; i < Dim; ++i)
            x[i] = A[i][Dim];
        return 0;
    };

    int GradDescent(
        int Nop
        , double *Value
        , double *Perb
        , int(*obf)(int, double*, double*)
        , double *FinalRemains
        , int MaxLoop        /*=40*/
        , int MaxRlx         /*=20*/
        , double MaxRemains	/*=1.E10*/
        , double MinRemains	/*=1.E-7*/
        , double MinPerb     /*=1.E-7*/
        , double MinGauss    /*=1.E-5*/
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

        double **dfdxT = (double**)malloc(sizeof(double*)*Nop);
        for (i = 0; i < Nop; ++i) {
            dfdxT[i] = (double*)malloc(sizeof(double)*(Nop + 1));
        }

        //Initialize
        remains1 = MaxRemains;
        looptimes = 0;
        if (ErrorObf = obf(Nop, Value, f)) {
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
                rlx *= 0.5;
                for (i = 0; i < Nop; ++i)
                    Value[i] += rlx * dx[i];
                if (ErrorObf = obf(Nop, Value, f)) {
                    ErrorNewt = 10 + ErrorObf;
                    goto EndNewtonRlx;
                }
                remains2 = 0.;
                for (i = 0; i < Nop; ++i)
                    remains2 += SQ(f[i]);
                //std::cout << "rlxtimes:" << rlxtimes << "    remains:" << remains2 << std::endl;
                if (remains2 < MinRemains) {
                    ErrorNewt = 0;//achieve the desired accuracy in rlx loop, return successful
                    goto EndNewtonRlx;
                }
            }
            if (remains2 <= remains1 * (1. - 0.25*rlx)) {
                remains1 = remains2;
            }
            for (i = 0; i < Nop; ++i) {
                Value[i] -= Perb[i];
                if (ErrorObf = obf(Nop, Value, f1)) {
                    ErrorNewt = 10 + ErrorObf;
                    goto EndNewtonRlx;
                }
                Value[i] += 2 * Perb[i];
                if (ErrorObf = obf(Nop, Value, f2)) {
                    ErrorNewt = 10 + ErrorObf;
                    goto EndNewtonRlx;
                }
                Value[i] -= Perb[i];
                for (j = 0; j < Nop; ++j)
                    dfdx[j][i] = (f2[j] - f1[j]) / (2 * Perb[i]);
                dfdx[i][Nop] = f[i];
            }

            for (i = 0; i < Nop; ++i) {
                for (j = 0; j < Nop + 1; ++j) {
                    dfdxT[i][j] = dfdx[i][j];
                }
            }

            if (Gauss(Nop, dfdx, dx, MinGauss)) {
                ErrorNewt = 5;
                std::cout << "dfdxT:" << std::endl;
                for (i = 0; i < Nop; ++i) {
                    for (j = 0; j < Nop + 1; ++j) {
                        std::cout << dfdxT[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "dfdx:" << std::endl;
                for (i = 0; i < Nop; ++i) {
                    for (j = 0; j < Nop + 1; ++j) {
                        std::cout << dfdx[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                goto EndNewtonRlx;
            }
            SumValue = 0.;
            SumValueD = 0.;
            for (i = 0; i < Nop; ++i) {
                Value[i] -= dx[i];
                SumValue += fabs(Value[i]);
                SumValueD += fabs(dx[i]);
            }
            if (ErrorObf = obf(Nop, Value, f)) {
                ErrorNewt = 10 + ErrorObf;
                goto EndNewtonRlx;
            }
            //remains1 = remains2;
            remains2 = 0.;
            for (i = 0; i < Nop; ++i)
                remains2 += SQ(f[i]);
            //std::cout << "looptimes:" << looptimes << "    remains:" << remains2 << std::endl;
    //         std::cout << "dfdx:" << std::endl;
    //         for ( i = 0; i < Nop; ++i ) {
    //             for ( j = 0; j < Nop + 1; ++j ) {
    //                 std::cout << dfdxT[i][j] << " ";
    //             }
    //             std::cout << std::endl;
    //         }
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
            free(dfdxT[i]);
        }
        free(dfdxT);
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

    int JacobiMtxFDM(
        int Nx
        , int Nf
        , double *TrimVal
        , int(*obf)(int, int, double*, double*)
        , double **JacobiMtx
        , double *PerbFwd
        , double *PerbAft/* = NULL*/
    ) {
        int i, j, Err = 0;
        double *fFwd, *fAft, *fMid;

        fFwd = (double*)malloc(sizeof(double)*Nf);
        fMid = (double*)malloc(sizeof(double)*Nf);
        fAft = (double*)malloc(sizeof(double)*Nf);

        if (Err = obf(Nx, Nf, TrimVal, fMid)) {
            goto EndJacobi;
        }

        for (i = 0; i < Nx; ++i) {
            TrimVal[i] += PerbFwd[i];
            if (Err = obf(Nx, Nf, TrimVal, fFwd)) {
                TrimVal[i] -= PerbFwd[i];
                goto EndJacobi;
            }
            TrimVal[i] -= PerbFwd[i];

            if (PerbAft == NULL) {
                for (j = 0; j < Nf; ++j) {
                    JacobiMtx[j][i] = (fFwd[j] - fMid[j]) / PerbFwd[i];
                }
            }
            else {
                TrimVal[i] += PerbAft[i];
                if (Err = obf(Nx, Nf, TrimVal, fAft)) {
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