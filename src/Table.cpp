/*
 * @Author:  L.F.Wang
 * @Date: 2021-08-01 21:05:07
 * @Last Modified by:   weiyuan
 * @Last Modified time: 2021-08-01 21:05:07
 */
 /*
 Name: Table.cpp
 Copyright:
 Author: L.F.Wang
 Date: 16/04/21 10:05
 Description:

 */

#include "Table.h"
namespace MATH {
    using std::cout;

    Table::Table() {
        ifInit = false;
    }
    Table::~Table() {
        if (ifInit == true)
            RlsSpace();
    }

    double Table::GetData(int* loc) const {
        int loc1, i;
        if (data == NULL) {
            cout << "data is null" << std::endl;
            return 0.;
        }
        for (i = 0; i < Ndim; ++i) {
            if (loc[i] < 0 || loc[i] >= Ntitle[i]) {
                cout << "location is not in range" << std::endl;
                return 0.;
            }
        }
        loc1 = 0;
        for (i = 0; i < Ndim; ++i) {
            loc1 += loc[i] * dimWide[i];
        }

        return this->data[loc1];
    }

    void Table::RlsSpace() {
        delete[] data;
        for (int i = 0; i < Ndim; ++i) {
            if (title[i] != NULL)
                delete[] title[i];
        }
        delete[] title;
        delete[] dimWide;
        delete[] Ntitle;
    }

    int Table::Init(
        int _Ndim,
        const int *_Ntitle,
        const double **_title,
        const double *_data
    ) {
        int i, j;

        if (ifInit == true)
            RlsSpace();

        Ndim = _Ndim;
        Ntitle = new int[Ndim];
        dimWide = new int[Ndim];
        title = new double*[Ndim];

        for (i = 0; i < Ndim; ++i)
            Ntitle[i] = _Ntitle[i];
        int Ndata = 1;
        for (i = 0; i < Ndim; ++i) {
            title[i] = new double[Ntitle[i]];
            Ndata *= Ntitle[i];
        }
        data = new double[Ndata];

        for (i = 0; i < Ndim; ++i) {
            for (j = 0; j < Ntitle[i]; ++j) {
                title[i][j] = _title[i][j];
            }
        }
        for (i = 0; i < Ndata; ++i)
            data[i] = _data[i];

        for (i = 0; i < Ndim; ++i) {
            dimWide[i] = 1;
            for (j = 0; j < i; ++j)
                dimWide[i] *= Ntitle[j];
        }

        ifInit = true;
        return 0;
    }

    int Table::Init(
        int _Ndim,
        int *_Ntitle,
        double **_title,
        double *_data
    ) {
        int i, j;

        if (ifInit == true)
            RlsSpace();

        Ndim = _Ndim;
        Ntitle = new int[Ndim];
        dimWide = new int[Ndim];
        title = new double*[Ndim];

        for (i = 0; i < Ndim; ++i)
            Ntitle[i] = _Ntitle[i];
        int Ndata = 1;
        for (i = 0; i < Ndim; ++i) {
            title[i] = new double[Ntitle[i]];
            Ndata *= Ntitle[i];
        }
        data = new double[Ndata];

        for (i = 0; i < Ndim; ++i) {
            for (j = 0; j < Ntitle[i]; ++j) {
                title[i][j] = _title[i][j];
            }
        }
        for (i = 0; i < Ndata; ++i)
            data[i] = _data[i];

        for (i = 0; i < Ndim; ++i) {
            dimWide[i] = 1;
            for (j = 0; j < i; ++j)
                dimWide[i] *= Ntitle[j];
        }

        ifInit = true;
        return 0;
    }

    int Table::Init(const char *outFile) {
        using namespace std;
        ifstream ifs(outFile);
        if (!ifs.is_open()) {
            cout << "OPEN FILE FAIL!" << endl;
            return 1;
        }

        int i, j;

        if (ifInit == true)
            RlsSpace();

        ifs >> Ndim;
        Ntitle = new int[Ndim];
        dimWide = new int[Ndim];
        title = new double*[Ndim];

        for (i = 0; i < Ndim; ++i)
            ifs >> Ntitle[i];
        int Ndata = 1;
        for (i = 0; i < Ndim; ++i) {
            title[i] = new double[Ntitle[i]];
            Ndata *= Ntitle[i];
        }
        data = new double[Ndata];

        for (i = 0; i < Ndim; ++i) {
            for (j = 0; j < Ntitle[i]; ++j) {
                ifs >> title[i][j];
            }
        }

        for (i = 0; i < Ndata; ++i)
            ifs >> data[i];

        ifs.close();

        for (i = 0; i < Ndim; ++i) {
            dimWide[i] = 1;
            for (j = 0; j < i; ++j)
                dimWide[i] *= Ntitle[j];
        }

        ifInit = true;
        return 0;
    }

    double Table::LookUp(const double *pos)const {
        int **loc;
        int i, j, k, vertexNum;
        char *bi;
        double *vertexData, ratio;
        double val;
        // 	for(i=0;i<Ndim;++i){
        // 		if(pos[i]<title[i][0] || pos[i]>title[i][Ntitle[i]-1]){
        // 			os<<"caution: Dim ("<<i+1<<") of the data point is not in table range"<<std::endl;
        // 			os<<"position is:";
        // 			for(j=0;j<Ndim;++j)
        // 				os<<pos[i]<<",";
        // 			os<<std::endl;
        // 		}
        // 	}
        vertexNum = (int)pow(2.0, Ndim);
        loc = new int*[vertexNum];
        for (i = 0; i < vertexNum; ++i) {
            loc[i] = new int[Ndim];
        }
        bi = new char[Ndim];
        vertexData = new double[vertexNum];
        for (i = 0; i < Ndim; ++i) {
            loc[0][i] = TmpBinarySearch(0, Ntitle[i] - 1, title[i], pos[i]);
        }
        for (i = 0; i < vertexNum; ++i) {
            toBinaryNum(i, Ndim, bi);
            for (j = 0; j < Ndim; ++j)
                loc[i][j] = loc[0][j] + (int)bi[j];
        }
        for (i = 0; i < vertexNum; ++i) {
            vertexData[i] = GetData(loc[i]);
        }
        for (i = Ndim - 1; i >= 0; --i) {
            ratio = (pos[i] - title[i][loc[0][i]]) / (title[i][loc[vertexNum - 1][i]] - title[i][loc[0][i]]);
            k = (int)pow(2.0, i);
            for (j = 0; j < k; ++j)
                vertexData[j] += ratio*(vertexData[j + k] - vertexData[j]);
        }
        val = vertexData[0];
        delete[] vertexData;
        delete[] bi;
        for (i = 0; i < vertexNum; ++i) {
            delete[] loc[i];
        }
        delete[] loc;

        return val;
    }
}