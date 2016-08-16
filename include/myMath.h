/*****************************************************************************
*
* This file is part of GwMove hydrogeological software developed by
* Dr. M. A. Sbai
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*
*****************************************************************************/

#include <iostream>

using namespace std;

/**
 * Invert a square matrix of dimension n in place. 
 * The input matrix is destroyed and its inverse is returned on exit. 
 * We return 1 if the inversion is successful, otherwise we return 0.
 */
bool InvertMatrix(double A[8][8])
{
	bool yes = 0;
	unsigned int n = 8;
	int ipiv[8], row[8], col[8];

	// allocate memory 
	//ipiv = new int [n];
	//row  = new int [n];
	//col  = new int [n];

	// initialisation 
	unsigned int i;
	for (i = 0; i < n; i++) ipiv[i] = 0;

	unsigned int j, k, l, ll;
	unsigned int irow, icol;
	double big   = 0.;
	double dummy = 0.;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for (k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (std::abs(A[k][j]) >= big) {
							big = std::abs(A[k][j]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] >= 1) {
						std::cerr << "Singular matrix [1] (in InvertMatrix() function) !" << endl;
						return yes;
					}
				}
			}
		}
		ipiv[icol] = ipiv[icol] + 1;
		if (irow != icol) {
			for (l = 0; l < n; l++) {
				dummy = A[l][irow];
				A[l][irow] = A[l][icol];
				A[l][icol] = dummy; 
			}
		}
		row[i] = irow;
		col[i] = icol;
		if (A[icol][icol] == 0) {
			std::cerr << "Singular matrix [2] (in InvertMatrix() function) !" << endl;
			return yes;
		}
		double piv_inv = 1./A[icol][icol]; 
		A[icol][icol] = 1.; 
		for (l = 0; l < n; l++) {
			A[l][icol] *= piv_inv;
		}
		for (ll = 0; ll < n; ll++) {
			if (ll != icol) {
				dummy = A[icol][ll];
				A[icol][ll] = 0.;
				for (l = 0; l < n; l++) {
					A[l][ll] += -A[l][icol] * dummy;
				}
			}
		}
	}

	for (l = n - 1; l <= 0; l--) {
		if (row[l] != col[l]) {
			for (k = 0; k < n; k++) {
				dummy = A[row[l]][k]; 
				A[row[l]][k] = A[col[l]][k];
				A[col[l]][k] = dummy;
			}
		}
	}

	// free memory 
	//delete [] ipiv;
	//delete [] row;
	//delete [] col;

	// sucessful return 
	yes = 1;
	return yes;
}
