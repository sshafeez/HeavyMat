// Copyright (C) 2009 foam
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

// from http://www.public.iastate.edu/~dicook/JSS/paper/paper.html

/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  this code and associated documentation,  provided       *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *     Copyright (c) held by Dianne Cook                    *
 *  All Rights Reserved.                                    *
 *                                                          *
 *  Questions and comments are welcome, and I request       *
 *  that you share any modifications with me.               *
 *                                                          *
 *                Dianne Cook                               *
 *             dicook@iastate.edu                           *
 *                                                          *
 ************************************************************/

/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <list>
#include "SVD.h"

using namespace std;

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x,y) ((x)>(y)?(x):(y))

///////////////////////////////////////////////////////

void normalizeGrid(Matrix<float> &grid) {
    for (int i = 0; i < grid.GetRows(); ++i) {
        for (int j = 0; j < grid.GetCols(); ++j) {
            if (grid[i][j] < 0) {
                grid[i][j] = 0;
            }
            if (grid[i][j] > 255) {
                grid[i][j] = 255;
            }
        }
    }
}

int calculateMaxError(Matrix<float> &first, Matrix<float> &second) {
    float* data1 = first.GetRawData();
    float* data2 = second.GetRawData();
    int max = 0;
    for (int i = 0; i < first.GetRows() * first.GetCols(); ++i) {
        if (abs(data1[i] - data2[i]) > max) {
            max = abs(data1[i] - data2[i]);
            //cout << data1[i] << " " << data2[i] << " " << abs(data1[i] - data2[i]) << "\n";
        }
        //cout << data1[i] << " " << data2[i] << "\n";
    }
    return max;
}

int calculateCumulativeError(Matrix<float> &first, Matrix<float> &second) {
    float* data1 = first.GetRawData();
    float* data2 = second.GetRawData();
    int sum = 0;
    for (int i = 0; i < first.GetRows() * first.GetCols(); ++i) {
        sum += abs(data1[i] - data2[i]);
    }
    return sum;
}

int multRatio(int sizeL, int colsR, double r) {
    int num = r*sizeL*sizeL*colsR;
    double denom = 2*sizeL*colsR+sizeL;
    cout << "Multiplication ratio " << r << endl;
    return sizeL - int(num/denom);
}

// returns # vecs that should
int numMultSavings(int rows, int cols, int savings) {
    int naive = rows*rows*cols;
    int svd = naive - savings;
    int denom = rows+2*rows*cols;
    double vec = (double)svd/denom;
    cout << vec << " " << vec-(int)vec << endl;
    if ((vec-(int)vec)>0.5) vec = (int)vec+1;
    cout << vec << endl;
    return (rows - int(vec));
}

double calcCompressionRatio(int rows, int vecs) {
    int numEntries = rows * rows;
    int numVecs = rows - vecs;
    int total = numVecs*2*rows + numVecs;
    double ratio = double(total) / numEntries;
    return 1-ratio;
}

int calcNumMults(int rows, int vecs, int cols) {
    return rows*vecs*(2*cols+1);
}

int calcNumAdds(int rows, int vecs, int cols) {
    return ((rows-1)*cols*vecs + (vecs-1)*rows*cols);
}

void eliminateBasis(Matrix<float> &Sg, float error, int cols) {
    int vecs = 0;
    for (int i = Sg.GetRows()-1; i >= 0; --i) {
        if (Sg[i][i] < error) {
            Sg[i][i] = 0;
            ++vecs;
        }
        else {
            cout << "Eliminated " << vecs << " vectors\n";
            cout << "Achieved compression ratio of " << calcCompressionRatio(Sg.GetRows(), vecs) << endl;
            int mults = calcNumMults(Sg.GetRows(), Sg.GetRows()-vecs, cols);
            int adds = calcNumAdds(Sg.GetRows(), Sg.GetRows()-vecs, cols);
            int naive = Sg.GetRows()*Sg.GetRows()*cols;
            cout << " # Mults: " << mults << endl;
            cout << " # Adds: " << adds << endl;
            cout << "Actual # mults saved: " << naive - mults << endl;
            break;
        }
    }
}

// # vecs to eliminate
void eliminateBasis(Matrix<float> &Sg, int vecs, int cols) {
    for (int i = 0; i < vecs; ++i) {
        Sg[Sg.GetRows()-1-i][Sg.GetRows()-1-i] = 0;
    }
    cout << "Eliminated " << vecs << " vectors\n";
    cout << "Achieved compression ratio of " << calcCompressionRatio(Sg.GetRows(), vecs) << endl;
    int mults = calcNumMults(Sg.GetRows(), Sg.GetRows()-vecs, cols);
    int adds = calcNumAdds(Sg.GetRows(), Sg.GetRows()-vecs, cols);
    int naive = Sg.GetRows()*Sg.GetRows()*cols;
    cout << " # Mults: " << mults << endl;
    cout << " # Adds: " << adds << endl;
    cout << "Actual # mults saved: " << naive - mults << endl;
}

void eliminateBasis(Matrix<float> &Sg, double ratio, int cols) {
    ratio = 1-ratio;
    int numEntries = Sg.GetRows() * Sg.GetRows();
    double total = ratio * numEntries;
    int numVecs = total / (2*Sg.GetRows()+1);
    int vecs = Sg.GetRows() - numVecs;
    for (int i = 0; i < vecs; ++i) {
        Sg[Sg.GetRows()-1-i][Sg.GetRows()-1-i] = 0;
    }
    cout << "Eliminated " << vecs << " vectors\n";
    cout << "Achieved compression ratio of " << calcCompressionRatio(Sg.GetRows(), vecs) << endl;
    cout << " # Mults: " << calcNumMults(Sg.GetRows(), Sg.GetRows()-vecs, cols) << endl;
    cout << " # Adds: " << calcNumAdds(Sg.GetRows(), Sg.GetRows()-vecs, cols) << endl;
}

Matrix<float> formMatrix(Vector<float> Sg) {
    float* data = (float*) malloc(sizeof(float)*Sg.Size()*Sg.Size());
    for (int i = 0; i < Sg.Size(); ++i) {
        for (int j = 0; j < Sg.Size(); ++j) {
            if (i == j) {
                data[i*Sg.Size()+j] = Sg[i];
            }
            else {
                data[i*Sg.Size()+j] = 0;
            }
        }
    }
    Matrix<float> S(Sg.Size(), Sg.Size(), data);
    //S.Print();
    return S;
}

Vector<float> SVD(Matrix<float> &a, Matrix<float> &b, Matrix<float> &c)
{
    if (a.IsInf()) { cout << "true\n"; }
    Vector<float> w(a.GetRows());
	Matrix<float> v(a.GetRows(),a.GetCols());
	dsvd(a, a.GetRows(), a.GetCols(), w.GetRawData(), v);
    Vector<float> w1 = w;
	a.SortCols(w);
    v.SortCols(w1);
    b = formMatrix(w);
    c = v;
	return w;
}

///////////////////////////////////////////////////////

static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}

int dsvd(Matrix<float> &a, int m, int n, float *w, Matrix<float> &v)
{
    if (v.IsInf()) { cout << "isinf1\n"; }
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (float)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (float)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (float)((double)a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (float)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (float)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (float)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
    //if (v.IsInf()) { cout << "isinf2\n"; }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++) {
                    /*if (((double)a[i][l] == 0) || (isnan((double)a[i][l]))) { cout << "error1\n"; }*/
                    v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
                    if (isnan(v[j][i])) { v[j][i] = 0; }
                    /* double division to avoid underflow */
                }
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) {
                        v[k][j] += (float)(s * (double)v[k][i]);
                        if (isnan(v[k][j])) { v[k][j] = 0; }
                    }
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
 //   if (v.IsInf()) { cout << "isinf2\n"; }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (float)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }
    //if (v.IsInf()) { cout << "isinf1\n"; }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                //if (v.IsInf()) { cout << "isinf4\n"; }
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    if (isnan(v[jj][j])) { v[jj][j] = 0; }
                    v[jj][i] = (float)(z * c - x * s);
                    if (isnan(v[jj][i])) { v[jj][i] = 0; }
                    //if (v.IsInf()) { cout << "isinf5\n"; }
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
   // if (a.IsInf()) { cout << "isinf2\n"; }
    //if (v.IsInf()) { cout << "isinf3\n"; }
    free((void*) rv1);
    return(1);
}
