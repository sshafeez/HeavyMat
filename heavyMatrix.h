/*
    Defines a class heavy_matrix which can cache linear dependencies.
*/
#include <vector>
#include <utility>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>
#include <math.h>
#include "matrix.h"

//Least Squares Regression
matrix LSR(matrix &A, matrix &b)
{
	matrix At = A;
	At.transpose();
	matrix AtA(0, 0);
	multiply(At, A, AtA);
	AtA.invert();
	matrix AtAAt(0, 0);
	multiply(AtA, At, AtAAt);
	matrix AtAAtb(0, 0);
	multiply(AtAAt, b, AtAAtb);
	return AtAAtb;
}

//Compute angle between A[v1] and vector b
double angle(const matrix &A, int v1, const matrix &b)
{
	//dot product
	double dot = 0;
	for (int i = 0; i < A.grid[v1].size(); ++i)
	{
		dot += A.grid[v1][i] * b.grid[i][0];
	}

	//vector magnitude
	double magv1 = 0;
	double magv2 = 0;
	for (int i = 0; i < A.grid[v1].size(); ++i)
	{
		magv1 += A.grid[v1][i] * A.grid[v1][i];
		magv2 += b.grid[i][0] * b.grid[i][0];
	}
	magv1 = sqrt(magv1);
	magv2 = sqrt(magv2);

	return abs(acos(dot / (magv1 * magv2)));
}

class heavy_matrix : public matrix
{
public:
	//Store the dependencies between vectors
	struct dep
	{
		int index;
		double scalar;
	};
	unordered_map<int, vector<dep>> rowDeps;
	unordered_map<int, vector<dep>> colDeps;

	//construct from scratch or std::vector<std::vector<>>
	heavy_matrix(int rows, int cols, bool initRandom = false) : matrix(rows, cols, initRandom) {}
	heavy_matrix(const heavy_matrix &other) : matrix(other) {}
	heavy_matrix(const vector<vector<double>> &data) : matrix(data) {}

	//find and cache linear dependencies
	void cache(double error = epsilon)
	{
		rowDeps.clear();
		colDeps.clear();

		//check rows then columns
		for (int dimension = 0; dimension < 2; ++dimension)
		{
			for (int i = 1; i < size().first; ++i)
			{
				//start with empty basis
				matrix basis(0, 0);
				vector<int> rowsInBasis;

				//vector to be reconstructed from basis
				matrix b(0, 0);
				b.append(grid[i]);
				//reconstruction error vector remaining
				matrix rem(0, 0);
				rem.append(grid[i]);
				//store the reconstruction dependencies
				vector<dep> weights;

				//compress until error threshold met
				int iters = 0;
				double cost = error + 1;
				while (cost > error && iters < i)
				{
					//find closest vector to reconstruction error
					int bestIndex = i - 1;
					double bestAngle = 361;
					for (int k = i - 1; k >= 0; --k)
					{
						//avoid duplicates (to maintain matrix invertibility)
						if (dimension == 0 && rowDeps.count(k) != 0)
							continue;
						if (dimension == 1 && colDeps.count(k) != 0)
							continue;
						if (std::find(rowsInBasis.begin(), rowsInBasis.end(), k) != rowsInBasis.end())
							continue;

						if (angle(*this, k, rem) < bestAngle)
						{
							bestAngle = angle(*this, k, rem);
							bestIndex = k;
						}
					}

					//add the chosen vector to the basis built so far
					basis.append(grid[bestIndex]);
					rowsInBasis.push_back(bestIndex);

					//get least squares and store the mapping
					matrix basisMapping = LSR(basis, b);
					weights.clear();
					for (int k = 0; k < basisMapping.grid.size(); ++k)
					{
						if (abs(basisMapping.grid[k][0]) > epsilon)
							weights.push_back({rowsInBasis[k], basisMapping.grid[k][0]});
					}

					//compute reconstruction/projection from the basis and mapping
					matrix proj(0, 0);
					multiply(basis, basisMapping, proj);
					proj.saturate();

					//recompute error and remaining error vector (projection error)
					cost = 0;
					for (int k = 0; k < b.size().first; ++k)
					{
						cost = max(cost, abs(b.grid[k][0] - proj.grid[k][0]));
						rem.grid[k][0] = b.grid[k][0] - proj.grid[k][0];
					}

					iters++;
				}

				//store the mapping/weights if it will reduce floating point operations
				if (cost < error && weights.size() < b.size().first)
				{
					if (dimension == 0)
						rowDeps[i] = weights;
					else
						colDeps[i] = weights;
				}
			}
			transpose();
		}
	}

	//compress image based on cached linear dependencies
	void writeback()
	{
		//decide dimension which will minimize floating point computations
		int rowSavings = 0;
		for (int i = 0; i < grid.size(); ++i)
		{
			if (rowDeps.count(i))
			{
				rowSavings += grid[0].size() - rowDeps[i].size();
			}
		}
		int colSavings = 0;
		for (int i = 0; i < grid[0].size(); ++i)
		{
			if (colDeps.count(i))
			{
				colSavings += grid.size() - colDeps[i].size();
			}
		}
		cout << "rowSavings: " << rowSavings << " colSavings: " << colSavings << endl;

		//rewrite rows or columns in terms of cached linear dependencies
		if (rowSavings > colSavings)
		{
			vector<vector<double>> compressed(grid.size(), vector<double>(grid[0].size(), 0));
			for (int i = 0; i < grid.size(); ++i)
			{
				if (rowDeps.count(i))
				{
					for (int k = 0; k < rowDeps[i].size(); ++k)
					{
						for (int j = 0; j < grid[0].size(); ++j)
						{
							compressed[i][j] += rowDeps[i][k].scalar * grid[rowDeps[i][k].index][j];
						}
					}
				}
				else
				{
					compressed[i] = grid[i];
				}
			}
			grid = compressed;
		}
		else
		{
			transpose();
			vector<vector<double>> compressed(grid.size(), vector<double>(grid[0].size(), 0));
			for (int i = 0; i < grid.size(); ++i)
			{
				if (colDeps.count(i))
				{
					for (int k = 0; k < colDeps[i].size(); ++k)
					{
						for (int j = 0; j < grid[0].size(); ++j)
						{
							compressed[i][j] += colDeps[i][k].scalar * grid[colDeps[i][k].index][j];
						}
					}
				}
				else
				{
					compressed[i] = grid[i];
				}
			}
			grid = compressed;
			transpose();
		}
		saturate();
	}
};

//multiply two matrices using cached linear dependencies.
void fastMultiply(heavy_matrix &left, heavy_matrix &right, heavy_matrix &dest)
{
	/// (mxn) * (nxq) = (mxq)
	if (left.size().second != right.size().first)
		throw "inner dimension mismatch";
	int rows = left.size().first;
	int inner = left.size().second;
	int cols = right.size().second;

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			dest.at(i, j) = 0;
			//use cached dependency if it exists for these coordinates
			if (left.rowDeps.count(i) || right.colDeps.count(j))
			{
				vector<heavy_matrix::dep> &deps = left.rowDeps.count(i) ? left.rowDeps[i] : right.colDeps[j];
				for (heavy_matrix::dep &comb : deps)
				{
					floating_mults++;
					floating_adds++;
					if (left.rowDeps.count(i))
						dest.at(i, j) += comb.scalar * dest.at(comb.index, j);
					else
						dest.at(i, j) += comb.scalar * dest.at(i, comb.index);
				}
				floating_adds--;
			}
			//otherwise compute normally
			else
			{
				for (int k = 0; k < inner; ++k)
				{
					floating_mults++;
					floating_adds++;
					dest.at(i, j) += left.at(i, k) * right.at(k, j);
				}
				floating_adds--;
			}
		}
	}
}

//compute total difference between two matrices
double calcError(matrix &left, matrix &right)
{
	if (left.size() != right.size())
		throw "mismatch";
	double maxErr = 0, total = 0;
	for (int i = 0; i < left.grid.size(); ++i)
	{
		for (int j = 0; j < left.grid[0].size(); ++j)
		{
			total += abs(left.grid[i][j] - right.grid[i][j]);
			maxErr = max(maxErr, abs(left.grid[i][j] - right.grid[i][j]));
		}
	}
	cout << "Max Error: " << maxErr << endl;
	cout << "Total Error: " << total << endl;
	return total;
}
