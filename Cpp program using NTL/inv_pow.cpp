#include <bits/stdc++.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

using namespace std;
using namespace NTL;

// ------------------ Determinant for small minors (mod p) ---------------------
long long determinant(const vector<vector<long long>>& M_in, long long p) {
    int n = (int)M_in.size();
    // work on a local copy and reduce mod p
    vector<vector<long long>> M(n, vector<long long>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) 
        {
            long long x = M_in[i][j] % p;
            if (x < 0) x += p;
            M[i][j] = x;
        }

    if (n == 1) return M[0][0] % p;
    if (n == 2)\
     {
        long long det = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % p;
        if (det < 0) det += p;
        return det;
    }

    long long det = 0;
    for (int col = 0; col < n; col++) 
    {
        vector<vector<long long>> sub(n - 1, vector<long long>(n - 1));
        for (int i = 1; i < n; i++) 
        {
            int cj = 0;
            for (int j = 0; j < n; j++) 
            {
                if (j == col) continue;
                sub[i - 1][cj++] = M[i][j];
            }
        }
        long long sign = (col % 2 == 0 ? 1 : -1);
        long long subdet = determinant(sub, p);
        long long term = ( ( (sign * M[0][col]) % p ) * subdet ) % p;
        det = (det + term) % p;
    }
    if (det < 0) det += p;
    return det;
}

bool isZeroMod(long long x, long long p) 
{
    x %= p;
    if (x < 0) x += p;
    return x == 0;
}

// ------------------ Fast matrix exponentiation for NTL ---------------------
mat_ZZ_p matPower(mat_ZZ_p base, long long exp) {
    long n = base.NumRows();
    mat_ZZ_p result;
    result.SetDims(n, n);

    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++)
            result[i][j] = (i == j ? ZZ_p(1) : ZZ_p(0));

    while (exp > 0) {
        if (exp & 1) result = result * base;
        base = base * base;
        exp >>= 1;
    }
    return result;
}

int main() 
{
    ofstream out("output_inversedemo.txt");
    if (!out.is_open())
     {
        cout << "ERROR: Could not open file output_inverse4.txt\n";
        return 0;
    }

    int n;
    long long p;
    cout << "Enter n: ";
    if (!(cin >> n)) return 0;
    cout << "Enter prime p: ";
    if (!(cin >> p)) return 0;

    if (n <= 0)
     {
        cout << "n must be positive\n";
        return 0;
    }
    if (p <= 1) 
    {
        cout << "p must be a prime > 1\n";
        return 0;
    }

    out << "Matrix Size n = " << n << "\n";
    out << "Prime p = " << p << "\n\n";

    ZZ_p::init(ZZ(p));

    vector<vector<long long>> A_ll(n, vector<long long>(n));
    mat_ZZ_p A;
    A.SetDims(n, n);

    cout << "Enter matrix elements (" << n * n << " numbers):\n";
    out << "Input Matrix:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long long v;
            cin >> v;
            A_ll[i][j] = v;
            A[i][j] = ZZ_p(v);
            out << v << " ";
        }
        out << "\n";
    }
    out << "\n";

    // ----------- Compute determinant & inverse using NTL ------------
    ZZ_p detA;
    mat_ZZ_p Ainv;
    bool invertible = true;

    try 
    {
        inv(detA, Ainv, A);
    } catch (...)
     {
        invertible = false;
    }

    out << "Determinant mod p = " << detA << "\n";

    if (!invertible || IsZero(detA)) 
    {
        out << "Matrix NOT invertible modulo " << p << "\n\n";
    } else
     {
        out << "\nInverse (A^-1 mod p):\n";
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++)
                out << Ainv[i][j] << " ";
            out << "\n";
        }
        out << "\n";
    }

    // ----------- Scan minors and write zero minors ----------
    out << "----------------------------------------------\n";
    out << "ZERO MINORS\n";
    out << "----------------------------------------------\n\n";

    bool found = false;

    // For each minor size k (1..n-1)
    for (int k = 1; k < n; ++k)
     {
        // create base masks in non-decreasing order: false...false,true...true
        vector<bool> baseR(n, false), baseC(n, false);
        fill(baseR.end() - k, baseR.end(), true); // important: put true at the end
        fill(baseC.end() - k, baseC.end(), true);

        vector<bool> selectR = baseR;
        // iterate over all row combinations
        do
         {
            // reset column mask for each row-selection
            vector<bool> selectC = baseC;
            do 
            {
                vector<int> rset, cset;
                for (int i = 0; i < n; ++i) if (selectR[i]) rset.push_back(i);
                for (int j = 0; j < n; ++j) if (selectC[j]) cset.push_back(j);

                // build k x k minor (use entries reduced mod p)
                vector<vector<long long>> M(k, vector<long long>(k));
                for (int i = 0; i < k; ++i)
                    for (int j = 0; j < k; ++j)
                     {
                        long long x = A_ll[rset[i]][cset[j]] % p;
                        if (x < 0) x += p;
                        M[i][j] = x;
                    }

                long long d = determinant(M, p);
                if (isZeroMod(d, p))
                 {
                    found = true;
                    out << "Zero minor of size " << k << "x" << k << "\n";
                    out << "Rows: ";
                    for (int x : rset) out << x << " ";
                    out << "\nCols: ";
                    for (int x : cset) out << x << " ";
                    out << "\nMatrix:\n";
                    for (auto &row : M) {
                        for (auto &x : row) out << x << " ";
                        out << "\n";
                    }
                    out << "Det = 0 (mod " << p << ")\n\n";
                }

            } 
            while (next_permutation(selectC.begin(), selectC.end()));
        } 
        while (next_permutation(selectR.begin(), selectR.end()));
    }

    if (!found)
        out << "No zero minors found.\n\n";

    // ------------ ODD POWER OPERATION ON INVERSE -------------------
    if (!IsZero(detA))
     {
        int m;
        cout << "Enter m (odd powers of inverse up to m): ";
        cin >> m;

        out << "----------------------------------------------\n";
        out << "ODD POWERS OF INVERSE MATRIX\n";
        out << "----------------------------------------------\n\n";

        for (int k = 1; k <= m; k += 2) {
            out << "(A^-1)^" << k << ":\n";
            mat_ZZ_p P = matPower(Ainv, k);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++)
                    out << P[i][j] << " ";
                out << "\n";
            }
            out << "\n";
        }
    }

    out.close();
    cout << "\nAll results saved to output_inversedemo.txt\n";
    return 0;
}
