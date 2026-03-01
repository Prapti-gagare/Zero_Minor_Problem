// minor_relation_file.cpp
// Build: g++ -O2 minor_relation_file.cpp -lntl -lgmp -o minor_relation_file
// Requires NTL + GMP
#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <numeric>
#include <chrono>

using namespace std;
using namespace NTL;

// Utility: print vector<int> indices (1-based)
string idxList(const vector<int>& idx) {
    string s = "{";
    for (size_t i = 0; i < idx.size(); ++i) {
        s += to_string(idx[i] + 1);
        if (i + 1 < idx.size()) s += ",";
    }
    s += "}";
    return s;
}

// Determinant using Gaussian elimination
ZZ_p detZZp(Mat<ZZ_p> M) {
    long n = M.NumRows();
    if (n == 0) return ZZ_p(1);
    ZZ_p det = ZZ_p(1);
    for (long i = 0; i < n; i++) {
        long piv = i;
        while (piv < n && IsZero(M[piv][i])) ++piv;
        if (piv == n) return ZZ_p(0);
        if (piv != i) {
            swap(M[i], M[piv]);
            det = -det;
        }
        ZZ_p aii = M[i][i];
        ZZ_p aii_inv = inv(aii);
        det *= aii;
        for (long j = i; j < n; j++) M[i][j] *= aii_inv;
        for (long r = i + 1; r < n; r++) {
            ZZ_p f = M[r][i];
            if (IsZero(f)) continue;
            for (long c = i; c < n; c++) M[r][c] -= f * M[i][c];
        }
    }
    return det;
}

// Extract submatrix
Mat<ZZ_p> submatrix(const Mat<ZZ_p>& M, const vector<int>& rows, const vector<int>& cols) {
    long k = rows.size();
    Mat<ZZ_p> S;
    S.SetDims(k, k);
    for (long i = 0; i < k; i++)
        for (long j = 0; j < k; j++)
            S[i][j] = M[rows[i]][cols[j]];
    return S;
}

// Matrix inversion (Gauss-Jordan)
bool invertZZp(const Mat<ZZ_p>& A, Mat<ZZ_p>& Inv) {
    long n = A.NumRows();
    Mat<ZZ_p> B;
    B.SetDims(n, 2 * n);
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < n; j++) B[i][j] = A[i][j];
        for (long j = 0; j < n; j++) B[i][n + j] = (i == j ? ZZ_p(1) : ZZ_p(0));
    }
    for (long i = 0; i < n; i++) {
        long piv = i;
        while (piv < n && IsZero(B[piv][i])) ++piv;
        if (piv == n) return false;
        if (piv != i) swap(B[i], B[piv]);
        ZZ_p aii = B[i][i];
        ZZ_p inv_aii = inv(aii);
        for (long j = 0; j < 2 * n; j++) B[i][j] *= inv_aii;
        for (long r = 0; r < n; r++) {
            if (r == i) continue;
            ZZ_p factor = B[r][i];
            if (IsZero(factor)) continue;
            for (long c = 0; c < 2 * n; c++) B[r][c] -= factor * B[i][c];
        }
    }
    Inv.SetDims(n, n);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++)
            Inv[i][j] = B[i][n + j];
    return true;
}

// Generate pattern matrix
bool generatePatternMatrix(int n, const ZZ& Pzz, Mat<ZZ_p>& A,
                           int max_global_restarts = 200, int max_row_attempts = 500) {
    long p_long;
    try {
        p_long = conv<long>(Pzz);
    } catch (...) {
        cerr << "Prime too large.\n";
        return false;
    }
    if (p_long - 1 < n * n) {
        cerr << "Prime too small.\n";
        return false;
    }

    vector<int> pool;
    for (int x = 1; x <= (int)(p_long - 1); ++x) pool.push_back(x);

    mt19937_64 rng(chrono::high_resolution_clock::now().time_since_epoch().count());

    for (int global_try = 0; global_try < max_global_restarts; ++global_try) {
        vector<int> remaining = pool;
        shuffle(remaining.begin(), remaining.end(), rng);
        A.SetDims(n, n);
        bool fail_global = false;

        for (int r = 0; r < n && !fail_global; ++r) {
            bool row_done = false;
            for (int attempt = 0; attempt < max_row_attempts && !row_done; ++attempt) {
                if ((int)remaining.size() < n) { fail_global = true; break; }
                shuffle(remaining.begin(), remaining.end(), rng);
                vector<int> chosen(remaining.begin(), remaining.begin() + (n - 1));
                long sum = accumulate(chosen.begin(), chosen.end(), 0L);
                long last = (p_long - 1) - sum;
                if (last <= 0 || last >= p_long) continue;
                if (find(chosen.begin(), chosen.end(), (int)last) != chosen.end()) continue;
                if (find(remaining.begin(), remaining.end(), (int)last) == remaining.end()) continue;
                chosen.push_back((int)last);
                unordered_set<int> remset(chosen.begin(), chosen.end());
                vector<int> newrem;
                for (int v : remaining)
                    if (!remset.count(v)) newrem.push_back(v);
                remaining.swap(newrem);
                for (int c = 0; c < n; c++) A[r][c] = ZZ_p(chosen[c]);
                row_done = true;
            }
            if (!row_done) fail_global = true;
        }
        if (!fail_global) return true;
    }
    return false;
}

// Generate combinations
void combinations(int n, int k, vector<vector<int>>& out) {
    vector<int> bitmask(k, 1);
    bitmask.resize(n, 0);
    do {
        vector<int> comb;
        for (int i = 0; i < n; i++) if (bitmask[i]) comb.push_back(i);
        out.push_back(comb);
    } while (prev_permutation(bitmask.begin(), bitmask.end()));
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ofstream fout("OUTPUT.txt");
    if (!fout) { cerr << "Failed to open output2.txt\n"; return 1; }

    int n;
    long P_long;
    cout << "Enter n (matrix size n x n): ";
    if (!(cin >> n)) return 0;
    cout << "Enter prime P (must satisfy P-1 >= n*n): ";
    if (!(cin >> P_long)) return 0;

    ZZ Pzz = to_ZZ(P_long);
    ZZ_p::init(Pzz);

    fout << "===== Pattern Matrix Generation & Minor Relation Verification =====\n";
    fout << "Matrix size n = " << n << "\nPrime P = " << P_long << "\n\n";

    Mat<ZZ_p> A;
    cout << "Generating pattern matrix A...\n";
    bool ok = generatePatternMatrix(n, Pzz, A);
    if (!ok) { fout << "Failed to generate pattern matrix.\n"; return 1; }

    fout << "Pattern Matrix A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            fout << conv<long>(rep(A[i][j])) << "\t";
        fout << "\n";
    }

    fout << "\nRow sums (should equal P-1):\n";
    for (int i = 0; i < n; i++) {
        long sum = 0;
        for (int j = 0; j < n; j++) sum += conv<long>(rep(A[i][j]));
        fout << "Row " << i + 1 << ": " << sum << "\n";
    }

    ZZ_p detA = detZZp(A);
    fout << "\nDet(A) mod P = " << conv<long>(rep(detA)) << "\n";

    Mat<ZZ_p> Ainv;
    bool inv_ok = invertZZp(A, Ainv);
    if (!inv_ok) { fout << "Matrix not invertible.\n"; return 1; }

    fout << "\nInverse Matrix A^-1:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            fout << conv<long>(rep(Ainv[i][j])) << "\t";
        fout << "\n";
    }

    for (int k = 1; k <= n; k++) {
        fout << "\n========== k = " << k << " minors ==========\n";
        vector<vector<int>> rowCombs; combinations(n, k, rowCombs);
        vector<vector<int>> colCombs; combinations(n, k, colCombs);

        int total = 0, matches = 0;
        for (const auto& I : rowCombs) {
            for (const auto& J : colCombs) {
                total++;
                Mat<ZZ_p> subA = submatrix(A, I, J);
                Mat<ZZ_p> subAinv = submatrix(Ainv, I, J);
                ZZ_p det_subAinv = detZZp(subAinv);

                vector<int> Ic, Jc;
                for (int t = 0; t < n; t++) {
                    if (find(I.begin(), I.end(), t) == I.end()) Ic.push_back(t);
                    if (find(J.begin(), J.end(), t) == J.end()) Jc.push_back(t);
                }
                Mat<ZZ_p> comp = submatrix(A, Jc, Ic);
                ZZ_p det_comp = detZZp(comp);

                long sumI = 0, sumJ = 0;
                for (int v : I) sumI += (v + 1);
                for (int v : J) sumJ += (v + 1);
                ZZ_p sign = ((sumI + sumJ) % 2 == 0) ? ZZ_p(1) : ZZ_p(-1);
                ZZ_p rhs = sign * det_comp * inv(detA);

                fout << "I=" << idxList(I)
                     << "  J=" << idxList(J)
                     << "  det(A^-1)_{I,J}=" << conv<long>(rep(det_subAinv))
                     << "  RHS=" << conv<long>(rep(rhs))
                     << (det_subAinv == rhs ? " OK" : " Mismatch")
                     << "\n";
                if (det_subAinv == rhs) matches++;
            }
        }
        fout << "k=" << k << ": verified " << matches << " / " << total << " minors.\n";
    }

    fout << "\nComputation complete. Results saved successfully.\n";
    fout.close();

    cout << "Computation finished. Results saved to 'output2-1.txt'.\n";
    return 0;
}
// minor_relation_file.cpp
// Build: g++ -O2 minor_relation_file.cpp -lntl -lgmp -o minor_relation_file
// Requires NTL + GMP

