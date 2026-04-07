#include <bits/stdc++.h>
using namespace std;

typedef vector<vector<int>> Matrix;

int MOD;

// ─── Matrix operations ────────────────────────────────────────────────────────

Matrix multiply(const Matrix &A, const Matrix &B, int n) {
    Matrix C(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) {
            if (A[i][k] == 0) continue;
            for (int j = 0; j < n; j++)
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % MOD;
        }
    return C;
}

Matrix identity(int n) {
    Matrix I(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) I[i][i] = 1;
    return I;
}

Matrix power(Matrix A, long long exp, int n) {
    Matrix result = identity(n);
    while (exp > 0) {
        if (exp & 1) result = multiply(result, A, n);
        A = multiply(A, A, n);
        exp >>= 1;
    }
    return result;
}

bool isIdentity(const Matrix &A, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if ((i == j && A[i][j] != 1) || (i != j && A[i][j] != 0))
                return false;
    return true;
}

// ─── Parser for zero_minor_matrices.txt ──────────────────────────────────────

bool parseNextMatrix(ifstream &fin, Matrix &mat, int &n, int &p) {
    string line;
    while (getline(fin, line)) {
        if (line.rfind("Matrix (mod", 0) == 0) {
            size_t lp = line.find('(');
            size_t rp = line.find(')');
            string inner = line.substr(lp + 1, rp - lp - 1);
            istringstream ss(inner);
            string kw; int pVal;
            ss >> kw >> pVal;
            p = pVal;

            mat.clear();
            while (getline(fin, line)) {
                if (line.find("---") != string::npos) break;
                istringstream rs(line);
                vector<int> row;
                int v;
                while (rs >> v) row.push_back(v);
                if (!row.empty()) mat.push_back(row);
            }
            if (mat.empty()) continue;
            n = (int)mat.size();
            return true;
        }
    }
    return false;
}

// ─── Print one permutation with SAFE/UNSAFE label ────────────────────────────

void writePermutation(ofstream &fout, const Matrix &mat, int n, int p,
                      int permNum, bool safe) {
    fout << "  Permutation " << permNum << " : "
         << (safe ? "SAFE" : "UNSAFE") << "\n";
    for (int i = 0; i < n; i++) {
        fout << "  ";
        for (int j = 0; j < n; j++)
            fout << mat[i][j] << " ";
        fout << "\n";
    }
    fout << "\n";
}

// ─── Generate all n! permutations (Johnson-Trotter), check & print each ───────

void processPermutations(const Matrix &mat, int n, int p,
                         ofstream &fout,
                         int &safeCount, int &unsafeCount, int &total,
                         int matrixNum) {

    int fact = 1;
    for (int i = 1; i <= n; i++) fact *= i;

    fout << "=============================================================\n";
    fout << "Source Matrix #" << matrixNum << "  (mod " << p << ")\n";
    fout << "Original rows:\n";
    for (int i = 0; i < n; i++) {
        fout << "  ";
        for (int j = 0; j < n; j++) fout << mat[i][j] << " ";
        fout << "\n";
    }
    fout << "\nAll " << n << "! = " << fact << " row-permutations:\n\n";

    vector<int> perm(n), dir(n, -1);
    for (int i = 0; i < n; i++) perm[i] = i;

    int permNum = 0;
    int localSafe = 0;

    while (true) {
        Matrix temp(n);
        for (int i = 0; i < n; i++) temp[i] = mat[perm[i]];

        Matrix Ap = power(temp, (long long)(p - 1), n);
        bool safe = isIdentity(Ap, n);

        permNum++;
        total++;
        if (safe) { safeCount++;  localSafe++; }
        else        unsafeCount++;

        writePermutation(fout, temp, n, p, permNum, safe);

        // Johnson-Trotter: find largest mobile element
        int mobile = -1, mobileIdx = -1;
        for (int i = 0; i < n; i++) {
            int next = i + dir[i];
            if (next >= 0 && next < n && perm[i] > perm[next]) {
                if (perm[i] > mobile) {
                    mobile = perm[i];
                    mobileIdx = i;
                }
            }
        }

        if (mobileIdx == -1) break;

        int swapIdx = mobileIdx + dir[mobileIdx];
        swap(perm[mobileIdx], perm[swapIdx]);
        swap(dir[mobileIdx], dir[swapIdx]);

        for (int i = 0; i < n; i++)
            if (perm[i] > mobile)
                dir[i] *= -1;
    }

    fout << "  >> Source Matrix #" << matrixNum
         << " : SAFE = " << localSafe
         << " | UNSAFE = " << (permNum - localSafe)
         << " (out of " << permNum << " permutations)\n";
    fout << "\n\n";
}

// ─── Main ─────────────────────────────────────────────────────────────────────

int main() {
    ifstream fin("zero_minor_matrices.txt");
    if (!fin) { cerr << "Cannot open zero_minor_matrices.txt\n"; return 1; }

    ofstream fout("combination.txt");

    int safeCount = 0, unsafeCount = 0, total = 0, matrixCount = 0;

    Matrix mat;
    int n, p;
    while (parseNextMatrix(fin, mat, n, p)) {
        matrixCount++;
        MOD = p;
        processPermutations(mat, n, p, fout, safeCount, unsafeCount, total, matrixCount);
    }

    fout << "=============================================================\n";
    fout << "FINAL SUMMARY\n";
    fout << "=============================================================\n";
    fout << "Source matrices read     : " << matrixCount  << "\n";
    fout << "Total permutations tested: " << total        << "\n";
    fout << "SAFE   (A^(p-1) = I)    : " << safeCount    << "\n";
    fout << "UNSAFE                   : " << unsafeCount  << "\n";
    fout << "=============================================================\n";

    fin.close();
    fout.close();

    cout << "Done! Results written to combination.txt\n";
    cout << "Source matrices : " << matrixCount << "\n";
    cout << "Total tested    : " << total       << "\n";
    cout << "Safe            : " << safeCount   << "\n";
    cout << "Unsafe          : " << unsafeCount << "\n";
    return 0;
}