#include <bits/stdc++.h>
using namespace std;

// Print matrix
void writeMatrix(ofstream &fout, const vector<vector<int>> &mat, int n, int p) {
    fout << n << " " << p << "\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fout << mat[i][j] << " ";
        }
        fout << "\n";
    }
    fout << "\n";
}

// Generate permutations using adjacent swaps (Johnson–Trotter)
void generatePermutations(vector<vector<int>> mat, int n, int p, ofstream &fout) {
    vector<int> perm(n), dir(n, -1); // -1 = left, 1 = right

    for (int i = 0; i < n; i++) perm[i] = i;

    while (true) {
        // Print current permutation
        vector<vector<int>> temp(n);
        for (int i = 0; i < n; i++) {
            temp[i] = mat[perm[i]];
        }
        writeMatrix(fout, temp, n, p);

        // Find largest mobile element
        int mobile = -1, mobileIndex = -1;
        for (int i = 0; i < n; i++) {
            int next = i + dir[i];
            if (next >= 0 && next < n && perm[i] > perm[next]) {
                if (perm[i] > mobile) {
                    mobile = perm[i];
                    mobileIndex = i;
                }
            }
        }

        if (mobileIndex == -1) break;

        int swapIndex = mobileIndex + dir[mobileIndex];
        swap(perm[mobileIndex], perm[swapIndex]);
        swap(dir[mobileIndex], dir[swapIndex]);

        // Reverse direction of larger elements
        for (int i = 0; i < n; i++) {
            if (perm[i] > mobile) {
                dir[i] *= -1;
            }
        }
    }
}

int main() {
    ifstream fin("shuffle1.txt");
    ofstream fout("output.txt");

    int n, p;

    while (fin >> n >> p) {
        vector<vector<int>> mat(n, vector<int>(n));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                fin >> mat[i][j];

        generatePermutations(mat, n, p, fout);
    }

    fin.close();
    fout.close();

    cout << "Output written to output.txt\n";
    return 0;
}  