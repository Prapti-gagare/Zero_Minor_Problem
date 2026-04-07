#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

// Safe modulo function (handles negative values, long long input)
long long mod(long long x, long long p) {
    long long r = x % p;
    if (r < 0) r += p;
    return r;
}

// Parse a matrix row from the format "  |  1   2  13 |"
// Returns true if exactly 3 values were parsed
bool parseMatrixRow(const string& line, int row[3]) {
    size_t start = line.find('|');
    size_t end   = line.rfind('|');
    if (start == string::npos || end == string::npos || start == end)
        return false;

    string inner = line.substr(start + 1, end - start - 1);
    istringstream ss(inner);
    int val;
    int count = 0;
    while (ss >> val && count < 3)
        row[count++] = val;

    return count == 3;
}

// Check ALL 2x2 minors under mod p
bool hasZeroMinor(int m[3][3], int p) {
    for (int r1 = 0; r1 < 3; r1++) {
        for (int r2 = r1 + 1; r2 < 3; r2++) {
            for (int c1 = 0; c1 < 3; c1++) {
                for (int c2 = c1 + 1; c2 < 3; c2++) {
                    long long a = m[r1][c1], b = m[r1][c2];
                    long long c = m[r2][c1], d = m[r2][c2];
                    if (mod(a * d - b * c, p) == 0)
                        return true;
                }
            }
        }
    }
    return false;
}

// Print matrix
void printMatrix(int m[3][3], ostream& out) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            out << m[i][j] << " ";
        out << "\n";
    }
}

int main() {

    int p;
    cout << "Enter modulus p: ";
    cin >> p;

    ifstream fin("nonsingular.txt");
    ofstream fout("zero_minor_matrices.txt");

    if (!fin) {
        cout << "File not found!\n";
        return 1;
    }

    string line;
    int mat[3][3];
    int rowsParsed = 0;
    int total = 0, count = 0;

    while (getline(fin, line)) {

        // Only process data rows: contain '|' but are NOT border lines like "+-----------+"
        if (line.find('|') == string::npos) continue;
        if (line.find("+-") != string::npos) continue;

        int row[3];
        if (!parseMatrixRow(line, row)) continue;

        // Store this data row
        for (int j = 0; j < 3; j++)
            mat[rowsParsed][j] = row[j];
        rowsParsed++;

        // Once we have all 3 rows, process the matrix immediately
        if (rowsParsed == 3) {
            total++;
            if (hasZeroMinor(mat, p)) {
                count++;
                cout << "Matrix with zero minor (mod " << p << "):\n";
                printMatrix(mat, cout);
                cout << "-------------------\n";

                fout << "Matrix (mod " << p << "):\n";
                printMatrix(mat, fout);
                fout << "-------------------\n";
            }
            rowsParsed = 0;
        }
    }

    cout << "\nTotal matrices processed = " << total << "\n";
    cout << "Matrices with zero minor (mod " << p << ") = " << count << "\n";

    return 0;
}