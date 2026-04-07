/*
 * Pattern Matrix Generator using NTL (Number Theory Library)
 *
 * A "Pattern Matrix" of dimension n built on prime p satisfies:
 *   1. All elements are unique within the matrix
 *   2. All elements are in range [1, p-1]  (greater than 0, less than p)
 *   3. Each row sum ≡ (p-1) mod p  (i.e., row sum mod p == p-1)
 *
 * Output is written to BOTH the terminal AND a .txt file automatically.
 * File is named:  pattern_p<P>_n<N>.txt
 */

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vec_ZZ_p.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <set>
#include <climits>

using namespace NTL;
using namespace std;

// ---- dual output: terminal + file -----------------------------------------

ofstream gFile;

void out(const string& s) {
    cout << s;
    if (gFile.is_open()) gFile << s;
}

void outln(const string& s = "") { out(s + "\n"); }

// ---- helpers ---------------------------------------------------------------

bool isPrime(long n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (long i = 3; i * i <= n; i += 2)
        if (n % i == 0) return false;
    return true;
}

// ---- globals ---------------------------------------------------------------

long P_global, N_global, TARGET, TOTAL_CELLS;
vector<vector<long>> current;
vector<long>         rowSum;
std::set<long>       used;
long                 matrixCount;
long                 MAX_PRINT;
bool                 VERBOSE;

// ---- print one matrix ------------------------------------------------------

void printMatrix(const vector<vector<long>>& M, long p) {
    int n = M.size();
    int w = (int)to_string(p - 1).size() + 1;
    string sep = "  +" + string(n * w, '-') + "+";
    outln(sep);
    for (int i = 0; i < n; i++) {
        string row = "  |";
        long s = 0;
        for (int j = 0; j < n; j++) {
            string num = to_string(M[i][j]);
            row += string(w - (int)num.size(), ' ') + num;
            s += M[i][j];
        }
        row += " |   row_sum mod p = " + to_string(s % p);
        outln(row);
    }
    outln(sep);
}

// ---- backtracking ----------------------------------------------------------

void solve(int pos) {
    if (pos == TOTAL_CELLS) {
        matrixCount++;
        if (VERBOSE || matrixCount <= MAX_PRINT) {
            outln();
            outln("  Matrix #" + to_string(matrixCount) + ":");
            printMatrix(current, P_global);
        }
        if (!VERBOSE && matrixCount == MAX_PRINT + 1)
            outln("\n  [Remaining matrices suppressed in terminal — all saved to file.]");
        return;
    }

    int row = pos / N_global;
    int col = pos % N_global;

    for (long v = 1; v < P_global; v++) {
        if (used.count(v)) continue;
        long newSum = rowSum[row] + v;
        // Last cell in row: enforce row-sum constraint exactly
        if (col == N_global - 1 && newSum % P_global != TARGET) continue;

        used.insert(v);
        current[row][col] = v;
        rowSum[row] += v;
        solve(pos + 1);
        used.erase(v);
        current[row][col] = 0;
        rowSum[row] -= v;
    }
}

// ---- main ------------------------------------------------------------------

int main() {
    cout << "=======================================================\n";
    cout << "       Pattern Matrix Generator  (NTL-powered)\n";
    cout << "=======================================================\n\n";

    long p, n;
    cout << "Enter prime p (e.g. 5, 7, 11, 13): ";
    cin >> p;
    if (!isPrime(p)) { cerr << "Error: " << p << " is not prime.\n"; return 1; }

    cout << "Enter matrix dimension n  [n*n must be <= p-1 = " << p-1 << "]: ";
    cin >> n;
    if (n <= 0) { cerr << "Error: n must be positive.\n"; return 1; }
    if (n * n > p - 1) {
        cerr << "Error: n*n=" << n*n << " > p-1=" << p-1
             << ". Not enough distinct values in {1..p-1}.\n";
        return 1;
    }

    cout << "How many matrices to show on screen? (0=count only, -1=all): ";
    cin >> MAX_PRINT;
    VERBOSE = (MAX_PRINT == -1);
    if (VERBOSE) MAX_PRINT = LLONG_MAX;

    // Open output file
    string filename = "output_pattern_17.txt";
    gFile.open(filename);
    if (!gFile.is_open()) {
        cerr << "Warning: Cannot open '" << filename << "' for writing.\n";
    } else {
        cout << "\n  >> Output will ALSO be saved to: " << filename << "\n";
        gFile << "Pattern Matrix Generator — NTL\n";
        gFile << "p = " << p << ", n = " << n << "\n";
        gFile << string(55, '=') << "\n";
    }

    // NTL modular setup
    ZZ_p::init(ZZ(p));

    // Init globals
    P_global    = p;
    N_global    = n;
    TARGET      = p - 1;
    TOTAL_CELLS = n * n;
    matrixCount = 0;
    current.assign(n, vector<long>(n, 0));
    rowSum.assign(n, 0);

    outln();
    outln(string(55, '-'));
    outln("  p           = " + to_string(p));
    outln("  dimension   = " + to_string(n) + " x " + to_string(n));
    outln("  elements    in {1, 2, ..., " + to_string(p-1) + "}");
    outln("  uniqueness  : all " + to_string(n*n) + " entries are distinct");
    outln("  row sum     = " + to_string(p-1) + "  (mod " + to_string(p) + ")");
    outln(string(55, '-'));
    outln();
    outln("Searching for all pattern matrices...");

    solve(0);

    outln();
    outln(string(55, '='));
    outln("  Total pattern matrices found: " + to_string(matrixCount));
    outln(string(55, '='));

    if (gFile.is_open()) {
        gFile.close();
        cout << "\n  >> All " << matrixCount << " matrices saved to: " << filename << "\n\n";
    }

    return 0;
}