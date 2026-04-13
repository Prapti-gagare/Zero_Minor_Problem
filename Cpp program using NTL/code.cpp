/*
 * pipeline.cpp
 *
 * Unified orchestrator for the Pattern Matrix pipeline.
 *
 * Pipeline:
 *   For each prime p in primes.txt:
 *     Step 1: check2by2  → output_pattern_<p>.txt  (split at 500 MB)
 *     Step 2: seperate   → nonsingular_<p>.txt / singular_<p>.txt  (split at 500 MB)
 *     Step 3: zero       → zero_minor_<p>.txt
 *     Step 4: dist       → distinct_<p>.txt
 *     Step 5: seperate   → nonsingular_dist_<p>.txt / singular_dist_<p>.txt
 *     Step 6: zero       → zero_minor_dist_<p>.txt
 *   Final:   final_summary.txt  (one row per prime)
 *
 * Build:
 *   g++ -O2 -std=c++17 pipeline.cpp -lNTL -lgmp -o pipeline
 *
 * Notes:
 *   - NTL is only linked for Step 1 (the pattern-matrix generator).
 *   - All other steps are self-contained and embedded below.
 *   - Max file size before splitting: 500 MB.
 */

// ─────────────────────────────────────────────────────────────────────────────
//  Standard headers
// ─────────────────────────────────────────────────────────────────────────────

#include <NTL/ZZ_p.h>          // only used in step 1

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
//  Constants
// ─────────────────────────────────────────────────────────────────────────────

static const long long MAX_FILE_BYTES = 500LL * 1024 * 1024;   // 500 MB

// ─────────────────────────────────────────────────────────────────────────────
//  Shared utility: split-aware output writer
// ─────────────────────────────────────────────────────────────────────────────

struct SplitWriter {
    string base;          // e.g. "output_pattern_17"
    string ext  = ".txt"; // always .txt
    int    part = 0;      // 0 → no suffix yet
    long long bytes = 0;
    ofstream f;

    // Open the first file
    SplitWriter(const string& baseName) : base(baseName) { openNext(); }

    void openNext() {
        if (f.is_open()) f.close();
        string name = part == 0 ? base + ext
                                : base + "_" + to_string(part) + ext;
        f.open(name);
        if (!f.is_open()) throw runtime_error("Cannot open " + name);
        bytes = 0;
        ++part;
    }

    // Return the current filename (for later reading)
    string currentName() const {
        int p = part - 1;
        return p <= 1 ? base + ext
                      : base + "_" + to_string(p) + ext;
    }

    // All filenames written so far
    vector<string> allNames;

    void write(const string& s) {
        if (bytes + (long long)s.size() > MAX_FILE_BYTES) {
            allNames.push_back(currentName());
            openNext();
        }
        f << s;
        bytes += (long long)s.size();
    }

    void writeln(const string& s = "") { write(s + "\n"); }

    void close() {
        if (f.is_open()) {
            allNames.push_back(currentName());
            f.close();
        }
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Utility helpers
// ─────────────────────────────────────────────────────────────────────────────

bool isPrime(long n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (long i = 3; i * i <= n; i += 2)
        if (n % i == 0) return false;
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
//  STEP 1 — Pattern matrix generator (check2by2 logic)
// ─────────────────────────────────────────────────────────────────────────────

namespace Step1 {

long P, N, TARGET, TOTAL_CELLS;
vector<vector<long>> current;
vector<long>         rowSum;
set<long>            used;
long                 matrixCount;
SplitWriter*         gWriter;

void printMatrix(const vector<vector<long>>& M, long p) {
    int n = M.size();
    int w = (int)to_string(p - 1).size() + 1;
    string sep = "  +" + string(n * w, '-') + "+";
    gWriter->writeln(sep);
    for (int i = 0; i < n; i++) {
        string row = "  |";
        long s = 0;
        for (int j = 0; j < n; j++) {
            string num = to_string(M[i][j]);
            row += string(w - (int)num.size(), ' ') + num;
            s += M[i][j];
        }
        row += " |   row_sum mod p = " + to_string(s % p);
        gWriter->writeln(row);
    }
    gWriter->writeln(sep);
}

void solve(int pos) {
    if (pos == TOTAL_CELLS) {
        matrixCount++;
        gWriter->writeln();
        gWriter->writeln("  Matrix #" + to_string(matrixCount) + ":");
        printMatrix(current, P);
        return;
    }

    int row = pos / N;
    int col = pos % N;

    for (long v = 1; v < P; v++) {
        if (used.count(v)) continue;
        long newSum = rowSum[row] + v;
        if (col == N - 1 && newSum % P != TARGET) continue;
        // Pruning: remaining cells can't possibly reach target
        int remaining = N - col - 1;
        if (col < N - 1) {
            long minAdd = 0, maxAdd = 0;
            int cnt = 0;
            for (long c = 1; c < P && cnt < remaining; c++)
                if (!used.count(c) && c != v) { minAdd += c; cnt++; }
            cnt = 0;
            for (long c = P - 1; c >= 1 && cnt < remaining; c--)
                if (!used.count(c) && c != v) { maxAdd += c; cnt++; }
            long lo = (newSum + minAdd) % P;
            long hi = (newSum + maxAdd) % P;
            (void)lo; (void)hi;  // basic pruning only on last cell
        }

        used.insert(v);
        current[row][col] = v;
        rowSum[row] += v;
        solve(pos + 1);
        used.erase(v);
        current[row][col] = 0;
        rowSum[row] -= v;
    }
}

// Returns total matrices found, fills writerFiles with all output filenames
long run(long p, long n, const string& baseOut, vector<string>& writerFiles) {
    NTL::ZZ_p::init(NTL::ZZ(p));

    P           = p;
    N           = n;
    TARGET      = p - 1;
    TOTAL_CELLS = n * n;
    matrixCount = 0;
    current.assign(n, vector<long>(n, 0));
    rowSum.assign(n, 0);
    used.clear();

    SplitWriter writer(baseOut);
    gWriter = &writer;

    writer.writeln("Pattern Matrix Generator");
    writer.writeln(string(55, '='));
    writer.writeln("  p           = " + to_string(p));
    writer.writeln("  dimension   = " + to_string(n) + " x " + to_string(n));
    writer.writeln("  elements    in {1, 2, ..., " + to_string(p - 1) + "}");
    writer.writeln("  uniqueness  : all " + to_string(n * n) + " entries are distinct");
    writer.writeln("  row sum     = " + to_string(p - 1) + "  (mod " + to_string(p) + ")");
    writer.writeln(string(55, '-'));
    writer.writeln();

    solve(0);

    writer.writeln();
    writer.writeln(string(55, '='));
    writer.writeln("  Total pattern matrices found: " + to_string(matrixCount));
    writer.writeln(string(55, '='));
    writer.close();

    writerFiles = writer.allNames;
    return matrixCount;
}

} // namespace Step1

// ─────────────────────────────────────────────────────────────────────────────
//  STEP 2 — Singular / non-singular separator (seperate.cpp logic)
// ─────────────────────────────────────────────────────────────────────────────

namespace Step2 {

struct Matrix {
    int id;
    vector<vector<double>> data;
};

long long modpow(long long base, long long exp, long long mod) {
    long long result = 1; base %= mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

long long modinv(long long a, long long p) {
    return modpow(((a % p) + p) % p, p - 2, p);
}

long long determinantModP(vector<vector<double>> matD, long long p) {
    int n = matD.size();
    vector<vector<long long>> mat(n, vector<long long>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mat[i][j] = (((long long)matD[i][j] % p) + p) % p;

    long long det = 1;
    for (int col = 0; col < n; col++) {
        int pivot = -1;
        for (int row = col; row < n; row++)
            if (mat[row][col] != 0) { pivot = row; break; }
        if (pivot == -1) return 0;
        if (pivot != col) { swap(mat[pivot], mat[col]); det = (p - det) % p; }
        det = det * mat[col][col] % p;
        long long inv = modinv(mat[col][col], p);
        for (int row = col + 1; row < n; row++) {
            long long factor = mat[row][col] * inv % p;
            for (int k = col; k < n; k++)
                mat[row][k] = (mat[row][k] - factor * mat[col][k] % p + p) % p;
        }
    }
    return det;
}

vector<Matrix> parseMatrices(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) return {};

    vector<Matrix> matrices;
    string line;
    regex matrixHeader(R"(Matrix\s+#(\d+):)");
    regex matrixRow(R"(\|\s*((?:\d+\s*)+)\|)");
    Matrix current;
    bool inMatrix = false;

    while (getline(file, line)) {
        smatch m;
        if (regex_search(line, m, matrixHeader)) {
            if (inMatrix && !current.data.empty()) matrices.push_back(current);
            current = Matrix(); current.id = stoi(m[1]); current.data.clear();
            inMatrix = true; continue;
        }
        if (inMatrix && regex_search(line, m, matrixRow)) {
            istringstream iss(m[1].str());
            vector<double> row; double val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) current.data.push_back(row);
        }
    }
    if (inMatrix && !current.data.empty()) matrices.push_back(current);
    return matrices;
}

void writeMatrixToWriter(SplitWriter& sw, const Matrix& mat, long long det, long long p) {
    int maxVal = 0;
    for (const auto& row : mat.data) for (double v : row) if ((int)v > maxVal) maxVal = (int)v;
    int width = to_string(maxVal).size() + 1;
    int cols = mat.data.empty() ? 0 : mat.data[0].size();
    string border = "  +-" + string(cols * width, '-') + "-+";

    sw.writeln("  Matrix #" + to_string(mat.id) + ":");
    sw.writeln(border);
    for (const auto& row : mat.data) {
        string line = "  |";
        for (double v : row) {
            string s = to_string((int)v);
            line += string(width - s.size(), ' ') + s + " ";
        }
        line += "|";
        sw.writeln(line);
    }
    sw.writeln(border);
    sw.writeln("  Determinant (mod " + to_string(p) + ") = " + to_string(det));
    sw.writeln();
}

// Returns {singCount, nonSingCount}
// Fills singFiles and nonSingFiles with actual filenames written
pair<long,long> run(const vector<string>& inputFiles,
                    const string& singBase,
                    const string& nonSingBase,
                    long long p,
                    vector<string>& singFiles,
                    vector<string>& nonSingFiles) {

    SplitWriter singWriter(singBase);
    SplitWriter nonSingWriter(nonSingBase);

    singWriter.writeln("===========================================");
    singWriter.writeln("  SINGULAR MATRICES (det = 0)");
    singWriter.writeln("===========================================");
    singWriter.writeln();

    nonSingWriter.writeln("===========================================");
    nonSingWriter.writeln("  NON-SINGULAR MATRICES (det != 0)");
    nonSingWriter.writeln("===========================================");
    nonSingWriter.writeln();

    long singCount = 0, nonSingCount = 0;

    for (const auto& fname : inputFiles) {
        auto matrices = parseMatrices(fname);
        for (const auto& mat : matrices) {
            long long det = determinantModP(mat.data, p);
            if (det == 0) {
                writeMatrixToWriter(singWriter, mat, det, p);
                singCount++;
            } else {
                writeMatrixToWriter(nonSingWriter, mat, det, p);
                nonSingCount++;
            }
        }
    }

    singWriter.writeln("-------------------------------------------");
    singWriter.writeln("Total singular matrices     : " + to_string(singCount));
    singWriter.close();

    nonSingWriter.writeln("-------------------------------------------");
    nonSingWriter.writeln("Total non-singular matrices : " + to_string(nonSingCount));
    nonSingWriter.close();

    singFiles    = singWriter.allNames;
    nonSingFiles = nonSingWriter.allNames;
    return {singCount, nonSingCount};
}

} // namespace Step2

// ─────────────────────────────────────────────────────────────────────────────
//  STEP 3 — Zero minor checker (zero.cpp logic)
// ─────────────────────────────────────────────────────────────────────────────

namespace Step3 {

long long mod(long long x, long long p) {
    long long r = x % p; if (r < 0) r += p; return r;
}

bool parseMatrixRow(const string& line, int row[3]) {
    size_t start = line.find('|');
    size_t end   = line.rfind('|');
    if (start == string::npos || end == string::npos || start == end) return false;
    string inner = line.substr(start + 1, end - start - 1);
    istringstream ss(inner);
    int val, count = 0;
    while (ss >> val && count < 3) row[count++] = val;
    return count == 3;
}

bool hasZeroMinor(int m[3][3], long long p) {
    for (int r1 = 0; r1 < 3; r1++)
        for (int r2 = r1 + 1; r2 < 3; r2++)
            for (int c1 = 0; c1 < 3; c1++)
                for (int c2 = c1 + 1; c2 < 3; c2++) {
                    long long a = m[r1][c1], b = m[r1][c2];
                    long long c = m[r2][c1], d = m[r2][c2];
                    if (mod(a * d - b * c, p) == 0) return true;
                }
    return false;
}

// Returns count of zero-minor matrices
long run(const vector<string>& nonSingFiles, const string& outBase, long long p) {
    SplitWriter writer(outBase);
    writer.writeln("Zero Minor Matrices (mod " + to_string(p) + ")");
    writer.writeln(string(45, '='));
    writer.writeln();

    long total = 0, count = 0;
    for (const auto& fname : nonSingFiles) {
        ifstream fin(fname);
        if (!fin.is_open()) continue;
        string line;
        int mat[3][3];
        int rowsParsed = 0;

        while (getline(fin, line)) {
            if (line.find('|') == string::npos) continue;
            if (line.find("+-") != string::npos) continue;
            int row[3];
            if (!parseMatrixRow(line, row)) continue;
            for (int j = 0; j < 3; j++) mat[rowsParsed][j] = row[j];
            rowsParsed++;
            if (rowsParsed == 3) {
                total++;
                if (hasZeroMinor(mat, p)) {
                    count++;
                    writer.writeln("Matrix with zero minor (mod " + to_string(p) + "):");
                    for (int i = 0; i < 3; i++) {
                        writer.writeln(to_string(mat[i][0]) + " " +
                                       to_string(mat[i][1]) + " " +
                                       to_string(mat[i][2]));
                    }
                    writer.writeln("-------------------");
                }
                rowsParsed = 0;
            }
        }
    }

    writer.writeln();
    writer.writeln("Total matrices processed = " + to_string(total));
    writer.writeln("Matrices with zero minor = " + to_string(count));
    writer.close();
    return count;
}

} // namespace Step3

// ─────────────────────────────────────────────────────────────────────────────
//  STEP 4 — Distinct matrix filter (dist.cpp logic)
// ─────────────────────────────────────────────────────────────────────────────

namespace Step4 {

using Row    = vector<int>;
using Matrix = vector<Row>;

Matrix canonical(Matrix m) {
    for (auto& row : m) sort(row.begin(), row.end());
    sort(m.begin(), m.end());
    return m;
}

vector<int> extractInts(const string& s) {
    vector<int> nums;
    istringstream iss(s);
    string token;
    while (iss >> token) {
        size_t start = token.find_first_of("0123456789");
        if (start == string::npos) continue;
        try { nums.push_back(stoi(token.substr(start))); } catch (...) {}
    }
    return nums;
}

vector<Matrix> parseMatrices(const vector<string>& files, long p, long n) {
    vector<Matrix> result;
    regex matrixHeader(R"(Matrix\s+#(\d+):)");

    for (const auto& filename : files) {
        ifstream in(filename);
        if (!in.is_open()) continue;
        string line;
        Matrix current;
        bool inMatrix = false;

        while (getline(in, line)) {
            if (line.find('|') != string::npos) {
                size_t first = line.find('|');
                size_t last  = line.rfind('|');
                if (first != last && last > first + 1) {
                    string inner = line.substr(first + 1, last - first - 1);
                    auto nums = extractInts(inner);
                    if (!nums.empty()) {
                        if (!inMatrix) inMatrix = true;
                        current.push_back(nums);
                        if (n > 0 && (int)current.size() == (int)n) {
                            result.push_back(current);
                            current.clear();
                            inMatrix = false;
                        }
                    }
                }
            } else if (inMatrix && !current.empty()) {
                smatch m;
                if (regex_search(line, m, matrixHeader) ||
                    line.find("+-") != string::npos) {
                    current.clear(); inMatrix = false;
                }
            }
        }
        if (!current.empty()) result.push_back(current);
    }
    return result;
}

// Returns count of distinct matrices; writes them to outBase
long run(const vector<string>& inputFiles, const string& outBase, long p, long n) {
    auto matrices = parseMatrices(inputFiles, p, n);

    set<Matrix> seen;
    vector<Matrix> distinct;
    for (const auto& m : matrices) {
        Matrix key = canonical(m);
        if (seen.insert(key).second) distinct.push_back(m);
    }

    SplitWriter writer(outBase);
    writer.writeln("Distinct Matrix Filter");
    writer.writeln(string(55, '='));
    writer.writeln("  p           = " + to_string(p));
    writer.writeln("  dimension   = " + to_string(n) + " x " + to_string(n));
    writer.writeln("  total input matrices  = " + to_string(matrices.size()));
    writer.writeln("  distinct matrices     = " + to_string(distinct.size()));
    writer.writeln(string(55, '-'));
    writer.writeln();

    int idx = 1;
    for (const auto& m : distinct) {
        int maxVal = 0;
        for (const auto& row : m) for (int v : row) if (v > maxVal) maxVal = v;
        int width = to_string(maxVal).length() + 1;
        int cols  = m.empty() ? 0 : (int)m[0].size();
        int borderW = cols * (width + 1) + 1;
        string border(borderW, '-');
        if (!border.empty()) { border[0] = '+'; border.back() = '+'; }

        writer.writeln("  Matrix #" + to_string(idx++) + ":");
        writer.writeln("  " + border);
        for (const auto& row : m) {
            string line = "  |";
            int s = 0;
            for (int v : row) {
                string sv = to_string(v);
                line += string(width - sv.size(), ' ') + sv + " ";
                s += v;
            }
            line += "|   row_sum mod " + to_string(p) + " = " + to_string(s % p);
            writer.writeln(line);
        }
        writer.writeln("  " + border);
        writer.writeln();
    }

    writer.writeln(string(55, '='));
    writer.writeln("End of distinct matrices.");
    writer.close();

    return (long)distinct.size();
}

} // namespace Step4

// ─────────────────────────────────────────────────────────────────────────────
//  RAM usage helper (Linux /proc/self/status)
// ─────────────────────────────────────────────────────────────────────────────

long long getRamUsageKB() {
    ifstream status("/proc/self/status");
    if (!status.is_open()) return 0;
    string line;
    while (getline(status, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            istringstream iss(line);
            string label; long long kb;
            iss >> label >> kb;
            return kb;
        }
    }
    return 0;
}

string formatRAM(long long kb) {
    if (kb <= 0) return "N/A";
    if (kb < 1024)
        return to_string(kb) + " KB";
    if (kb < 1024 * 1024)
        return to_string(kb / 1024) + " MB (" + to_string(kb) + " KB)";
    return to_string(kb / (1024 * 1024)) + " GB (" + to_string(kb / 1024) + " MB)";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Summary row
// ─────────────────────────────────────────────────────────────────────────────

struct SummaryRow {
    long prime;
    long totalMatrices;
    long nonSingular;
    long zeroMinor;
    long distinctTotal;
    long distinctNonSingular;
    long distinctZeroMinor;
    long long ramPeakKB = 0;   // RAM (RSS) at end of this prime's processing
};

// ─────────────────────────────────────────────────────────────────────────────
//  Main orchestrator
// ─────────────────────────────────────────────────────────────────────────────

int main() {
    // ── Read config ──────────────────────────────────────────────────────────
    cout << "========================================\n";
    cout << "   Pattern Matrix Pipeline Orchestrator\n";
    cout << "========================================\n\n";

    string primesFile;
    cout << "Enter primes input file (e.g. primes.txt): ";
    cin >> primesFile;

    long n;
    cout << "Enter matrix dimension n (for check2by2, n x n): ";
    cin >> n;

    // ── Load primes ───────────────────────────────────────────────────────────
    ifstream pfile(primesFile);
    if (!pfile.is_open()) { cerr << "Cannot open " << primesFile << "\n"; return 1; }

    vector<long> primes;
    long v;
    while (pfile >> v) {
        if (!isPrime(v)) { cerr << "Warning: " << v << " is not prime, skipping.\n"; continue; }
        if (n * n > v - 1) {
            cerr << "Warning: for p=" << v << ", n*n=" << n*n
                 << " > p-1=" << v-1 << ", skipping.\n";
            continue;
        }
        primes.push_back(v);
    }
    pfile.close();

    if (primes.empty()) { cerr << "No valid primes found.\n"; return 1; }

    cout << "\nFound " << primes.size() << " valid prime(s): ";
    for (long pr : primes) cout << pr << " ";
    cout << "\n\n";

    // ── Process each prime ────────────────────────────────────────────────────
    vector<SummaryRow> summary;

    for (long p : primes) {
        cout << "────────────────────────────────────────────\n";
        cout << "  Processing prime p = " << p << "\n";
        cout << "────────────────────────────────────────────\n";

        SummaryRow row{};
        row.prime = p;

        // Step 1: generate pattern matrices
        cout << "  [Step 1] Generating pattern matrices...\n";
        string s1base = "output_pattern_" + to_string(p);
        vector<string> s1files;
        row.totalMatrices = Step1::run(p, n, s1base, s1files);
        cout << "           → " << row.totalMatrices << " matrices in "
             << s1files.size() << " file(s)\n";

        // Step 2: separate singular / non-singular (from all matrices)
        cout << "  [Step 2] Classifying singular/non-singular (all matrices)...\n";
        string s2singBase    = "singular_" + to_string(p);
        string s2nonSingBase = "nonsingular_" + to_string(p);
        vector<string> singFiles, nonSingFiles;
        pair<long,long> s2result = Step2::run(
            s1files, s2singBase, s2nonSingBase, (long long)p, singFiles, nonSingFiles);
        long singCnt    = s2result.first;
        long nonSingCnt = s2result.second;
        row.nonSingular = nonSingCnt;
        cout << "           → singular=" << singCnt
             << "  non-singular=" << nonSingCnt << "\n";

        // Step 3: zero minors among non-singular
        cout << "  [Step 3] Finding zero minors (all non-singular)...\n";
        string s3base = "zero_minor_" + to_string(p);
        row.zeroMinor = Step3::run(nonSingFiles, s3base, (long long)p);
        cout << "           → " << row.zeroMinor << " zero-minor matrices\n";

        // Step 4: distinct matrices from all pattern matrices
        cout << "  [Step 4] Finding distinct matrices...\n";
        string s4base = "distinct_" + to_string(p);
        row.distinctTotal = Step4::run(s1files, s4base, p, n);
        cout << "           → " << row.distinctTotal << " distinct matrices\n";

        // Step 5: separate singular / non-singular among distinct matrices
        cout << "  [Step 5] Classifying distinct matrices...\n";
        string s5singBase    = "singular_dist_" + to_string(p);
        string s5nonSingBase = "nonsingular_dist_" + to_string(p);
        vector<string> distSingFiles, distNonSingFiles;
        // Re-parse distinct output for step 5
        vector<string> distFiles = { s4base + ".txt" };
        // Add split parts if they exist
        for (int part = 2; ; part++) {
            string f = s4base + "_" + to_string(part) + ".txt";
            if (!ifstream(f).good()) break;
            distFiles.push_back(f);
        }
        pair<long,long> s5result = Step2::run(
            distFiles, s5singBase, s5nonSingBase, (long long)p,
            distSingFiles, distNonSingFiles);
        long dSingCnt    = s5result.first;
        long dNonSingCnt = s5result.second;
        row.distinctNonSingular = dNonSingCnt;
        cout << "           → distinct singular=" << dSingCnt
             << "  distinct non-singular=" << dNonSingCnt << "\n";

        // Step 6: zero minors among distinct non-singular
        cout << "  [Step 6] Zero minors among distinct non-singular...\n";
        string s6base = "zero_minor_dist_" + to_string(p);
        row.distinctZeroMinor = Step3::run(distNonSingFiles, s6base, (long long)p);
        cout << "           → " << row.distinctZeroMinor << " zero-minor matrices\n";

        row.ramPeakKB = getRamUsageKB();
        cout << "  [RAM]  " << formatRAM(row.ramPeakKB) << " used after processing p=" << p << "\n\n";

        summary.push_back(row);

        // ── Write per-prime summary file ─────────────────────────────────────
        {
            string perPrimeFile = "final_summary_" + to_string(p) + ".txt";
            cout << "  Writing " << perPrimeFile << "...\n\n";
            ofstream out(perPrimeFile);
            if (!out.is_open()) {
                cerr << "Cannot open " << perPrimeFile << "\n";
            } else {
                out << "Pattern Matrix Pipeline — Per-Prime Summary\n";
                out << string(60, '=') << "\n";
                out << "  Prime p = " << p << "\n";
                out << string(60, '-') << "\n\n";

                out << "  [Matrix Counts]\n";
                out << "    Total pattern matrices generated : " << row.totalMatrices << "\n";
                out << "    Non-singular matrices            : " << row.nonSingular << "\n";
                out << "    Singular matrices                : " << (row.totalMatrices - row.nonSingular) << "\n";
                out << "    Zero-minor (non-singular)        : " << row.zeroMinor << "\n\n";

                out << "  [Distinct Matrix Counts]\n";
                out << "    Distinct matrices total          : " << row.distinctTotal << "\n";
                out << "    Distinct non-singular            : " << row.distinctNonSingular << "\n";
                out << "    Distinct singular                : " << (row.distinctTotal - row.distinctNonSingular) << "\n";
                out << "    Distinct zero-minor              : " << row.distinctZeroMinor << "\n\n";

                out << "  [Resource Usage]\n";
                out << "    RAM used (RSS at completion)     : " << formatRAM(row.ramPeakKB) << "\n\n";

                out << string(60, '=') << "\n";
                out << "\nColumn Descriptions:\n";
                out << "  Total pattern matrices  : All matrices from check2by2 for this prime\n";
                out << "  Non-singular            : Matrices with det != 0 (mod p)\n";
                out << "  Zero-minor              : Non-singular matrices with a zero 2x2 minor (mod p)\n";
                out << "  Distinct                : Unique up to row/column/element permutations\n";
                out << "  RAM (RSS)               : Resident Set Size — physical RAM held by the\n";
                out << "                            process at the moment this prime finished.\n";
                out.close();
            }
        }
    }

    cout << "\n========================================\n";
    cout << "  Pipeline complete.\n";
    cout << "  Per-prime summaries written:\n";
    for (const auto& r : summary)
        cout << "    final_summary_" << r.prime << ".txt\n";
    cout << "========================================\n";

    return 0;
}