/*
 * distinct_matrices.cpp
 *
 * Reads pattern matrices from a file (NTL-style format), then finds all
 * "distinct" matrices — i.e. removes duplicates that are merely row-permutations,
 * column-permutations, or within-row element-permutations of each other,
 * while keeping the per-row sum invariant.
 *
 * Two matrices A and B are considered the SAME (non-distinct) if B can be
 * obtained from A by any combination of:
 *   1. Permuting rows
 *   2. Permuting columns
 *   3. Permuting elements within any individual row
 *
 * Because column permutations and within-row permutations together simply mean
 * "permute the multiset of elements in each row in any way across any row",
 * the canonical representative used here is:
 *   - Sort the elements within each row.
 *   - Sort the rows lexicographically.
 * Two matrices that produce the same canonical form are equivalent.
 *
 * Usage:
 *   ./distinct_matrices <input_file> <output_file>
 *   (defaults: input.txt  output.txt)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <iomanip>
#include <regex>
#include <stdexcept>

// ─── Data types ───────────────────────────────────────────────────────────────

using Row    = std::vector<int>;
using Matrix = std::vector<Row>;

struct ParsedInput {
    int p   = 0;   // prime modulus
    int n   = 0;   // dimension
    int target_sum = 0;   // expected row sum mod p
    std::vector<Matrix> matrices;
};

// ─── Canonical form ───────────────────────────────────────────────────────────

/*
 * Canonical form:
 *   1. Sort each row independently (elements within a row may be freely reordered).
 *   2. Sort all rows lexicographically.
 * This single key captures equivalence under any row permutation AND any
 * independent column permutation (within each row).
 */
Matrix canonical(Matrix m) {
    for (auto& row : m)
        std::sort(row.begin(), row.end());
    std::sort(m.begin(), m.end());
    return m;
}

// ─── Parsing ──────────────────────────────────────────────────────────────────

/*
 * Extracts all integers from a string.
 */
std::vector<int> extractInts(const std::string& s) {
    std::vector<int> nums;
    std::istringstream iss(s);
    std::string token;
    while (iss >> token) {
        try {
            // strip non-numeric prefix/suffix characters
            size_t start = token.find_first_of("0123456789");
            if (start == std::string::npos) continue;
            nums.push_back(std::stoi(token.substr(start)));
        } catch (...) {}
    }
    return nums;
}

ParsedInput parseFile(const std::string& filename) {
    std::ifstream in(filename);
    if (!in.is_open())
        throw std::runtime_error("Cannot open input file: " + filename);

    ParsedInput result;
    std::string line;
    Matrix current;
    bool inMatrix = false;

    while (std::getline(in, line)) {
        // ── header: p = 17
        if (result.p == 0) {
            auto pos = line.find("p =");
            if (pos != std::string::npos) {
                auto nums = extractInts(line.substr(pos + 3));
                if (!nums.empty()) result.p = nums[0];
            }
        }
        // ── header: dimension = 3 x 3  OR  n = 3
        if (result.n == 0) {
            auto pos = line.find("dimension");
            if (pos != std::string::npos) {
                auto nums = extractInts(line);
                if (!nums.empty()) result.n = nums[0];
            } else {
                auto pos2 = line.find("n =");
                if (pos2 != std::string::npos) {
                    auto nums = extractInts(line.substr(pos2 + 3));
                    if (!nums.empty()) result.n = nums[0];
                }
            }
        }
        // ── header: row sum = 16
        if (result.target_sum == 0) {
            auto pos = line.find("row sum");
            if (pos != std::string::npos) {
                auto nums = extractInts(line);
                if (!nums.empty()) result.target_sum = nums[0];
            }
        }

        // ── matrix row: lines containing '|'  (e.g.  |  1  2 13 |)
        if (line.find('|') != std::string::npos) {
            // extract numbers between the pipes
            size_t first = line.find('|');
            size_t last  = line.rfind('|');
            if (first != last && last > first + 1) {
                std::string inner = line.substr(first + 1, last - first - 1);
                auto nums = extractInts(inner);
                if (!nums.empty()) {
                    if (!inMatrix) inMatrix = true;
                    current.push_back(nums);
                    if (result.n > 0 && (int)current.size() == result.n) {
                        result.matrices.push_back(current);
                        current.clear();
                        inMatrix = false;
                    }
                }
            }
        } else if (inMatrix && (line.find("+-") != std::string::npos ||
                                line.find("Matrix") != std::string::npos)) {
            // reset partial accumulation on unexpected boundary
            if (!current.empty()) {
                current.clear();
                inMatrix = false;
            }
        }
    }
    // push any leftover partial matrix
    if (!current.empty())
        result.matrices.push_back(current);

    if (result.p == 0) result.p = 17;   // sane defaults
    if (result.n == 0) result.n = 3;

    return result;
}

// ─── Distinct-filter ──────────────────────────────────────────────────────────

std::vector<Matrix> findDistinct(const std::vector<Matrix>& matrices) {
    std::set<Matrix> seen;
    std::vector<Matrix> distinct;
    for (const auto& m : matrices) {
        Matrix key = canonical(m);
        if (seen.insert(key).second)   // true  ⇒ not seen before
            distinct.push_back(m);
    }
    return distinct;
}

// ─── Validation ───────────────────────────────────────────────────────────────

bool validateMatrix(const Matrix& m, int p, int target_sum) {
    for (const auto& row : m) {
        int s = 0;
        for (int v : row) s += v;
        if ((s % p) != (target_sum % p)) return false;
    }
    return true;
}

// ─── Pretty printing ──────────────────────────────────────────────────────────

void printMatrix(std::ostream& out, const Matrix& m, int idx, int p, int target_sum) {
    // measure column width
    int maxVal = 0;
    for (const auto& row : m)
        for (int v : row)
            if (v > maxVal) maxVal = v;
    int width = std::to_string(maxVal).length() + 1;

    // border width
    int n      = (int)m.size();
    int cols   = n > 0 ? (int)m[0].size() : 0;
    int borderW = cols * (width + 1) + 1;
    std::string border(borderW, '-');
    border[0] = '+'; border.back() = '+';

    out << "  Matrix #" << idx << ":\n";
    out << "  " << border << "\n";
    for (const auto& row : m) {
        out << "  |";
        for (int v : row)
            out << std::setw(width) << v << " ";
        out << "|";
        int s = 0;
        for (int v : row) s += v;
        out << "   row_sum mod " << p << " = " << (s % p) << "\n";
    }
    out << "  " << border << "\n";
}

void writeOutput(const std::string& filename,
                 const ParsedInput& input,
                 const std::vector<Matrix>& allMatrices,
                 const std::vector<Matrix>& distinct) {

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error("Cannot open output file: " + filename);

    out << "Distinct Matrix Filter\n";
    out << "=======================================================\n";
    out << "-------------------------------------------------------\n";
    out << "  p           = " << input.p << "\n";
    out << "  dimension   = " << input.n << " x " << input.n << "\n";
    out << "  total input matrices  = " << allMatrices.size() << "\n";
    out << "  distinct matrices     = " << distinct.size() << "\n";
    out << "  row sum     = " << input.target_sum << "  (mod " << input.p << ")\n";
    out << "-------------------------------------------------------\n\n";

    out << "Distinct matrices (canonical representatives):\n\n";
    for (int i = 0; i < (int)distinct.size(); ++i)
        printMatrix(out, distinct[i], i + 1, input.p, input.target_sum);

    out << "\n=======================================================\n";
    out << "End of distinct matrices.\n";

    std::cout << "Written " << distinct.size()
              << " distinct matrices to \"" << filename << "\".\n";
}

// ─── Main ─────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    std::string inFile  = (argc > 1) ? argv[1] : "output_pattern_1.txt";
    std::string outFile = (argc > 2) ? argv[2] : "distinct_19.txt";

    try {
        // 1. Parse
        std::cout << "Reading matrices from \"" << inFile << "\"...\n";
        ParsedInput input = parseFile(inFile);
        std::cout << "  p = " << input.p
                  << ", n = " << input.n
                  << ", target row-sum = " << input.target_sum << "\n";
        std::cout << "  Total matrices parsed: " << input.matrices.size() << "\n";

        // 2. Validate (warn on bad matrices but don't drop them)
        int bad = 0;
        for (const auto& m : input.matrices)
            if (!validateMatrix(m, input.p, input.target_sum)) ++bad;
        if (bad > 0)
            std::cerr << "  Warning: " << bad
                      << " matrix(ces) fail the row-sum check.\n";

        // 3. Find distinct
        std::vector<Matrix> distinct = findDistinct(input.matrices);
        std::cout << "  Distinct matrices: " << distinct.size() << "\n";

        // 4. Print to stdout
        std::cout << "\n--- Distinct Matrices ---\n\n";
        for (int i = 0; i < (int)distinct.size(); ++i)
            printMatrix(std::cout, distinct[i], i + 1, input.p, input.target_sum);

        // 5. Write output file
        writeOutput(outFile, input, input.matrices, distinct);

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}