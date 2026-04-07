#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <iomanip>
#include <cmath>
#include <iostream>

// Represents one matrix: its ID and its 2D data
struct Matrix {
    int id;
    std::vector<std::vector<double>> data;
};

// Parse all matrices from the NTL output file
std::vector<Matrix> parseMatrices(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) std::exit(1);

    std::vector<Matrix> matrices;
    std::string line;

    std::regex matrixHeader(R"(Matrix\s+#(\d+):)");
    std::regex matrixRow(R"(\|\s*((?:\d+\s*)+)\|)");

    Matrix current;
    bool inMatrix = false;

    while (std::getline(file, line)) {
        std::smatch m;

        if (std::regex_search(line, m, matrixHeader)) {
            if (inMatrix && !current.data.empty())
                matrices.push_back(current);
            current = Matrix();
            current.id = std::stoi(m[1]);
            current.data.clear();
            inMatrix = true;
            continue;
        }

        if (inMatrix && std::regex_search(line, m, matrixRow)) {
            std::istringstream iss(m[1].str());
            std::vector<double> row;
            double val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) current.data.push_back(row);
        }
    }

    if (inMatrix && !current.data.empty())
        matrices.push_back(current);

    return matrices;
}

// Modular exponentiation: base^exp mod p
long long modpow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

// Modular inverse using Fermat's little theorem (p must be prime)
long long modinv(long long a, long long p) {
    return modpow(((a % p) + p) % p, p - 2, p);
}

// Compute determinant using Gaussian elimination mod p (p must be prime)
// Returns determinant in range [0, p-1], or -1 if matrix is singular mod p
long long determinantModP(std::vector<std::vector<double>> matD, long long p) {
    int n = matD.size();

    // Convert to integer matrix mod p
    std::vector<std::vector<long long>> mat(n, std::vector<long long>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mat[i][j] = (((long long)matD[i][j] % p) + p) % p;

    long long det = 1;

    for (int col = 0; col < n; ++col) {
        // Find pivot row with non-zero entry mod p
        int pivot = -1;
        for (int row = col; row < n; ++row) {
            if (mat[row][col] != 0) { pivot = row; break; }
        }

        if (pivot == -1) return 0; // singular mod p

        if (pivot != col) {
            std::swap(mat[pivot], mat[col]);
            det = (p - det) % p; // multiply by -1 mod p
        }

        det = det * mat[col][col] % p;

        long long inv = modinv(mat[col][col], p);
        for (int row = col + 1; row < n; ++row) {
            long long factor = mat[row][col] * inv % p;
            for (int k = col; k < n; ++k)
                mat[row][k] = (mat[row][k] - factor * mat[col][k] % p + p) % p;
        }
    }

    return det;
}

// Write a matrix in NTL-style format to an output stream
void writeMatrix(std::ofstream& out, const Matrix& mat, long long det, long long p) {
    int maxVal = 0;
    for (const auto& row : mat.data)
        for (double v : row)
            if ((int)v > maxVal) maxVal = (int)v;
    int width = std::to_string(maxVal).size() + 1;

    int cols = mat.data.empty() ? 0 : mat.data[0].size();
    std::string border = "  +-" + std::string(cols * width, '-') + "-+";

    out << "  Matrix #" << mat.id << ":\n";
    out << border << "\n";
    for (const auto& row : mat.data) {
        out << "  |";
        for (double v : row) {
            std::string s = std::to_string((int)v);
            out << std::string(width - s.size(), ' ') << s << " ";
        }
        out << "|\n";
    }
    out << border << "\n";
    out << "  Determinant (mod " << p << ") = " << det << "\n\n";
}

int main(int argc, char* argv[]) {
    std::string inputFile  = "output_pattern_1.txt";
    std::string singFile   = "singular.txt";
    std::string nonSingFile= "nonsingular.txt";
    std::string summaryFile= "summary.txt";
    long long p = 17; // default prime modulus

    if (argc >= 2) inputFile   = argv[1];
    if (argc >= 3) singFile    = argv[2];
    if (argc >= 4) nonSingFile = argv[3];
    if (argc >= 5) summaryFile = argv[4];

    // Ask user for p
    std::cout << "Enter prime modulus p: ";
    std::cin >> p;
    std::cout << "Using p = " << p << "\n";

    std::vector<Matrix> matrices = parseMatrices(inputFile);

    std::ofstream singOut(singFile);
    std::ofstream nonSingOut(nonSingFile);
    std::ofstream sumOut(summaryFile);

    if (!singOut.is_open() || !nonSingOut.is_open() || !sumOut.is_open())
        return 1;

    // Headers
    singOut    << "===========================================\n";
    singOut    << "  SINGULAR MATRICES (det = 0)\n";
    singOut    << "===========================================\n\n";

    nonSingOut << "===========================================\n";
    nonSingOut << "  NON-SINGULAR MATRICES (det != 0)\n";
    nonSingOut << "===========================================\n\n";

    sumOut     << "===========================================\n";
    sumOut     << "  Per-Matrix Classification Report\n";
    sumOut     << "  Input: " << inputFile << "\n";
    sumOut     << "  Prime modulus p = " << p << "\n";
    sumOut     << "===========================================\n\n";
    sumOut     << std::left
               << std::setw(14) << "Matrix ID"
               << std::setw(20) << "Det (mod p)"
               << "Classification\n";
    sumOut     << "-------------------------------------------\n";

    int singCount = 0, nonSingCount = 0;

    for (const auto& mat : matrices) {
        long long det = determinantModP(mat.data, p);
        bool isSingular = (det == 0);

        if (isSingular) {
            writeMatrix(singOut, mat, det, p);
            singCount++;
        } else {
            writeMatrix(nonSingOut, mat, det, p);
            nonSingCount++;
        }

        // Per-matrix line in summary
        sumOut << std::left
               << std::setw(14) << ("#" + std::to_string(mat.id))
               << std::setw(20) << det
               << (isSingular ? "SINGULAR" : "Non-Singular") << "\n";
    }

    // Footers for singular / nonsingular files
    singOut << "-------------------------------------------\n";
    singOut << "Total singular matrices     : " << singCount << "\n";
    singOut << "Total matrices in file      : " << matrices.size() << "\n";

    nonSingOut << "-------------------------------------------\n";
    nonSingOut << "Total non-singular matrices : " << nonSingCount << "\n";
    nonSingOut << "Total matrices in file      : " << matrices.size() << "\n";

    // Summary footer
    sumOut << "-------------------------------------------\n";
    sumOut << "Total matrices parsed       : " << matrices.size() << "\n";
    sumOut << "Non-singular (saved to " << nonSingFile << ") : " << nonSingCount << "\n";
    sumOut << "Singular     (saved to " << singFile    << ") : " << singCount    << "\n";
    sumOut << "===========================================\n";

    singOut.close();
    nonSingOut.close();
    sumOut.close();

    return 0;
}