#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

// Function to compute (a - b) mod p, always positive
ZZ modDiff(const ZZ &a, const ZZ &b, const ZZ &p) {
    ZZ diff = (a - b) % p;
    if (diff < 0) diff += p;
    return diff;
}

int main() {
    string filename;
    cout << "Enter input file name: ";
    cin >> filename;

    ifstream infile(filename);
    if(!infile) {
        cerr << "Error opening file " << filename << "\n";
        return 1;
    }

    vector<vector<ZZ>> A;
    string line;

    // Read matrix from file
    while(getline(infile, line)) {
        if(line.empty()) continue;
        istringstream iss(line);
        ZZ num;
        vector<ZZ> row;
        while(iss >> num)
            row.push_back(num);
        if(!row.empty())
            A.push_back(row);
    }
    infile.close();

    if(A.empty()) {
        cerr << "No data found in file.\n";
        return 1;
    }

    int n = A.size();
    int m = A[0].size();

    ZZ p;
    cout << "Enter prime modulus p: ";
    cin >> p;

    ofstream outfile("pattern_property.txt");
    if(!outfile) {
        cerr << "Error creating output file.\n";
        return 1;
    }

    outfile << "Matrix size: " << n << "x" << m << "\n\n";

    // Row sums
    outfile << "Row sums modulo " << p << ":\n";
    for(int i = 0; i < n; i++) {
        ZZ sum = ZZ(0);
        for(int j = 0; j < m; j++)
            sum += A[i][j];
        sum %= p;
        if(sum < 0) sum += p;
        outfile << "Row " << i+1 << " sum mod p = " << sum << "\n";
    }

    // Row-difference sums with detailed format
    outfile << "\nRow-difference sums modulo " << p << ":\n";
    for(int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            // Print element-wise differences
            outfile << "R" << j+1 << " − R" << i+1 << ": [";
            for(int k = 0; k < m; k++) {
                outfile << A[j][k] << "−" << A[i][k];
                if(k != m-1) outfile << ",";
            }
            outfile << "] = [";

            vector<ZZ> diffMod;
            ZZ sumDiff = ZZ(0);
            for(int k = 0; k < m; k++) {
                ZZ d = modDiff(A[j][k], A[i][k], p);
                diffMod.push_back(d);
                sumDiff += d;
            }
            sumDiff %= p;

            for(int k = 0; k < m; k++) {
                outfile << diffMod[k];
                if(k != m-1) outfile << ",";
            }
            outfile << "] mod " << p << "\n";

            // Print sum calculation
            outfile << "Sum = ";
            ZZ sumCheck = ZZ(0);
            for(int k = 0; k < m; k++) {
                sumCheck += diffMod[k];
                outfile << diffMod[k];
                if(k != m-1) outfile << "+";
            }
            outfile << "=" << sumCheck << " → " << sumCheck << " mod " << p << " = " << sumDiff << "\n\n";
        }
    }

    outfile.close();
    cout << "Output written to pattern_property.txt in detailed format.\n";

    return 0;
}
