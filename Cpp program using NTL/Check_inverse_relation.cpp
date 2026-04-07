#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace NTL;

// Fast matrix exponentiation
mat_ZZ_p matrixPower(mat_ZZ_p A, long k)
{
    long n = A.NumRows();
    mat_ZZ_p result;
    ident(result, n);   // Identity matrix

    while (k > 0)
    {
        if (k % 2 == 1)
            result = result * A;
        A = A * A;
        k /= 2;
    }
    return result;
}

int main()
{
    long n, p, K;

    cout << "Enter prime p: ";
    cin >> p;
    ZZ_p::init(ZZ(p));

    cout << "Enter size of pattern matrix n: ";
    cin >> n;

    mat_ZZ_p A;
    A.SetDims(n, n);

    cout << "Enter pattern matrix elements (row-wise):\n";
    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++)
        {
            long x;
            cin >> x;
            A[i][j] = ZZ_p(x);
        }

    cout << "Enter value of K: ";
    cin >> K;

    ofstream fout("inverse_square.txt");
    if (!fout)
    {
        cout << "Unable to create output file.\n";
        return 0;
    }

    fout << "Prime p = " << p << endl;
    fout << "Matrix size n = " << n << endl;
    fout << "Maximum K = " << K << endl << endl;

    fout << "Matrix A:\n" << A << endl;

    ZZ_p detA = determinant(A);
    if (IsZero(detA))
    {
        fout << "Matrix A is NOT invertible modulo " << p << endl;
        cout << "Matrix A is not invertible. Check output.txt\n";
        fout.close();
        return 0;
    }

    mat_ZZ_p Ainv;
    inv(Ainv, A);
    fout << "Inverse of A (A^-1):\n" << Ainv << endl;


    fout << "==============================" << endl;
    fout << "Powers of (A^-1)^i  for i = 1 to K" << endl;
    fout << "==============================" << endl;

    for (long i = 1; i <= K; i++)
    {
        mat_ZZ_p temp = matrixPower(Ainv, i);
        fout << "(A^-1)^" << i << ":\n";
        fout << temp << endl;
    }

    fout << "==============================" << endl;
    fout << "Inverses of (A^i) for i = 1 to K" << endl;
    fout << "==============================" << endl;

    for (long i = 1; i <= K; i++)
    {
        mat_ZZ_p Ai = matrixPower(A, i);

        ZZ_p detAi = determinant(Ai);
        if (IsZero(detAi))
        {
            fout << "(A^" << i << ") is NOT invertible modulo " << p << endl << endl;
            continue;
        }

        mat_ZZ_p invAi;
        inv(invAi, Ai);

        fout << "(A^" << i << ")^-1:\n";
        fout << invAi << endl;
    }

    fout << "==============================" << endl;
    fout << "Verification for each power i" << endl;
    fout << "==============================" << endl;

    for (long i = 1; i <= K; i++)
    {
        mat_ZZ_p left = matrixPower(Ainv, i);  // (A^-1)^i
        mat_ZZ_p Ai = matrixPower(A, i);       // A^i

        ZZ_p detAi = determinant(Ai);
        if (IsZero(detAi))
        {
            fout << "For i = " << i << " : A^" << i << " is not invertible, so relation not applicable.\n";
            continue;
        }

        mat_ZZ_p right;
        inv(right, Ai);  // (A^i)^-1

        if (left == right)
            fout << "For i = " << i << " : (A^-1)^" << i << " = (A^" << i << ")^-1  --> VERIFIED\n";
        else
            fout << "For i = " << i << " : (A^-1)^" << i << " != (A^" << i << ")^-1  --> NOT VERIFIED\n";
    }

    fout.close();
    cout << "All matrices and results are stored in output.txt\n";

    return 0;
}
