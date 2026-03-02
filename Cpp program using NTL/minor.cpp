#include <iostream>

int countMinors(int n) 
{
    int totalMinors = 0;
    for (int k = 1; k <= n; ++k) {
        int count = (n - k + 1) * (n - k + 1);
        totalMinors += count;
        std::cout << "Number of " << k << "x" << k << " minors: " << count << std::endl;
    }
    return totalMinors;
}
int main() {
    int n;
    std::cout << "Enter the size of the square matrix (n x n): ";
    std::cin >> n;
    if (n <= 0) {
        std::cout << "Matrix size must be a positive integer." << std::endl;
        return 1;
    }
    int total = countMinors(n);
    std::cout << "\nTotal number of square minors in a " << n << "x" << n << " matrix: " << total << std::endl;
    return 0;
}
