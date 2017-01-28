// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 2016 Vegard Antun
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include "hadamard.h"
#include "fxt/walsh/walshseq.h"
#include "fxt/walsh/walshpal.h"
#include "Eigen/Dense"


void insertRow(Eigen::MatrixXi & A, int *x, const int row );
void insertCol(Eigen::MatrixXi & A, int *x, const int col );
Eigen::MatrixXi createMatrix(const unsigned int N, 
                             void (*transform)(int *, const uint32_t));

// Test of support functions
bool test_insertRowCol(bool verbose=true);

// Test of transforms
bool test_hadamard_sequency_sign_changes(bool verbose=true);
bool test_hadamardOrdinary_16(bool verbose=true);
bool test_hadamardPaley_16(bool verbose=true);
bool test_hadamardPaley_32(bool verbose=true);
bool test_sign_change_WAL( bool verbose=true);
bool test_powDyadic( bool verbose=true);

int main(int argc, char *argv[]) {
    
    const unsigned long nu = 3;
    const unsigned long N = powDyadic(nu);
    
    test_powDyadic();
    //test_sign_change_WAL();
    test_hadamardOrdinary_16();
    test_hadamardPaley_16();
    test_hadamardPaley_32();
    test_hadamard_sequency_sign_changes();
}

bool test_powDyadic( bool verbose) {
    bool success = true;
    unsigned long n = 1;
    for (unsigned long i = 0; i < 5; i++) {
        unsigned long p = powDyadic(i);
        success = success and (p==n);
        n *= 2;
    }

    if (verbose) {
        const char * result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
        std::cout << result <<  ": powDyadic" << std::endl;
    }
    return success;
}

bool test_sign_change_WAL( bool verbose) {
    const char *result;

    bool success = true;

    for (unsigned int N = 4; N < 1023; N *= 2) {
        // Test the number of sign changes
        int *x = new int [N];

        Eigen::MatrixXi X = Eigen::MatrixXi::Zero(N,N);

        for (unsigned int n = 0; n < N; n++) {
            for (unsigned int t = 0; t < N; t++) {
                x[t] = WAL(N, n, t);
            }

            insertRow(X, x, n);
        }

        for (int row = 0; row < N; row++) {
            int nbr_sign_changes = 0;
            for (int i = 0; i < N-1; i++) {
                if (X(row, i) != X(row, i+1)) {
                    nbr_sign_changes++;
                }
            }
            //std::cout << nbr_sign_changes << std::endl;
            success = success and (nbr_sign_changes == row);
        }

        if (verbose) {
            result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
            std::cout << result <<  ": Sequency order sign changes N = " << N  << std::endl;
        }

        // Test of orthogonal transform
        Eigen::MatrixXi I1 = X*X.transpose();
        Eigen::MatrixXi I2 = X.transpose()*X;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    success = (success and (I1(i,j) == 0));
                    success = (success and (I2(i,j) == 0));
                } else {
                    success = (success and (I1(i,j) == N));
                    success = (success and (I2(i,j) == N));
                }
            }
        }

        if (verbose) {
            result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
            std::cout << result <<  ": Orthogonal N = " << N  << std::endl;
        }

        delete [] x;
    }

    return success;
}


Eigen::MatrixXi createMatrix(const unsigned int N, 
                             void (*transform)(int *, const uint32_t)) 
{
        Eigen::MatrixXi A = Eigen::MatrixXi::Zero(N,N);

        int *x = new int [N];
        for (int i = 0; i < N; i++) {
            std::memset(x, 0, sizeof(int)*N);

            x[i] = 1;
            transform(x, N);
            insertCol(A, x, i);
        }
        
        delete [] x;
        return A;
}


/*

Create the sequency ordered Hadamard matrix and testes that the correct number
of sign changes exists. In addition it tests that the resulting matrix is
orthogonal.

This is done for N = 4, 8, ..., 1024.

*/
bool test_hadamard_sequency_sign_changes(bool verbose) {
    const char *result;

    bool success = true;
    for (int N = 4; N < 1025; N *= 2) {
        int *x = new int [N];

        Eigen::MatrixXi X = Eigen::MatrixXi::Zero(N,N);

        for (int i = 0; i < N; i++) {
            std::memset(x, 0, sizeof(int)*N);

            x[i] = 1;
            hadamardSequency<int>(x, N);
            insertCol(X, x, i);
        }

        for (int row = 0; row < N; row++) {
            int nbr_sign_changes = 0;
            for (int i = 0; i < N-1; i++) {
                if (X(row, i) != X(row, i+1)) {
                    nbr_sign_changes++;
                }
            }
            success = success and (nbr_sign_changes == row);
        }

        if (verbose) {
            result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
            std::cout << result <<  ": Sequency order sign changes N = " << N  << std::endl;
        }

        // Test of orthogonal transform
        Eigen::MatrixXi I1 = X*X.transpose();
        Eigen::MatrixXi I2 = X.transpose()*X;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    success = (success and (I1(i,j) == 0));
                    success = (success and (I2(i,j) == 0));
                } else {
                    success = (success and (I1(i,j) == N));
                    success = (success and (I2(i,j) == N));
                }
            }
        }

        if (verbose) {
            result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
            std::cout << result <<  ": Orthogonal N = " << N  << std::endl;
        }

        delete [] x;
    }

    return success;
}

bool test_hadamardOrdinary_16(bool verbose) {
    const char *result;
    const unsigned int N = 16;
    Eigen::MatrixXi mat1(N,N);

    mat1 << 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
     1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,
     1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,
     1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,
     1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,
     1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
     1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,
     1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,
     1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,
     1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,
     1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,
     1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
     1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,
     1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1,
     1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1,
     1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1;

    int x[N];

    Eigen::MatrixXi X = Eigen::MatrixXi::Zero(N,N);

    for (int i = 0; i < N; i++) {
        std::memset(x, 0, sizeof(int)*N);
        x[i] = 1;

        hadamardOrdinary<int>(x, N);
        insertCol(X, x, i);
    }

    Eigen::MatrixXi Y = mat1-X;

    bool success = (Y.norm() < 1e-10);

    if ( verbose ) {
        result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
        std::cout << result <<  ": hadamardOrdinary16" << std::endl;
    }

    return success;
}


bool test_hadamardPaley_16(bool verbose) {
    const char *result;
    const unsigned int N = 16;
    Eigen::MatrixXi mat1(N,N);

    mat1 << 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
     1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,
     1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,
     1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,
     1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,
     1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,
     1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,
     1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1,
     1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,
     1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,
     1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
     1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1,
     1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,
     1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
     1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,
     1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1;

    int x[N];

    Eigen::MatrixXi X = Eigen::MatrixXi::Zero(N,N);

    for (int i = 0; i < N; i++) {
        std::memset(x, 0, sizeof(int)*N);
        x[i] = 1;

        hadamardPaley<int>(x, N);
        insertCol(X, x, i);
    }

    Eigen::MatrixXi Y = mat1-X;

    bool success = (Y.norm() < 1e-10);

    if ( verbose ) {
        result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
        std::cout << result <<  ": hadamardPaley16" << std::endl;
    }

    return success;
}


bool test_hadamardPaley_32(bool verbose) {
    const char *result;
    const unsigned int N = 32;
    Eigen::MatrixXi mat1(N,N);
    mat1 << 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
     1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,
     1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,
     1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,
     1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1,
     1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,
     1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1,
     1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,
     1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1,
     1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,
     1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1,  1,  1, -1, -1,
     1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,
     1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1,
     1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1,
     1,  1, -1, -1, -1, -1,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  1,  1,
     1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,
     1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1,
     1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,
     1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1,  1, -1,  1, -1,
     1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
     1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1,
     1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1,
     1, -1,  1, -1, -1,  1, -1,  1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1,  1, -1,  1, -1, -1,  1, -1,  1,
     1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,
     1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1,
     1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
     1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1,
     1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,
     1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,
     1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,
     1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1;

    int x[N];

    Eigen::MatrixXi X = Eigen::MatrixXi::Zero(N,N);

    for (int i = 0; i < N; i++) {
        std::memset(x, 0, sizeof(int)*N);
        x[i] = 1;

        hadamardPaley<int>(x, N);
        insertCol(X, x, i);
    }

    Eigen::MatrixXi Y = mat1-X;

    bool success = (Y.norm() < 1e-10);

    if ( verbose ) {
        result = (success) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
        std::cout << result <<  ": hadamardPaley32" << std::endl;
    }

    return success;
}

bool test_insertRowCol(bool verbose) {
    unsigned int N = 5;
    const char *result;

    Eigen::MatrixXi X = Eigen::MatrixXi::Zero(N,N);

    // Test insertRow()
    int x[N];
    for (int n = 0; n < N; n++) {
        for (int i = 0; i < N; i++) {
            x[i] = i;
        }

        insertRow(X, x, n);
    }

    bool success1 = true;

    for (int n = 0; n < N; n++) {
        for (int i = 0; i < N; i++) {
            success1 = success1 and ( (X(n,i) - i) == 0 );
        }
    }

    if (verbose) {
        result = (success1) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
        std::cout << result <<  ": insertRow();" << std::endl;
    }

    // Test insertCol()

    for (int n = 0; n < N; n++) {
        for (int i = 0; i < N; i++) {
            x[i] = i;
        }

        insertCol(X, x, n);
    }

    bool success2 = true;

    for (int n = 0; n < N; n++) {
        for (int i = 0; i < N; i++) {
            success2 = success2 and ( (X(i,n) - i) == 0 );
        }
    }

    if (verbose) {
        result = (success2) ? "\e[32mPassed\e[0m" : "\e[31mFailed\e[0m";
        std::cout << result <<  ": insertCol();" << std::endl;
    }

    return (success1 and success2);

}

void insertRow(Eigen::MatrixXi & A, int *x  , const int row ) {
    const int N = A.cols();

    for (int i = 0; i < N; i++) {
        A(row,i) = x[i];
    }
}

void insertCol(Eigen::MatrixXi & A, int *x, const int col ) {
    const int N = A.rows();

    for (int i = 0; i < N; i++) {
        A(i,col) = x[i];
    }
}

