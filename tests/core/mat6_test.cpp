/**********************************************************************
 * mat6_test.cpp
 *
 * Copyright 2024 W. Greg Visotski
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 *         limitations under the License.
 **********************************************************************/

#include "astro/core/mat6.h"

#include <array>
#include <cmath>
#include <cstddef>

#include "astro/core/rotation.h"
#include "test_util.h"

using namespace astro;

namespace {

using Mat6L = std::array<std::array<long double, 6>, 6>;

/// A matrix with element (i, j) = 10 i + j + offset, so that every element names its own row and column.
Mat6 numbered(double offset = 0.0) {
    Mat6 m{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            m[i][j] = 10.0 * static_cast<double>(i) + static_cast<double>(j) + offset;
        }
    }
    return m;
}

/// a * b, summed by brute force in extended precision.
Mat6L productL(const Mat6L& a, const Mat6L& b) {
    Mat6L c{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            for (std::size_t k = 0; k < 6; ++k) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}

Mat6L widened(const Mat6& m) {
    Mat6L w{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            w[i][j] = m[i][j];
        }
    }
    return w;
}

Mat6L transposedL(const Mat6L& m) {
    Mat6L t{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            t[i][j] = m[j][i];
        }
    }
    return t;
}

/// A symmetric positive-definite covariance with position variances of order 1e4 m^2, velocity variances of
/// order 1e-2 m^2/s^2 and correlations between all six components.
Mat6 covariance() {
    // clang-format off
    const double sd[6]   = {100.0, 80.0, 120.0, 0.1, 0.15, 0.08};
    const double rho[6][6] = {{ 1.0,  0.3, -0.2,  0.5,  0.1, -0.1},
                              { 0.3,  1.0,  0.1, -0.2,  0.4,  0.2},
                              {-0.2,  0.1,  1.0,  0.1, -0.3,  0.6},
                              { 0.5, -0.2,  0.1,  1.0,  0.2, -0.1},
                              { 0.1,  0.4, -0.3,  0.2,  1.0,  0.3},
                              {-0.1,  0.2,  0.6, -0.1,  0.3,  1.0}};
    // clang-format on
    Mat6 p{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            p[i][j] = rho[i][j] * sd[i] * sd[j];
        }
    }
    return p;
}

} // namespace

TEST(Mat6, IdentityTransposeAndSum) {
    const Mat6 a = numbered(), t = transpose(a), s = add(a, identity6());
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_EQ(identity6()[i][j], i == j ? 1.0 : 0.0) << i << j;
            EXPECT_EQ(t[i][j], a[j][i]) << i << j;
            EXPECT_EQ(s[i][j], a[i][j] + (i == j ? 1.0 : 0.0)) << i << j;
        }
    }
    // Both terms of the sum count: (10 i + j) + (10 i + j + 0.5) = 20 i + 2 j + 0.5, so 50.5 at (2, 5).
    EXPECT_EQ(add(a, numbered(0.5))[2][5], 50.5);
}

// Whole numbers below 2^53 are added and multiplied exactly, so the brute-force product in long double is the
// exact product and the comparison is for equality.
TEST(Mat6, ProductOfWholeNumbersIsExact) {
    const Mat6  a = numbered(), b = transpose(numbered(-30.0));
    const Mat6  c = multiply(a, b);
    const Mat6L e = productL(widened(a), widened(b));
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_EQ(static_cast<long double>(c[i][j]), e[i][j]) << i << j;
        }
    }
    // By a permutation matrix, P[k][(k + 1) % 6] = 1: P a moves each row up by one, a P each column right by one.
    Mat6 p{};
    for (std::size_t k = 0; k < 6; ++k) {
        p[k][(k + 1) % 6] = 1.0;
    }
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_EQ(multiply(p, a)[i][j], a[(i + 1) % 6][j]) << i << j;
            EXPECT_EQ(multiply(a, p)[i][j], a[i][(j + 5) % 6]) << i << j;
        }
    }
    EXPECT_EQ(multiply(a, identity6()), a);
    EXPECT_EQ(multiply(identity6(), a), a);
}

// congruence(a, p) is a p a^T, and mat6.h promises it exactly symmetric. Truth: the same product by brute
// force in long double. Tolerance relative to the size of the elements: each is a sum of 36 products, with a
// few rounding errors of 1.1e-16 each; the worst seen was 6.9e-17 of the largest element.
TEST(Mat6, CongruenceIsTheExactlySymmetricProduct) {
    Mat6 a = identity6(); // a state transition matrix for 600 s of free flight, with some coupling added
    for (std::size_t i = 0; i < 3; ++i) {
        a[i][i + 3] = 600.0;
    }
    a[0][1]             = 0.01;
    a[4][2]             = -3e-6;
    a[5][3]             = 2e-3;
    const Mat6  p       = covariance();
    const Mat6  c       = congruence(a, p);
    const Mat6L e       = productL(productL(widened(a), widened(p)), transposedL(widened(a)));
    long double largest = 0.0L;
    for (const auto& row : e) {
        for (const long double x : row) {
            largest = std::max(largest, std::abs(x));
        }
    }
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_EQ(c[i][j], c[j][i]) << i << j; // bit for bit
            EXPECT_LE(static_cast<double>(std::abs(c[i][j] - e[i][j]) / largest), 5e-16) << i << j;
        }
    }
}

TEST(Mat6, BlockPicksTheThreeByThreeBlock) {
    const Mat6 a = numbered();
    for (std::size_t br = 0; br < 2; ++br) {
        for (std::size_t bc = 0; bc < 2; ++bc) {
            const Mat3 b = block(a, br, bc);
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    EXPECT_EQ(b(i, j), 10.0 * static_cast<double>(3 * br + i) + static_cast<double>(3 * bc + j));
                }
            }
        }
    }
}

// rotateCovariance(P, R) is diag(R, R) P diag(R, R)^T. Truth: that product by brute force in long double,
// with the full 6x6 matrices. Tolerance relative to the largest element of each block (position, cross,
// velocity), whose sizes differ by 1e6: a few rounding errors; the worst seen was 1.0e-16. The two
// off-diagonal blocks are exact transposes of each other, as mat6.h says.
TEST(Mat6, RotateCovarianceIsTheBlockDiagonalCongruence) {
    const Mat3 r = R3(0.7) * R1(-0.4) * R2(1.9);
    Mat6       d{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            d[i][j]         = r(i, j);
            d[i + 3][j + 3] = r(i, j);
        }
    }
    const Mat6  p = covariance();
    const Mat6  q = rotateCovariance(p, r);
    const Mat6L e = productL(productL(widened(d), widened(p)), transposedL(widened(d)));
    for (std::size_t br = 0; br < 2; ++br) {
        for (std::size_t bc = 0; bc < 2; ++bc) {
            long double largest = 0.0L;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    largest = std::max(largest, std::abs(e[3 * br + i][3 * bc + j]));
                }
            }
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    const std::size_t row = 3 * br + i, col = 3 * bc + j;
                    EXPECT_LE(static_cast<double>(std::abs(q[row][col] - e[row][col]) / largest), 5e-16) << row << col;
                }
            }
        }
    }
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(q[i + 3][j], q[j][i + 3]) << i << j;
        }
    }
}
