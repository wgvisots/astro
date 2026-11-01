/**********************************************************************
 * test_util.h
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

#pragma once

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "astro/core/vec.h"

namespace astro::test {

inline std::string dataPath(const std::string& name) {
    return std::string(ASTRO_TEST_DATA_DIR) + "/" + name;
}

inline void expectVecNear(const Vec3& a, const Vec3& b, double tol, const char* what = "") {
    EXPECT_NEAR(a.x, b.x, tol) << what << " x";
    EXPECT_NEAR(a.y, b.y, tol) << what << " y";
    EXPECT_NEAR(a.z, b.z, tol) << what << " z";
}

inline void expectMatNear(const Mat3& a, const Mat3& b, double tol, const char* what = "") {
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) EXPECT_NEAR(a(i, j), b(i, j), tol) << what << " (" << i << "," << j << ")";
}

inline void expectOrthonormal(const Mat3& a, double tol = 1e-14) {
    expectMatNear(a * a.transpose(), Mat3::identity(), tol, "A*A^T");
    EXPECT_NEAR(a.determinant(), 1.0, tol);
}

/// Reads a CSV file, skipping blank lines and lines starting with '#'.
/// The first non-comment line is a header and is returned separately.
struct Csv {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;

    [[nodiscard]] std::size_t col(const std::string& name) const {
        for (std::size_t i = 0; i < header.size(); ++i)
            if (header[i] == name) return i;
        throw std::runtime_error("no CSV column " + name);
    }
};

inline Csv readCsv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    Csv csv;
    std::string line;
    bool haveHeader = false;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string f;
        while (std::getline(ss, f, ',')) fields.push_back(f);
        if (!haveHeader) {
            csv.header = std::move(fields);
            haveHeader = true;
        } else {
            csv.rows.push_back(std::move(fields));
        }
    }
    return csv;
}

} // namespace astro::test
