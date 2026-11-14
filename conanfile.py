###########################################################################
# conanfile.py
#
# Copyright 2024 W. Greg Visotski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###########################################################################
#
# Conan 2 recipe for astro. From a checkout:
#
#   conan create . --build=missing                        # static library
#   conan create . --build=missing -o "astro/*:shared=True"
#
# and in a consuming project's conanfile: requires = "astro/<version>". The package gives CMake
# `find_package(astro)` and the target `astro::astro`, as an installed astro does. The test in
# test_package builds the same consumer program CI builds against the Debian package
# (tests/package/main.cpp).

import os
import re

from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.files import copy, load, rmdir

required_conan_version = ">=2.0"


class AstroConan(ConanFile):
    name = "astro"
    description = "Time, geodesy, reference frame, ephemeris and orbit utilities"
    license = "Apache-2.0"
    url = "https://github.com/wgvisots/astro"
    homepage = "https://github.com/wgvisots/astro"
    topics = ("astrodynamics", "astronomy", "geodesy", "orbit", "sgp4", "time-scales")
    package_type = "library"

    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True, "boost/*:header_only": True}

    # Only what the library build reads: no tests, examples, studies or test data.
    exports_sources = (
        "CMakeLists.txt",
        "cmake/*",
        "include/*",
        "src/*",
        "docs/Doxyfile.in",
        "LICENSE",
        "THIRD_PARTY_NOTICES.md",
        "README.md",
        "CHANGELOG.md",
    )

    def set_version(self):
        # The version lives in one place: the project() call of CMakeLists.txt.
        cmake = load(self, os.path.join(self.recipe_folder, "CMakeLists.txt"))
        match = re.search(r"project\(\s*astro\s+VERSION\s+(\d+\.\d+\.\d+)", cmake)
        if not match:
            raise ValueError("conanfile.py: no 'project(astro VERSION x.y.z)' in CMakeLists.txt")
        self.version = match.group(1)

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        # Boost.Describe, headers only; its headers appear in astro's public headers.
        self.requires("boost/[>=1.77 <2]", transitive_headers=True)

    def validate(self):
        check_min_cppstd(self, 20)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["ASTRO_BUILD_TESTS"] = False
        tc.cache_variables["ASTRO_BUILD_EXAMPLES"] = False
        tc.generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))
        copy(self, "THIRD_PARTY_NOTICES.md", self.source_folder, os.path.join(self.package_folder, "licenses"))
        CMake(self).install()
        # Consumers get Conan's generated CMake files (CMakeDeps), not astro's own.
        rmdir(self, os.path.join(self.package_folder, "lib", "cmake"))
        rmdir(self, os.path.join(self.package_folder, "share"))

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "astro")
        self.cpp_info.set_property("cmake_target_name", "astro::astro")
        self.cpp_info.libs = ["astro"]
        self.cpp_info.requires = ["boost::headers"]
        if not self.options.shared:
            # A static astro defines this for its consumers (see CMakeLists.txt), so that the
            # export macro in astro/export.h expands to nothing.
            self.cpp_info.defines = ["ASTRO_STATIC_DEFINE"]
        if self.settings.os in ("Linux", "FreeBSD"):
            self.cpp_info.system_libs = ["m"]
