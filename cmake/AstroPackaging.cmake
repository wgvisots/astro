# Binary packages with CPack: a Debian package (.deb) or an RPM package (.rpm) of the installed
# library, headers and CMake package files, for use with find_package(astro).
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release -DASTRO_BUILD_TESTS=OFF -DASTRO_BUILD_EXAMPLES=OFF
#   cmake --build build -j
#   cpack --config build/CPackConfig.cmake -G DEB      # or -G RPM (needs rpmbuild)
#
# The package installs under /usr. It holds a development package's contents: the static library
# (or the shared one, with -DBUILD_SHARED_LIBS=ON), the headers and the CMake package files. Its only
# dependency is the Boost headers (Boost.Describe), which the headers include.

set(CPACK_PACKAGE_NAME "astro")
set(CPACK_PACKAGE_VENDOR "${ASTRO_PACKAGE_VENDOR}")
set(CPACK_PACKAGE_CONTACT "${ASTRO_PACKAGE_CONTACT}")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_DESCRIPTION
    "C++20 library for time scales, geodesy, reference frames, ephemerides and orbits: IAU 2006/2000A \
frames, SGP4, numerical and analytic propagation with covariance, Lambert's problem, conjunction \
assessment and impact prediction.")
set(CPACK_PACKAGE_HOMEPAGE_URL "${ASTRO_PACKAGE_HOMEPAGE}")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
set(CPACK_STRIP_FILES ON)
set(CPACK_SOURCE_GENERATOR "")
if(NOT CPACK_GENERATOR)
    set(CPACK_GENERATOR "DEB")
endif()

# Debian: dpkg-shlibdeps adds the C++ run-time libraries a shared build needs.
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_PACKAGE_SECTION "libdevel")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-dev (>= 1.77)")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# RPM: rpmbuild finds the run-time libraries by itself.
set(CPACK_RPM_FILE_NAME RPM-DEFAULT)
set(CPACK_RPM_PACKAGE_LICENSE "Apache-2.0")
set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
set(CPACK_RPM_PACKAGE_REQUIRES "boost-devel >= 1.77")
set(CPACK_RPM_PACKAGE_URL "${ASTRO_PACKAGE_HOMEPAGE}")
# Directories that every system already owns must not be claimed by the package.
set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST_ADDITION /usr/lib/cmake /usr/share/doc)

include(CPack)
