#====================================================
# Third-party paths and shared libraries
#====================================================

set(NTIC_THIRD_PARTY_ROOT "${CMAKE_SOURCE_DIR}/3rdParties")
set(NTIC_THIRD_PARTY_LIB_DIR "${NTIC_THIRD_PARTY_ROOT}/lib")

include(
    ${CMAKE_CURRENT_LIST_DIR}/gklib.cmake
)

include(
    ${CMAKE_CURRENT_LIST_DIR}/metis.cmake
)


#====================================================
# Third-party notices and corresponding source bundle
#====================================================

# Preserve upstream files and ship the integration recipes with the sources.
# The bundle is independently configurable after copying to a writable tree.
set(NTIC_THIRD_PARTY_SOURCE_DEST "share/nticLBM/third-party-sources")

configure_file(
    "${CMAKE_CURRENT_LIST_DIR}/standalone.cmake.in"
    "${CMAKE_BINARY_DIR}/3rdParties/source_bundle/CMakeLists.txt"
    COPYONLY
)

install(FILES "${CMAKE_BINARY_DIR}/3rdParties/source_bundle/CMakeLists.txt"
    DESTINATION "${NTIC_THIRD_PARTY_SOURCE_DEST}"
)
install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/"
    DESTINATION "${NTIC_THIRD_PARTY_SOURCE_DEST}/cmake/3rdparties"
)
install(DIRECTORY "${NTIC_THIRD_PARTY_ROOT}/GKlib/"
    DESTINATION "${NTIC_THIRD_PARTY_SOURCE_DEST}/3rdParties/GKlib"
    PATTERN lib EXCLUDE PATTERN .git EXCLUDE
)
install(DIRECTORY "${NTIC_THIRD_PARTY_ROOT}/METIS/src/"
    DESTINATION "${NTIC_THIRD_PARTY_SOURCE_DEST}/3rdParties/METIS/src"
)
install(FILES "${NTIC_THIRD_PARTY_ROOT}/METIS/SOURCE_VERSION.txt"
    DESTINATION "${NTIC_THIRD_PARTY_SOURCE_DEST}/3rdParties/METIS"
)
install(FILES "${NTIC_THIRD_PARTY_ROOT}/THIRD_PARTY_NOTICES.md"
    DESTINATION "${NTIC_THIRD_PARTY_SOURCE_DEST}/3rdParties"
)
install(DIRECTORY "${NTIC_THIRD_PARTY_ROOT}/licenses/"
    DESTINATION "${NTIC_THIRD_PARTY_SOURCE_DEST}/3rdParties/licenses"
)
