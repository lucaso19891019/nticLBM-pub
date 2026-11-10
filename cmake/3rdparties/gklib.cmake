#====================================================
# GKlib shared-library configuration
#====================================================

# Scope generic upstream options so they do not configure nticLBM itself.
function(ntic_add_gklib)
    set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
    set(SHARED ON)
    set(GKLIB_BUILD_APPS OFF)
    set(OPENMP OFF)
    if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$")
        set(NO_X86 ON)
    endif()
    add_subdirectory(
        "${NTIC_THIRD_PARTY_ROOT}/GKlib"
        "${CMAKE_BINARY_DIR}/3rdParties/GKlib"
        EXCLUDE_FROM_ALL
    )
    set_target_properties(
        GKlib
        PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS ON
        LIBRARY_OUTPUT_DIRECTORY "${NTIC_THIRD_PARTY_LIB_DIR}"
        RUNTIME_OUTPUT_DIRECTORY "${NTIC_THIRD_PARTY_LIB_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${NTIC_THIRD_PARTY_LIB_DIR}"
    )
    foreach(NTIC_CONFIG DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
        set_target_properties(
            GKlib
            PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY_${NTIC_CONFIG} "${NTIC_THIRD_PARTY_LIB_DIR}"
            RUNTIME_OUTPUT_DIRECTORY_${NTIC_CONFIG} "${NTIC_THIRD_PARTY_LIB_DIR}"
            ARCHIVE_OUTPUT_DIRECTORY_${NTIC_CONFIG} "${NTIC_THIRD_PARTY_LIB_DIR}"
        )
    endforeach()
    # Upstream standalone install paths use its own project root. Keep its
    # install rules excluded and install the shared target explicitly here.
    install(TARGETS GKlib LIBRARY DESTINATION lib RUNTIME DESTINATION lib)
endfunction()

ntic_add_gklib()
