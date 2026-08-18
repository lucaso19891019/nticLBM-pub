#include "stl_reader.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ntic::lbm::stl {

namespace {

//=============================================================================
// Constants
//=============================================================================

constexpr std::uint64_t BINARY_HEADER_SIZE = 80;
constexpr std::uint64_t BINARY_COUNT_SIZE  = 4;
constexpr std::uint64_t BINARY_FACET_SIZE  = 50;
constexpr std::uint64_t BINARY_BASE_SIZE =
    BINARY_HEADER_SIZE + BINARY_COUNT_SIZE;


//=============================================================================
// Utility
//=============================================================================

[[nodiscard]]
bool isFinite(const STLVector& v) noexcept
{
    return std::isfinite(v[0]) &&
           std::isfinite(v[1]) &&
           std::isfinite(v[2]);
}


[[noreturn]]
void throwError(
    const std::filesystem::path& path,
    const std::string& message)
{
    throw std::runtime_error(
        "Invalid STL file '" + path.string() + "': " + message);
}


//=============================================================================
// Little-endian binary readers
//=============================================================================

[[nodiscard]]
std::uint16_t readUInt16LE(std::istream& input)
{
    unsigned char bytes[2];

    input.read(
        reinterpret_cast<char*>(bytes),
        sizeof(bytes));

    if (!input) {
        throw std::runtime_error(
            "Unexpected end of binary STL file.");
    }

    return
        static_cast<std::uint16_t>(bytes[0]) |
        (static_cast<std::uint16_t>(bytes[1]) << 8);
}


[[nodiscard]]
std::uint32_t readUInt32LE(std::istream& input)
{
    unsigned char bytes[4];

    input.read(
        reinterpret_cast<char*>(bytes),
        sizeof(bytes));

    if (!input) {
        throw std::runtime_error(
            "Unexpected end of binary STL file.");
    }

    return
        static_cast<std::uint32_t>(bytes[0]) |
        (static_cast<std::uint32_t>(bytes[1]) << 8) |
        (static_cast<std::uint32_t>(bytes[2]) << 16) |
        (static_cast<std::uint32_t>(bytes[3]) << 24);
}


[[nodiscard]]
float readFloat32LE(std::istream& input)
{
    const std::uint32_t bits = readUInt32LE(input);

    float value;

    static_assert(
        sizeof(float) == sizeof(std::uint32_t),
        "STL reader requires 32-bit float.");

    std::memcpy(&value, &bits, sizeof(value));

    return value;
}


[[nodiscard]]
STLVector readBinaryVector(std::istream& input)
{
    return {
        static_cast<double>(readFloat32LE(input)),
        static_cast<double>(readFloat32LE(input)),
        static_cast<double>(readFloat32LE(input))
    };
}


//=============================================================================
// Binary STL detection
//=============================================================================

[[nodiscard]]
bool isBinarySTL(
    const std::filesystem::path& path,
    const std::uint64_t fileSize)
{
    if (fileSize < BINARY_BASE_SIZE) {
        return false;
    }

    std::ifstream input(path, std::ios::binary);

    if (!input) {
        throw std::runtime_error(
            "Unable to open STL file '" + path.string() + "'.");
    }

    input.seekg(
        static_cast<std::streamoff>(BINARY_HEADER_SIZE),
        std::ios::beg);

    if (!input) {
        return false;
    }

    const std::uint32_t facetCount = readUInt32LE(input);

    if (facetCount >
        (std::numeric_limits<std::uint64_t>::max() -
         BINARY_BASE_SIZE) /
            BINARY_FACET_SIZE) {
        return false;
    }

    const std::uint64_t expectedSize =
        BINARY_BASE_SIZE +
        static_cast<std::uint64_t>(facetCount) *
            BINARY_FACET_SIZE;

    return expectedSize == fileSize;
}


//=============================================================================
// Binary STL
//=============================================================================

[[nodiscard]]
STLData readBinary(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);

    if (!input) {
        throw std::runtime_error(
            "Unable to open STL file '" + path.string() + "'.");
    }

    // Skip 80-byte binary STL header.
    input.seekg(
        static_cast<std::streamoff>(BINARY_HEADER_SIZE),
        std::ios::beg);

    if (!input) {
        throwError(path, "unable to read binary header.");
    }

    std::uint32_t facetCount;

    try {
        facetCount = readUInt32LE(input);
    }
    catch (const std::exception&) {
        throwError(path, "unable to read binary facet count.");
    }

    if (facetCount == 0) {
        throwError(path, "binary STL contains no facets.");
    }

    STLData data;
    data.format = STLFormat::Binary;
    data.facets.reserve(facetCount);

    for (std::uint32_t i = 0; i < facetCount; ++i) {

        STLFacet facet;

        try {
            facet.normal = readBinaryVector(input);

            for (std::size_t v = 0; v < 3; ++v) {
                facet.vertices[v] =
                    readBinaryVector(input);
            }

            // Attribute byte count.
            static_cast<void>(readUInt16LE(input));
        }
        catch (const std::exception&) {
            throwError(
                path,
                "unexpected end of file while reading facet " +
                std::to_string(i) + ".");
        }

        if (!isFinite(facet.normal)) {
            throwError(
                path,
                "non-finite normal in facet " +
                std::to_string(i) + ".");
        }

        for (std::size_t v = 0; v < 3; ++v) {
            if (!isFinite(facet.vertices[v])) {
                throwError(
                    path,
                    "non-finite vertex in facet " +
                    std::to_string(i) + ".");
            }
        }

        data.facets.push_back(facet);
    }

    return data;
}


//=============================================================================
// ASCII parsing helpers
//=============================================================================

void expect(
    std::istream& input,
    const std::filesystem::path& path,
    const std::string& expected,
    const std::size_t facetIndex)
{
    std::string token;

    if (!(input >> token)) {
        throwError(
            path,
            "unexpected end of file in facet " +
            std::to_string(facetIndex) +
            "; expected '" + expected + "'.");
    }

    if (token != expected) {
        throwError(
            path,
            "expected '" + expected +
            "' in facet " +
            std::to_string(facetIndex) +
            ", but found '" + token + "'.");
    }
}


[[nodiscard]]
STLVector readASCIIVector(
    std::istream& input,
    const std::filesystem::path& path,
    const std::size_t facetIndex,
    const char* description)
{
    STLVector value{};

    if (!(input >> value[0] >> value[1] >> value[2])) {
        throwError(
            path,
            "unable to read " +
            std::string(description) +
            " in facet " +
            std::to_string(facetIndex) + ".");
    }

    if (!isFinite(value)) {
        throwError(
            path,
            "non-finite " +
            std::string(description) +
            " in facet " +
            std::to_string(facetIndex) + ".");
    }

    return value;
}


//=============================================================================
// ASCII STL
//=============================================================================

[[nodiscard]]
STLData readASCII(const std::filesystem::path& path)
{
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "Unable to open STL file '" + path.string() + "'.");
    }

    std::string token;

    if (!(input >> token) || token != "solid") {
        throwError(
            path,
            "ASCII STL must begin with 'solid'.");
    }

    // Consume the remainder of the first line.
    std::string line;
    std::getline(input, line);

    STLData data;
    data.format = STLFormat::ASCII;

    std::size_t facetIndex = 0;

    while (input >> token) {

        if (token == "endsolid") {
            // The optional solid name is ignored.
            std::getline(input, line);

            // Nothing except whitespace is allowed afterwards.
            std::string trailing;

            if (input >> trailing) {
                throwError(
                    path,
                    "unexpected data after 'endsolid'.");
            }

            if (data.facets.empty()) {
                throwError(
                    path,
                    "ASCII STL contains no facets.");
            }

            return data;
        }

        if (token != "facet") {
            throwError(
                path,
                "expected 'facet' or 'endsolid', but found '" +
                token + "'.");
        }

        expect(
            input,
            path,
            "normal",
            facetIndex);

        STLFacet facet;

        facet.normal =
            readASCIIVector(
                input,
                path,
                facetIndex,
                "normal");

        expect(
            input,
            path,
            "outer",
            facetIndex);

        expect(
            input,
            path,
            "loop",
            facetIndex);

        for (std::size_t v = 0; v < 3; ++v) {

            expect(
                input,
                path,
                "vertex",
                facetIndex);

            facet.vertices[v] =
                readASCIIVector(
                    input,
                    path,
                    facetIndex,
                    "vertex");
        }

        expect(
            input,
            path,
            "endloop",
            facetIndex);

        expect(
            input,
            path,
            "endfacet",
            facetIndex);

        data.facets.push_back(facet);

        ++facetIndex;
    }

    throwError(
        path,
        "ASCII STL is missing 'endsolid'.");
}


//=============================================================================
// File size
//=============================================================================

[[nodiscard]]
std::uint64_t getFileSize(
    const std::filesystem::path& path)
{
    std::error_code error;

    const auto size =
        std::filesystem::file_size(path, error);

    if (error) {
        throw std::runtime_error(
            "Unable to determine size of STL file '" +
            path.string() + "'.");
    }

    return static_cast<std::uint64_t>(size);
}

} // namespace


//=============================================================================
// Public STL reader
//=============================================================================

STLData read(const std::filesystem::path& path)
{
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error(
            "STL file does not exist: '" +
            path.string() + "'.");
    }

    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error(
            "STL path is not a regular file: '" +
            path.string() + "'.");
    }

    const std::uint64_t fileSize =
        getFileSize(path);

    if (fileSize == 0) {
        throwError(path, "file is empty.");
    }

    if (isBinarySTL(path, fileSize)) {
        return readBinary(path);
    }

    return readASCII(path);
}

} // namespace ntic::lbm::input::stl