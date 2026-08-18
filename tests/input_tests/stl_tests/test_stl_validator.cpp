#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    //-------------------------------------------------------------------------
    // Command-line arguments
    //-------------------------------------------------------------------------

    if (argc != 4) {
        std::cerr
            << "Usage: "
            << argv[0]
            << " <stl_file> <valid|invalid> <test|full>\n";

        return 1;
    }

    const std::string stlFile  = argv[1];
    const std::string expected = argv[2];
    const std::string mode     = argv[3];

    if (expected != "valid" && expected != "invalid") {
        std::cerr
            << "Error: expected result must be "
            << "\"valid\" or \"invalid\".\n";

        return 1;
    }

    if (mode != "test" &&
        mode != "full") {

        std::cerr
            << "Error: validation mode must be "
            << "\"test\" or \"full\".\n";

        return 1;
    }

    //-------------------------------------------------------------------------
    // Read and validate STL
    //-------------------------------------------------------------------------

    bool validationSucceeded = false;
    std::string errorMessage;

    try {
        auto data = ntic::lbm::stl::read(stlFile);

        ntic::lbm::stl::STLComponents components;

        ntic::lbm::stl::validate(data,mode,components);

        validationSucceeded = true;
    }
    catch (const std::exception& e) {
        validationSucceeded = false;
        errorMessage = e.what();
    }


    //-------------------------------------------------------------------------
    // Compare validation result with the expected result
    //-------------------------------------------------------------------------

    const bool expectedValid =
        (expected == "valid");

    if (validationSucceeded == expectedValid) {

        std::cout
            << "[PASS] "
            << stlFile
            << '\n';

        if (!validationSucceeded) {
            std::cout
                << "       Detected: "
                << errorMessage
                << '\n';
        }

        return 0;
    }


    //-------------------------------------------------------------------------
    // Test failure
    //-------------------------------------------------------------------------

    std::cerr
        << "[FAIL] "
        << stlFile
        << '\n';

    if (expectedValid) {

        std::cerr
            << "       Expected: valid\n"
            << "       Detected: invalid\n"
            << "       Error: "
            << errorMessage
            << '\n';
    }
    else {

        std::cerr
            << "       Expected: invalid\n"
            << "       Detected: valid\n";
    }

    return 1;
}
