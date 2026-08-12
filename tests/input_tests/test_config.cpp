#include "config.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace ntic::lbm::config;

namespace
{

bool
near(
    double a,
    double b,
    double tolerance = 1.0e-12)
{
    return std::abs(a - b)
        <= tolerance
        * std::max(
            1.0,
            std::max(
                std::abs(a),
                std::abs(b)));
}

void
testEmptyConfig()
{
    Config config;

    config.parse("");

    NTIC_EXPECT_EQ(
        config.builtinCount(),
        static_cast<std::size_t>(0));

    NTIC_EXPECT_EQ(
        config.userDefinedCount(),
        static_cast<std::size_t>(0));

    NTIC_EXPECT_EQ(
        config.userDefinedBoundCount(),
        static_cast<std::size_t>(0));

    NTIC_EXPECT_TRUE(
        !config.hasUnboundUserDefinedParameters());
}

void
testBuiltinParsingAndBinding()
{
    Config config;

    config.parse(
        R"(
# General
case name = cavity_case
precision = double
geometry file = "mesh.stl"
restart file = "restart.dat"

# Physical
rho = 1000 kg/m^3
viscosity = 1.0e-3 Pa*s
nu = 1.0e-6 m^2/s
reference length = 10 cm
reference velocity = 2 m/s
Reynolds Number = 20000
total physical time = 5 s
restart time = 1 s
dt = 0.01 s
dx = 2 mm
periodic x = true
periodic y = false
periodic z = yes
symmetry x = false
symmetric y = true
z symmetry = no

# Lattice
velocity set = D3Q19
collision operator = cumulant
tau = 0.8
omega = 1.25
Ma = 0.1
levels = 3
wall boundary model = bouzidi
inlet = zou_he
outlet = convective
nsteps = 5000

# Output
output dir = "results"
write interval = 0.1 s
restart interval = 1 s
)");

    NTIC_EXPECT_EQ(
        config.builtinCount(),
        static_cast<std::size_t>(33));

    NTIC_EXPECT_EQ(
        config.userDefinedCount(),
        static_cast<std::size_t>(0));

    NTIC_EXPECT_TRUE(
        config.contains("case_name"));

    NTIC_EXPECT_TRUE(
        config.contains("Case Name"));

    NTIC_EXPECT_TRUE(
        config.containsBuiltin("MU"));

    NTIC_EXPECT_TRUE(
        !config.containsUserDefined("MU"));

    std::string caseName;
    config.bindString(
        caseName,
        "case_name");

    NTIC_EXPECT_EQ(
        caseName,
        std::string("cavity_case"));

    std::string precision;
    config.bindString(
        precision,
        "Floating Point Precision");

    NTIC_EXPECT_EQ(
        precision,
        std::string("double"));

    char geometryFile[64] = {};
    config.bindString(
        geometryFile,
        sizeof(geometryFile),
        "geometry_file");

    NTIC_EXPECT_EQ(
        std::string(geometryFile),
        std::string("mesh.stl"));

    double density = 0.0;
    config.bindQuantity(
        density,
        "rho");

    NTIC_EXPECT_TRUE(
        near(
            density,
            1000.0));

    double dynamicViscosity = 0.0;
    config.bindQuantity(
        dynamicViscosity,
        "Dynamic Viscosity");

    NTIC_EXPECT_TRUE(
        near(
            dynamicViscosity,
            1.0e-3));

    float gridSpacing = 0.0f;
    config.bindQuantity(
        gridSpacing,
        "grid resolution");

    NTIC_EXPECT_TRUE(
        near(
            static_cast<double>(gridSpacing),
            2.0e-3,
            1.0e-6));

    int levels = 0;
    config.bindQuantity(
        levels,
        "number_of_levels");

    NTIC_EXPECT_EQ(
        levels,
        3);

    int maxSteps = 0;
    config.bindQuantity(
        maxSteps,
        "max steps");

    NTIC_EXPECT_EQ(
        maxSteps,
        5000);

    bool periodicX = false;
    config.bindBool(
        periodicX,
        "x periodic");

    NTIC_EXPECT_TRUE(
        periodicX);

    bool periodicY = true;
    config.bindBool(
        periodicY,
        "periodic_y");

    NTIC_EXPECT_TRUE(
        !periodicY);

    bool symmetricY = false;
    config.bindBool(
        symmetricY,
        "symmetry y");

    NTIC_EXPECT_TRUE(
        symmetricY);

    std::string latticeModel;
    config.bindString(
        latticeModel,
        "velocity set");

    NTIC_EXPECT_EQ(
        latticeModel,
        std::string("D3Q19"));

    double outputInterval = 0.0;
    config.bindQuantity(
        outputInterval,
        "write interval");

    NTIC_EXPECT_TRUE(
        near(
            outputInterval,
            0.1));
}

void
testUserDefinedBinding()
{
    Config config;

    config.parse(
        R"(
pressure = 101325 Pa
surface tension = 72 mN/m
custom flag = true
thermal model = enthalpy
iterations = 12
fractional iterations = 12.5
custom path = "/tmp/result.dat"
)");

    NTIC_EXPECT_EQ(
        config.builtinCount(),
        static_cast<std::size_t>(0));

    NTIC_EXPECT_EQ(
        config.userDefinedCount(),
        static_cast<std::size_t>(7));

    NTIC_EXPECT_EQ(
        config.userDefinedBoundCount(),
        static_cast<std::size_t>(0));

    NTIC_EXPECT_TRUE(
        config.hasUnboundUserDefinedParameters());

    NTIC_EXPECT_TRUE(
        config.containsUserDefined(
            "Pressure"));

    NTIC_EXPECT_TRUE(
        config.containsUserDefined(
            "surface-tension"));

    NTIC_EXPECT_TRUE(
        config.contains(
            "thermal_model"));

    double pressure = 0.0;

    config.bindUserDefinedQuantity(
        pressure,
        "PRESSURE",
        PhysicalDimension::Pressure());

    NTIC_EXPECT_TRUE(
        near(
            pressure,
            101325.0));

    float surfaceTension = 0.0f;

    config.bindUserDefinedQuantity(
        surfaceTension,
        "Surface-Tension",
        PhysicalDimension(
            ntic::lbm::common::Rational(0),
            ntic::lbm::common::Rational(1),
            ntic::lbm::common::Rational(-2),
            ntic::lbm::common::Rational(0)));

    NTIC_EXPECT_TRUE(
        near(
            static_cast<double>(surfaceTension),
            0.072,
            1.0e-6));

    bool customFlag = false;

    config.bindUserDefinedBool(
        customFlag,
        "custom_flag");

    NTIC_EXPECT_TRUE(
        customFlag);

    std::string thermalModel;

    config.bindUserDefinedString(
        thermalModel,
        "Thermal Model");

    NTIC_EXPECT_EQ(
        thermalModel,
        std::string("enthalpy"));

    int iterations = 0;

    config.bindUserDefinedQuantity(
        iterations,
        "iterations",
        PhysicalDimension::Dimensionless());

    NTIC_EXPECT_EQ(
        iterations,
        12);

    char customPath[64] = {};

    config.bindUserDefinedString(
        customPath,
        sizeof(customPath),
        "custom path");

    NTIC_EXPECT_EQ(
        std::string(customPath),
        std::string("/tmp/result.dat"));

    NTIC_EXPECT_EQ(
        config.userDefinedBoundCount(),
        static_cast<std::size_t>(6));

    NTIC_EXPECT_TRUE(
        config.hasUnboundUserDefinedParameters());

    std::ostringstream report;

    config.reportUnboundUserDefinedParameters(
        report);

    const std::string reportText =
        report.str();

    NTIC_EXPECT_TRUE(
        reportText.find(
            "fractional_iterations")
        != std::string::npos);

    NTIC_EXPECT_TRUE(
        reportText.find(
            "line")
        != std::string::npos);

    NTIC_EXPECT_TRUE(
        reportText.find(
            "Quantity")
        != std::string::npos);

    double fractionalIterations = 0.0;

    config.bindUserDefinedQuantity(
        fractionalIterations,
        "fractional iterations",
        PhysicalDimension::Dimensionless());

    NTIC_EXPECT_TRUE(
        near(
            fractionalIterations,
            12.5));

    NTIC_EXPECT_EQ(
        config.userDefinedBoundCount(),
        static_cast<std::size_t>(7));

    NTIC_EXPECT_TRUE(
        !config.hasUnboundUserDefinedParameters());

    std::ostringstream emptyReport;

    config.reportUnboundUserDefinedParameters(
        emptyReport);

    NTIC_EXPECT_TRUE(
        emptyReport.str().empty());
}

void
testCommentsAndWhitespace()
{
    Config config;

    config.parse(
        R"(
        # full-line comment

        case name     =     "demo # case"     # trailing comment
        density       :     998 kg/m^3         # density
        periodic x    =     ON                 # bool
        custom value  =     5 cm               # user value

        )");

    std::string caseName;
    config.bindString(
        caseName,
        "case_name");

    NTIC_EXPECT_EQ(
        caseName,
        std::string("demo # case"));

    double density = 0.0;
    config.bindQuantity(
        density,
        "density");

    NTIC_EXPECT_TRUE(
        near(
            density,
            998.0));

    bool periodicX = false;
    config.bindBool(
        periodicX,
        "periodic_x");

    NTIC_EXPECT_TRUE(
        periodicX);

    double customValue = 0.0;

    config.bindUserDefinedQuantity(
        customValue,
        "CUSTOM-VALUE",
        PhysicalDimension::Length());

    NTIC_EXPECT_TRUE(
        near(
            customValue,
            0.05));
}

void
testTypeErrors()
{
    {
        Config config;

        NTIC_EXPECT_THROW(
            config.parse(
                "density = true"),
            std::runtime_error);
    }

    {
        Config config;

        NTIC_EXPECT_THROW(
            config.parse(
                "periodic_x = 1"),
            std::runtime_error);
    }

    {
        Config config;

        NTIC_EXPECT_THROW(
            config.parse(
                "case_name = true"),
            std::runtime_error);
    }

    {
        Config config;

        config.parse(
            "custom = true");

        double value = 0.0;

        NTIC_EXPECT_THROW(
            config.bindUserDefinedQuantity(
                value,
                "custom",
                PhysicalDimension::Dimensionless()),
            std::runtime_error);

        NTIC_EXPECT_EQ(
            config.userDefinedBoundCount(),
            static_cast<std::size_t>(0));
    }

    {
        Config config;

        config.parse(
            "custom = 1 m");

        bool value = false;

        NTIC_EXPECT_THROW(
            config.bindUserDefinedBool(
                value,
                "custom"),
            std::runtime_error);

        NTIC_EXPECT_EQ(
            config.userDefinedBoundCount(),
            static_cast<std::size_t>(0));
    }

    {
        Config config;

        config.parse(
            "custom = false");

        std::string value;

        NTIC_EXPECT_THROW(
            config.bindUserDefinedString(
                value,
                "custom"),
            std::runtime_error);

        NTIC_EXPECT_EQ(
            config.userDefinedBoundCount(),
            static_cast<std::size_t>(0));
    }

    {
        Config config;

        config.parse(
            "density = 1000 kg/m^3");

        bool value = false;

        NTIC_EXPECT_THROW(
            config.bindBool(
                value,
                "density"),
            std::runtime_error);
    }
}

void
testDimensionErrors()
{
    {
        Config config;

        NTIC_EXPECT_THROW(
            config.parse(
                "density = 10 m"),
            std::runtime_error);
    }

    {
        Config config;

        config.parse(
            "custom = 10 m");

        double value = 0.0;

        NTIC_EXPECT_THROW(
            config.bindUserDefinedQuantity(
                value,
                "custom",
                PhysicalDimension::Time()),
            std::runtime_error);

        NTIC_EXPECT_EQ(
            config.userDefinedBoundCount(),
            static_cast<std::size_t>(0));
    }
}

void
testIntegerConversion()
{
    {
        Config config;

        config.parse(
            "max_steps = 100");

        int value = 0;

        config.bindQuantity(
            value,
            "max_steps");

        NTIC_EXPECT_EQ(
            value,
            100);
    }

    {
        Config config;

        config.parse(
            "max_steps = 100.5");

        int value = 0;

        NTIC_EXPECT_THROW(
            config.bindQuantity(
                value,
                "max_steps"),
            std::runtime_error);
    }

    {
        Config config;

        config.parse(
            "custom = 9.5");

        int value = 0;

        NTIC_EXPECT_THROW(
            config.bindUserDefinedQuantity(
                value,
                "custom",
                PhysicalDimension::Dimensionless()),
            std::runtime_error);

        NTIC_EXPECT_EQ(
            config.userDefinedBoundCount(),
            static_cast<std::size_t>(0));
    }
}

void
testStringBuffers()
{
    {
        Config config;

        config.parse(
            "case_name = abc");

        char buffer[4] = {};

        config.bindString(
            buffer,
            sizeof(buffer),
            "case_name");

        NTIC_EXPECT_EQ(
            std::string(buffer),
            std::string("abc"));
    }

    {
        Config config;

        config.parse(
            "case_name = abcde");

        char buffer[4] = {};

        NTIC_EXPECT_THROW(
            config.bindString(
                buffer,
                sizeof(buffer),
                "case_name"),
            std::runtime_error);
    }

    {
        Config config;

        config.parse(
            "custom = abcde");

        char buffer[4] = {};

        NTIC_EXPECT_THROW(
            config.bindUserDefinedString(
                buffer,
                sizeof(buffer),
                "custom"),
            std::runtime_error);

        NTIC_EXPECT_EQ(
            config.userDefinedBoundCount(),
            static_cast<std::size_t>(0));
    }

    {
        Config config;

        config.parse(
            "custom = abc");

        NTIC_EXPECT_THROW(
            config.bindUserDefinedString(
                nullptr,
                8,
                "custom"),
            std::invalid_argument);
    }
}

void
testDuplicateParameters()
{
    {
        Config config;

        NTIC_EXPECT_THROW(
            config.parse(
                R"(
density = 1000 kg/m^3
rho = 998 kg/m^3
)"),
            std::runtime_error);
    }

    {
        Config config;

        NTIC_EXPECT_THROW(
            config.parse(
                R"(
custom value = 1
custom_value = 2
)"),
            std::runtime_error);
    }
}

void
testMissingParameters()
{
    Config config;

    config.parse(
        "density = 1000 kg/m^3");

    double value = 0.0;

    NTIC_EXPECT_THROW(
        config.bindQuantity(
            value,
            "dynamic_viscosity"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        config.bindUserDefinedQuantity(
            value,
            "unknown_custom",
            PhysicalDimension::Dimensionless()),
        std::runtime_error);
}

void
testInvalidIdentifiers()
{
    Config config;

    NTIC_EXPECT_THROW(
        config.parse(
            "dynamic/viscosity = 1 Pa*s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        config.parse(
            "dynamic.viscosity = 1 Pa*s"),
        std::runtime_error);

    NTIC_EXPECT_THROW(
        config.parse(
            "3d parameter = 1"),
        std::runtime_error);
}

void
testTransactionalParse()
{
    Config config;

    config.parse(
        R"(
density = 1000 kg/m^3
custom = true
)");

    NTIC_EXPECT_EQ(
        config.builtinCount(),
        static_cast<std::size_t>(1));

    NTIC_EXPECT_EQ(
        config.userDefinedCount(),
        static_cast<std::size_t>(1));

    NTIC_EXPECT_THROW(
        config.parse(
            R"(
density = 10 m
another = false
)"),
        std::runtime_error);

    NTIC_EXPECT_EQ(
        config.builtinCount(),
        static_cast<std::size_t>(1));

    NTIC_EXPECT_EQ(
        config.userDefinedCount(),
        static_cast<std::size_t>(1));

    double density = 0.0;

    config.bindQuantity(
        density,
        "density");

    NTIC_EXPECT_TRUE(
        near(
            density,
            1000.0));
}

void
testParseFile()
{
    const std::string fileName =
        "test_config_input.tmp";

    {
        std::ofstream output(
            fileName);

        output
            << "case_name = file_case\n"
            << "density = 950 kg/m^3\n"
            << "custom flag = false\n";
    }

    Config config;

    config.parseFile(
        fileName);

    std::remove(
        fileName.c_str());

    std::string caseName;

    config.bindString(
        caseName,
        "case_name");

    NTIC_EXPECT_EQ(
        caseName,
        std::string("file_case"));

    double density = 0.0;

    config.bindQuantity(
        density,
        "density");

    NTIC_EXPECT_TRUE(
        near(
            density,
            950.0));

    bool customFlag = true;

    config.bindUserDefinedBool(
        customFlag,
        "custom flag");

    NTIC_EXPECT_TRUE(
        !customFlag);

    Config missingFileConfig;

    NTIC_EXPECT_THROW(
        missingFileConfig.parseFile(
            "file_that_should_not_exist.cfg"),
        std::runtime_error);
}

} // namespace

int main()
{
    testEmptyConfig();
    testBuiltinParsingAndBinding();
    testUserDefinedBinding();
    testCommentsAndWhitespace();
    testTypeErrors();
    testDimensionErrors();
    testIntegerConversion();
    testStringBuffers();
    testDuplicateParameters();
    testMissingParameters();
    testInvalidIdentifiers();
    testTransactionalParse();
    testParseFile();

    ntic::lbm::test::printSummary();

    return 0;
}
