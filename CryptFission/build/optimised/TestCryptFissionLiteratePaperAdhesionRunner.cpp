/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "projects/CryptFission/test/TestCryptFissionLiteratePaperAdhesion.hpp"

static TestCryptFissionLiteratePaper suite_TestCryptFissionLiteratePaper;

static CxxTest::List Tests_TestCryptFissionLiteratePaper = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestCryptFissionLiteratePaper( "projects/CryptFission/test/TestCryptFissionLiteratePaperAdhesion.hpp", 69, "TestCryptFissionLiteratePaper", suite_TestCryptFissionLiteratePaper, Tests_TestCryptFissionLiteratePaper );

static class TestDescription_TestCryptFissionLiteratePaper_TestEpithelialLayerUndergoingFission : public CxxTest::RealTestDescription {
public:
 TestDescription_TestCryptFissionLiteratePaper_TestEpithelialLayerUndergoingFission() : CxxTest::RealTestDescription( Tests_TestCryptFissionLiteratePaper, suiteDescription_TestCryptFissionLiteratePaper, 72, "TestEpithelialLayerUndergoingFission" ) {}
 void runTest() { suite_TestCryptFissionLiteratePaper.TestEpithelialLayerUndergoingFission(); }
} testDescription_TestCryptFissionLiteratePaper_TestEpithelialLayerUndergoingFission;

#include <cxxtest/Root.cpp>
