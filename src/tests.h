#ifndef TESTS_H
#define TESTS_H

#include "particle_filter.h"
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>

class ParticleFilterTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ParticleFilterTest );
  CPPUNIT_TEST( testInit );
  CPPUNIT_TEST( testPrediction );
  CPPUNIT_TEST( testDataAssociation );
  CPPUNIT_TEST( testFilterOnSensorRange );
  CPPUNIT_TEST( testConvertToMapCoordinates );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void testInit();
  void testPrediction();
  void testDataAssociation();
  void testUpdateWeights();

  void testFilterOnSensorRange();
  void testConvertToMapCoordinates();

private:
  ParticleFilter *mParticleFilter;
};


#endif
