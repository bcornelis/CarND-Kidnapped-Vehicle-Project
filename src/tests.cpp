#include "tests.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ParticleFilterTest);

LandmarkObs createLandMark(int id, double x, double y) {
  LandmarkObs newLandmark;
  newLandmark.x = x;
  newLandmark.y = y;
  newLandmark.id = id;
  return newLandmark;
}

Map::single_landmark_s createLandmarkS(int id, double x, double y) {
  Map::single_landmark_s newLandMark;
  newLandMark.id_i = id;
  newLandMark.x_f = x;
  newLandMark.y_f = y;
  return newLandMark;
}

Particle createParticle(double x, double y, double theta) {
  Particle particle;
  particle.x = x;
  particle.y = y;
  particle.theta = theta;
  return particle;
}

void ParticleFilterTest::testInit() {
  // perform the initialization
  double std[3] = { .5, .7, .9 };
  mParticleFilter->init(1000, 10., 10., M_PI, std);

  // validate
  CPPUNIT_ASSERT(mParticleFilter->num_particles == 1000);
  CPPUNIT_ASSERT(mParticleFilter->particles.size() == 1000);
  CPPUNIT_ASSERT(mParticleFilter->weights.size() == 1000);

  // make sure the values of the particle are initialized
  CPPUNIT_ASSERT(mParticleFilter->particles[100].x > 0.);
  CPPUNIT_ASSERT(mParticleFilter->particles[100].y > 0.);
  CPPUNIT_ASSERT(mParticleFilter->particles[100].theta > 0.);

  // and 'verify randomness'
  CPPUNIT_ASSERT(
      mParticleFilter->particles[100].x != mParticleFilter->particles[150].x);
  CPPUNIT_ASSERT(
      mParticleFilter->particles[100].y != mParticleFilter->particles[150].y);
  CPPUNIT_ASSERT(
      mParticleFilter->particles[100].theta
          != mParticleFilter->particles[150].theta);

  // should be marked as initialized
  CPPUNIT_ASSERT(mParticleFilter->is_initialized);
}

void ParticleFilterTest::testPrediction() {
  // create random particles
  Particle p1;
  p1.x = 5.;
  p1.y = 6.;
  p1.theta = M_PI;
  mParticleFilter->particles.push_back(p1);

  // perform the prediction step
  mParticleFilter->prediction(.1, NULL, 4., 5.);

  // and verify
  CPPUNIT_ASSERT(mParticleFilter->particles.size() == 1);
  Particle updated = mParticleFilter->particles[0];
  CPPUNIT_ASSERT(updated.x == 4.6164595691166372);
  CPPUNIT_ASSERT(updated.y == 5.9020660495122979);
  CPPUNIT_ASSERT(updated.theta == 3.6415926535897931);
}

void ParticleFilterTest::testDataAssociation() {
  // create 4 sample landmarks
  std::vector<LandmarkObs> predicted;
  predicted.push_back(createLandMark(0, 5, 5));
  predicted.push_back(createLandMark(1, 10, 5));
  predicted.push_back(createLandMark(2, 10, 10));
  predicted.push_back(createLandMark(3, 5, 10));

  // and some observations
  std::vector<LandmarkObs> observations;
  observations.push_back(createLandMark(-1, 0, 0));
  observations.push_back(createLandMark(-1, 4, 4));
  observations.push_back(createLandMark(-1, 9, 9));
  observations.push_back(createLandMark(-1, 10, 10));

  // perform the data association step
  mParticleFilter->dataAssociation(predicted, observations);

  // and verify the index updates
  CPPUNIT_ASSERT(observations[0].id == 0);
  CPPUNIT_ASSERT(observations[1].id == 0);
  CPPUNIT_ASSERT(observations[2].id == 2);
  CPPUNIT_ASSERT(observations[3].id == 2);
}

void ParticleFilterTest::testFilterOnSensorRange() {
  double sensor_range = 3;

  Particle particle;
  particle.x = 5;
  particle.y = 5;
  Map map_landmarks;
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 0, 0));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 1, 1));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 2, 2));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 3, 3));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 4, 4));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 5, 5));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 6, 6));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 7, 7));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 8, 8));
  map_landmarks.landmark_list.push_back(createLandmarkS(0, 9, 9));

  std::vector<LandmarkObs> filteredLandmarks;
  filterOnSensorRange(particle, map_landmarks, sensor_range, filteredLandmarks);

  CPPUNIT_ASSERT(filteredLandmarks.size() == 5);
}

void ParticleFilterTest::testConvertToMapCoordinates() {

  Particle particle = createParticle(5, 5, M_PI/4);

  std::vector<LandmarkObs> observations;
  observations.push_back(createLandMark(-1, 0, 0));
  observations.push_back(createLandMark(-1, .5, .5));
  observations.push_back(createLandMark(-1, -.9, -.2));

  std::vector<LandmarkObs> converted;
  convertToMapCoordinates(observations, particle, converted);

  // verify landmark 1
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0000, converted[0].x, .0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0000, converted[0].y, .0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0000, converted[1].x, .0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.7071, converted[1].y, .0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5050, converted[2].x, .0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.2221, converted[2].y, .0001);
}

void ParticleFilterTest::testUpdateWeights() {
  //mParticleFilter->updateWeights(sensor_range, std_landmark, observations, map_landmarks);
}

void ParticleFilterTest::setUp() {
  mParticleFilter = new ParticleFilter();
}

void ParticleFilterTest::tearDown() {
  delete mParticleFilter;
}
