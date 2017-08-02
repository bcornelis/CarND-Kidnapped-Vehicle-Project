/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

// for generating random numbers
default_random_engine gen;

void ParticleFilter::init(int numParticles, double x, double y, double theta,
    double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // set the number of particles to use
  num_particles = numParticles;

  // create a gaussian distributions for each of the parameters
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_psi(theta, std[2]);

  // initialize the particles
  for (int particleIdx = 0; particleIdx < num_particles; particleIdx++) {
    // create a new 'random' particle
    Particle newParticle;
    newParticle.id = particleIdx;
    newParticle.x = dist_x(gen);
    newParticle.y = dist_y(gen);
    newParticle.theta = dist_psi(gen);
    newParticle.weight = 1.;

    // store the particle in the filter
    particles.push_back(newParticle);

    // and store the weight
    weights.push_back(1.);
  }

  //mark the particle filter as initialized
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
    double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // random noise distributions
  // TODO!

  // update the state of all particles
  for (std::vector<Particle>::iterator iter = particles.begin();
      iter != particles.end(); ++iter) {

    // for easier calculation/reference
    double current_x = (*iter).x;
    double current_y = (*iter).y;
    double current_theta = (*iter).theta;

    // some intermediate calculations
    double theta_yaw_deltat = current_theta + (yaw_rate * delta_t);
    double velocity_over_yaw = velocity / yaw_rate;

    // calculate the new values
    double new_x = current_x
        + velocity_over_yaw * (sin(theta_yaw_deltat) - sin(current_theta));
    double new_y = current_y
        + velocity_over_yaw * (cos(current_theta) - cos(theta_yaw_deltat));
    double new_theta = theta_yaw_deltat;

    // update the particles state
    (*iter).x = new_x;
    (*iter).y = new_y;
    (*iter).theta = new_theta;

  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
    std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.

  // for each observation
  for (std::vector<LandmarkObs>::iterator observation_iter =
      observations.begin(); observation_iter != observations.end();
      ++observation_iter) {

    double closest_dist = -1.;
    LandmarkObs* closest_landmark;

    // find the predicted landmark closest to the observed one
    for (std::vector<LandmarkObs>::iterator predicted_iter = predicted.begin();
        predicted_iter != predicted.end(); ++predicted_iter) {
      // calculate the distance between both of them
      double current_dist = dist((*observation_iter).x, (*observation_iter).y,
          (*predicted_iter).x, (*predicted_iter).y);

      // is it the closest one so far?
      if (closest_dist < 0 || closest_dist > current_dist) {
        // save the current landmark as being the closest one
        closest_dist = current_dist;
        closest_landmark = &(*predicted_iter);
      }
    }

    // if the closest one has been found, copy it's state
    if (closest_dist > -1) {
      (*observation_iter).id = closest_landmark->id;
    }

  }

}

void filterOnSensorRange(const Particle& particle, const Map& map_landmarks,
    double sensor_range, std::vector<LandmarkObs>& filteredLandmarks) {

  for (int landmarkCtr = 0; landmarkCtr < map_landmarks.landmark_list.size();
      ++landmarkCtr) {
    Map::single_landmark_s c = map_landmarks.landmark_list[landmarkCtr];
    if (dist(particle.x, particle.y, c.x_f, c.y_f) <= sensor_range) {
      LandmarkObs newLmObs;
      newLmObs.id = c.id_i;
      newLmObs.x = c.x_f;
      newLmObs.y = c.y_f;
      filteredLandmarks.push_back(newLmObs);
    }
  }
}

void convertToMapCoordinates(const std::vector<LandmarkObs>& observations,
    const Particle& particle, std::vector<LandmarkObs>& converted) {
  for (int ctr = 0; ctr < observations.size(); ++ctr) {
    // reference to the original observation
    LandmarkObs existing = observations[ctr];

    // transform into a new one
    LandmarkObs newLandmark;
    newLandmark.x = (existing.x * cos(particle.theta))
        - (existing.y * sin(particle.theta)) + particle.x;
    newLandmark.y = (existing.x * sin(particle.theta))
        + (existing.y * cos(particle.theta)) + particle.y;
    newLandmark.id = existing.id;

    // and store the transformed one
    converted.push_back(newLandmark);
  }
}

const double multivariateGaussDensity(const LandmarkObs& obs,
    const LandmarkObs &pred, double std_landmark[]) {
  double normalizer = 2. * M_PI * std_landmark[0] * std_landmark[1];

  double diff_x = pow(obs.x - pred.x, 2);
  double diff_y = pow(obs.y - pred.y, 2);
  double x_normalizer = 2 * std_landmark[0] * std_landmark[0];
  double y_normalizer = 2 * std_landmark[1] * std_landmark[1];

  double total = diff_x / x_normalizer + diff_y / y_normalizer;

  return exp(-total) / normalizer;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
    std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html

  // for all particles...
  for (int particleCtr = 0; particleCtr < particles.size(); particleCtr++) {

    Particle p = particles[particleCtr];

    // perform the 'dataAssociation'. To perform this operation we need two things:
    // 1. get the landmarks in the sensor range
    std::vector<LandmarkObs> filteredLandmarks;
    filterOnSensorRange(p, map_landmarks, sensor_range, filteredLandmarks);

    // 2. convert the observations into the maps coordinate system
    std::vector<LandmarkObs> transformedObservations;
    convertToMapCoordinates(observations, p, transformedObservations);

    // perform the data association
    dataAssociation(filteredLandmarks, transformedObservations);

    // calculate the weight
    double newWeight = 1.;
    for (int ctr = 0; ctr < transformedObservations.size(); ++ctr) {
      // find the corresponding filteredLandmark
      for (int filteredLandmarksCtr = 0;
          filteredLandmarksCtr < filteredLandmarks.size();
          ++filteredLandmarksCtr) {
        if (filteredLandmarks[filteredLandmarksCtr].id
            == transformedObservations[ctr].id) {
          // TODO: I'm sure indexes are not correct.
          newWeight *= multivariateGaussDensity(transformedObservations[ctr],
              filteredLandmarks[filteredLandmarksCtr], std_landmark);
          break;
        }
      }
    }

    // and store the weight
    particles[particleCtr].weight = newWeight;
    weights[particleCtr] = newWeight;

  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  discrete_distribution<int> idx_distr(weights.begin(), weights.end());

  std: vector<Particle> updated_weights;

  // process all weights
  for (int particleCtr = 0; particleCtr < num_particles; ++particleCtr) {
    // generate a new random number
    updated_weights.push_back(particles[idx_distr(gen)]);
  }

  // and store the new list of particles
  particles = updated_weights;

}

Particle ParticleFilter::SetAssociations(Particle particle,
    std::vector<int> associations, std::vector<double> sense_x,
    std::vector<double> sense_y) {
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  //Clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();

  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

  return particle;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best) {
  vector<double> v = best.sense_x;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best) {
  vector<double> v = best.sense_y;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}
