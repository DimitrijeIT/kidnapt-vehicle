/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

using std::cout;
using std::endl;

static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  is_initialized= true;
  const double default_weight = 1;
  double sample_x, sample_y, sample_theta;
  std::normal_distribution<double> dist_x (x, std[0]);
  std::normal_distribution<double> dist_y (y, std[1]);
  std::normal_distribution<double> dist_theta( theta, std[2]);

  // cout << "Init" << endl;
  for(int i = 0; i< num_particles; ++i)
  {
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);
       
    Particle p;
    p.id = i;
    p.x = sample_x;
    p.y = sample_y;
    p.theta = sample_theta;
    p.weight = default_weight;

    particles.push_back(p);
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // cout << "------ Predict state " << endl;
  
    std::normal_distribution<double> noise_x (0, std_pos[0]);
    std::normal_distribution<double> noise_y (0, std_pos[1]);
    std::normal_distribution<double> noise_theta (0, std_pos[2]);

  for(Particle& p: particles)
  {

    if(fabs(yaw_rate) < 0.00001)
    {
      p.x = p.x + velocity * delta_t * cos(p.theta);
      p.y = p.y + velocity * delta_t * sin(p.theta);
    }
    else
    {
      p.x = p.x + (velocity/yaw_rate)*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
      p.y = p.y + (velocity/yaw_rate)*(cos(p.theta) - cos(p.theta + yaw_rate*delta_t));
      p.theta = p.theta + yaw_rate * delta_t;
    }

    p.x += noise_x(gen);
    p.y += noise_y(gen);
    p.theta += noise_theta(gen);
  }
  // cout << "------ Predict end " << endl;
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // cout << "Update weight state " << endl;
  for(Particle& p: particles)
  {
    //Find landmarks in sensor range from particles
   vector<Map::single_landmark_s> landmarks_in_sensor_range; 
   for(const Map::single_landmark_s& landmark: map_landmarks.landmark_list)
   {
     double distance_particle_landmark = dist(p.x, p.y, landmark.x_f, landmark.y_f);
     if(distance_particle_landmark < sensor_range){
       landmarks_in_sensor_range.push_back(landmark);
     }
   }

   if(landmarks_in_sensor_range.size() < observations.size()) 
   {
     //Partivle is way off and has less landmarks in rannge than observed
     // cout << " -------- Particle id = " << p.id << " has less landmarks on range than observed" << endl;  
     p.weight = 0;
     continue;
   }

   p.weight = 1;
   for(const LandmarkObs& obs: observations)
   {
      // Homogenous Transformation transform to map coordinates
      double x_map = p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);
      double y_map = p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);

      // Find association to landmark neerest to observation
      double min_distance = std::numeric_limits<double>::max();
      double distance = 0;
      Map::single_landmark_s associated_landmark;

      for(const Map::single_landmark_s& landmark: landmarks_in_sensor_range)
      {
        distance = dist(x_map,y_map, landmark.x_f, landmark.y_f);
        if(distance < min_distance)
        {
          min_distance = distance;
          associated_landmark = landmark;
        }
      }

      double observation_weight = multiv_prob(std_landmark[0], std_landmark[1], 
                                          x_map, y_map, 
                                          associated_landmark.x_f, associated_landmark.y_f);

      // cout << "   Observed weiht = " << observation_weight;                              
        p.weight *= observation_weight;
        // cout << "  particle weight = " << p.weight;
      // cout << endl;
   }
   // cout << "Weeight id =  " << p.id << " weight = " << p.weight << " Iteration "<< iteration_num << endl;
  }
  iteration_num++;
  // cout << "Update weight end " << endl;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  //Update weights vector
   // cout << "Resample state " << endl;
  std::vector<double> weights_local;
  double sum = 0;
  for(const Particle& p : particles)
  {
    sum += p.weight;
  }

  double normalized_sum = 0;
  for(const Particle& p: particles)
  {
    double normalized_weight = p.weight/sum;
    weights_local.push_back(normalized_weight);
    normalized_sum += normalized_weight;
  }

  // cout << " Normalized sum shouhdl be 1 = " << normalized_sum << endl;

  std::discrete_distribution<> distrbution(weights_local.begin(), weights_local.end());
  std::vector<Particle> new_particles;
  for(int i = 0; i < num_particles; ++i)
  {
    int index = distrbution(gen);
    Particle p = particles.at(index);
    // cout << "For loop distribution index = " << index << " particle id = " << p.id << " weight = " << p.weight << endl; 
    new_particles.push_back(particles.at(index));
  }

  // cout << "For end loop distribution gen  " << endl;
  particles = new_particles;
  // cout << "Resample end " << endl;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);


  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}