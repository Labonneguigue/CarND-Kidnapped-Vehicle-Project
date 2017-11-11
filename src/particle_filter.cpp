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
#include <cassert>
#include "particle_filter.h"
#include <chrono>
#include <functional>

#define ALMOST_ZERO 0.00001
using namespace std;

ParticleFilter::ParticleFilter(int nbParticles)
    : mNbParticles(nbParticles),
    mIsInitialized(false)
    {}


bool ParticleFilter::init(Config& config) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    // This line creates 3 normal (Gaussian) distributions for x, y and theta
    normal_distribution<float> dist_x(config.x, config.std_dev_x);
    normal_distribution<float> dist_y(config.y, config.std_dev_y);
    normal_distribution<float> dist_theta(config.theta, config.std_dev_theta);

    // Instantiation of a random nomber generator
    default_random_engine gen;
    mParticles.clear();

    for (int nbParticles = 0; nbParticles < mNbParticles; ++nbParticles) {
        float sample_x, sample_y, sample_theta;
        // I draw random samples from each gaussian distribution to
        // create the particles.
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
        sample_theta = dist_theta(gen);
        // I add the new particle to my vector of particles
        Particle newParticle(nbParticles, sample_x, sample_y, sample_theta);

        //newParticle.print();

        mParticles.push_back(newParticle);
        mWeights.push_back(1.0F);
    }
    mIsInitialized = true;
    mSavedConfig = config;
    return true;
}

void ParticleFilter::prediction(float delta_t, float std_pos[], float velocity, float yaw_rate) {

    // Instantiation of a random nomber generator
    default_random_engine gen;
    assert(mNbParticles == mParticles.size());

    // The yaw_rate angle sent from the simulator can be non-normalized
    NormalizeAngle(yaw_rate);

    for (int nbParticles = 0; nbParticles < mNbParticles; ++nbParticles) {
        float newXpos, newYpos, newTheta;

        // There are 2 different sets of equations that determine the evolution
        // of x, y, and theta depending on whether the yaw rate is almost zero
        // or not.
        if (yaw_rate > ALMOST_ZERO){

            newTheta = mParticles[nbParticles].theta + (delta_t * yaw_rate);
            NormalizeAngle(newTheta);
            newXpos = mParticles[nbParticles].x + ((velocity / yaw_rate) *
            (std::sin(newTheta) - std::sin(mParticles[nbParticles].theta)));

            newYpos = mParticles[nbParticles].y + ((velocity / yaw_rate) *
            (std::cos(mParticles[nbParticles].theta) - std::cos(newTheta)));

        }
        else
        {
            newXpos = mParticles[nbParticles].x + (delta_t * velocity * cos(mParticles[nbParticles].theta));
            newYpos = mParticles[nbParticles].y + (delta_t * velocity * sin(mParticles[nbParticles].theta));
            newTheta = mParticles[nbParticles].theta;
        }

        // This line creates 3 normal (Gaussian) distributions for x, y and theta
        // to sample from. This introduces random noise that helps the particles
        // spread over the map and fill every position locations for the car.
        normal_distribution<float> dist_x(newXpos, std_pos[0]);
        normal_distribution<float> dist_y(newYpos, std_pos[1]);
        normal_distribution<float> dist_theta(newTheta, std_pos[2]*2);

        mParticles[nbParticles].x = dist_x(gen);
        mParticles[nbParticles].y = dist_y(gen);
        mParticles[nbParticles].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // The idea is to assign the id of the closest predicted landmark to the id member of the observations
    // landmarks.

    float distance;
    float minDistance;
    int index;

    for (int obs = 0 ; obs < observations.size() ; obs++){
        minDistance = std::numeric_limits<float>::max();
        index = 0;
        for (int pre = 0 ; pre < predicted.size() ; pre++){
            distance = dist(observations[obs].x, observations[obs].y, predicted[pre].x, predicted[pre].y);
            if ( distance < minDistance ){
                minDistance = distance;
                index = pre;
            }
        }
        // index of the closest prediction is in the index variable.
        if (index < predicted.size()){
            observations[obs].id = index;
        }
        else
        {
            // error
            observations[obs].id = -1;
        }
    }
}

void ParticleFilter::updateWeights(float sensor_range, float std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // I first setup the multivariate gaussian constant parameters into a structure
    // The normalizing factor: the standard deviations for x and y
    // are given as std_landmark[0] and std_landmark[1]

    MultivariateGaussian multivariateGaussian(std_landmark[0], std_landmark[1]);

    for ( int particle = 0 ; particle < mNbParticles ; particle++){
        // For each particle, I first gather in a vector every landmarks that
        // are supposed to be detected (ones that lie within sensor_range)
        vector<LandmarkObs> reachableLandmarks;
        for (int landmark = 0 ; landmark < map_landmarks.landmark_list.size() ; landmark++){
            // Position of both the map landmarks and the particles are expressed
            // in the map coordinates.
            if (dist(map_landmarks.landmark_list[landmark].x_f,
                     map_landmarks.landmark_list[landmark].y_f,
                     mParticles[particle].x,
                     mParticles[particle].y) < sensor_range){
                reachableLandmarks.push_back(LandmarkObs(map_landmarks.landmark_list[landmark].id_i,
                                                         map_landmarks.landmark_list[landmark].x_f,
                                                         map_landmarks.landmark_list[landmark].y_f));
            }
        }

        assert(!reachableLandmarks.empty());

        // I now convert the car observations from the car coordinate system to the map one.
        vector<LandmarkObs> transformedObservations;
        for (int obs = 0 ; obs < observations.size() ; obs++){
            Eigen::MatrixXd T = getHomogenousTransformationMatrix(mParticles[particle].theta,
                                                                  mParticles[particle].x,
                                                                  mParticles[particle].y);
            // Vector of observation x and y in car coordinate system.
            Eigen::VectorXd obsVec_carC(3);
            obsVec_carC << observations[obs].x, observations[obs].y, 1.0F;
            Eigen::VectorXd obsVec_mapC = T * obsVec_carC;
            transformedObservations.push_back(LandmarkObs(observations[obs].id,
                                                          obsVec_mapC[0],
                                                          obsVec_mapC[1]));
        }

        // Association of the observations with their most likely landmark.
        // Both are in car coordinates
        dataAssociation(reachableLandmarks, transformedObservations);

        // A multi-variate gaussian distribution would now provide the probability
        // that the observation was effectivelly performed on the compared landmark
        // A multiplication of all of these multi-variate gaussians gives a probabiliy
        // that the studied particle occupies the actual car location.
        // I start with a weight probablity of 1.
        double_t currentWeight = 1.0F;

        for (int obs = 0 ; obs < transformedObservations.size() ; obs++){
            int mapObsId = transformedObservations[obs].id;
            if ( mapObsId != -1){
                currentWeight *= multivariateGaussian.correlation(transformedObservations[obs].x,
                                                                 reachableLandmarks[mapObsId].x,
                                                                 transformedObservations[obs].y,
                                                                 reachableLandmarks[mapObsId].y);
            }
            else
            {
                std::cout << "ERROR ";
                currentWeight = 0.0F;
            }
        }

        mWeights[particle] = currentWeight;
        mParticles[particle].weight = currentWeight;
    }
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    vector<Particle> newSetParticles;

    // A way to generate random numbers using the Mersenne Twister algorithm
    // found here : https://www.guyrutenberg.com/2014/05/03/c-mt19937-example/
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    static std::mt19937 rng(seed);

    // Creation of a distribution from which the probability to draw one is proportional
    // to the weight inputted.
    discrete_distribution<int> index(mWeights.begin(), mWeights.end());

    // I assume the weights add up to 1
    for (int p = 0 ; p < mNbParticles ; p++){
        Particle particle = mParticles[index(rng)];
        newSetParticles.push_back(particle);
    }
    // Replace the former set of particles by the newly created.
    mParticles = std::move(newSetParticles);
}


Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<float> sense_x, std::vector<float> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    // TODO: find a better way
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<float> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<float> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
