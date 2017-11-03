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

ParticleFilter::ParticleFilter()
    : mNbParticles(100),
    mIsInitialized(false)
    {}


bool ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    // This line creates 3 normal (Gaussian) distributions for x, y and theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    // Instantiation of a random nomber generator
    default_random_engine gen;

    for (int nbParticles = 0; nbParticles < mNbParticles; ++nbParticles) {
        double sample_x, sample_y, sample_theta;
        // I draw random samples from each gaussian distribution to
        // create the particles.
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
        sample_theta = dist_theta(gen);
        // I add the new particle to my vector of particles
        Particle newParticle(nbParticles, sample_x, sample_y, sample_theta);

        newParticle.print();

        mParticles.push_back(newParticle);
        mWeights.push_back(1.0F);
    }
    mIsInitialized = true;
    return true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    // Instantiation of a random nomber generator
    default_random_engine gen;
    assert(mNbParticles == mParticles.size());

    for (int nbParticles = 0; nbParticles < mNbParticles; ++nbParticles) {
        double newXpos, newYpos, newTheta;
        if (yaw_rate > ALMOST_ZERO){
            newTheta = mParticles[nbParticles].theta + delta_t * yaw_rate;

            newXpos = mParticles[nbParticles].x + ((velocity / yaw_rate) *
            (std::sin(newTheta) - std::sin(mParticles[nbParticles].theta)));

            newYpos = mParticles[nbParticles].y + ((velocity / yaw_rate) *
            (std::sin(mParticles[nbParticles].theta) - std::cos(newTheta)));
        }
        else
        {
            newXpos = mParticles[nbParticles].x + (delta_t * velocity * cos(mParticles[nbParticles].theta));
            newYpos = mParticles[nbParticles].y + (delta_t * velocity * sin(mParticles[nbParticles].theta));
            newTheta = mParticles[nbParticles].theta;
        }
        // This line creates 3 normal (Gaussian) distributions for x, y and theta
        // to sample from.
        normal_distribution<double> dist_x(newXpos, std_pos[0]);
        normal_distribution<double> dist_y(newYpos, std_pos[1]);
        normal_distribution<double> dist_theta(newTheta, std_pos[2]);

        mParticles[nbParticles].x = dist_x(gen);
        mParticles[nbParticles].y = dist_y(gen);
        mParticles[nbParticles].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    // The idea is to assign the id of the closest predicted landmark to the id member of the observations
    // landmarks.

    double distance;
    double minDistance;
    int index;

    for (int obs = 0 ; obs < observations.size() ; obs++){
        minDistance = std::numeric_limits<double>::max();
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
        // TODO: Check that this is possible
        // Improves efficiency
        //predicted.erase(predicted.begin() + index);
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

    // I first setup the multivariate gaussian constant parameters:
    // The normalizing factor: the standard deviations for x and y
    // are given as std_landmark[0] and std_landmark[1]
    double gausianNormalizingFactor = 1.0F / ( 2.0F * M_PI * std_landmark[0] * std_landmark[1]);

    double twoSigXSquared = 2.0F * std_landmark[0] * std_landmark[0];
    double twoSigYSquared = 2.0F * std_landmark[1] * std_landmark[1];

    // Sum of all assigned weight for the normalization step
    double cumulativeWeightsSum = 0.0F;

    std::cout << "New particles weights ... \n";

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
                double_t expTerm = (std::pow(transformedObservations[obs].x - reachableLandmarks[mapObsId].x, 2) / twoSigXSquared) +
                                 (std::pow(transformedObservations[obs].y - reachableLandmarks[mapObsId].y, 2) / twoSigYSquared);
                currentWeight *= gausianNormalizingFactor * std::exp(-expTerm);
            }
            else
            {
                std::cout << "ERROR ";
                currentWeight = 0.0F;
            }
        }

        std::cout << currentWeight << "\n";
        
        mWeights[particle] = currentWeight;
        mParticles[particle].weight = currentWeight;
        //cumulativeWeightsSum += currentWeight;
    }

    std::cout << "Sum of all weights : " << cumulativeWeightsSum << "\n";
/*
    // Normalization of the weights so that their sum is 1.0F
    for (int weight = 0; weight < mNbParticles ; weight++){
        mWeights[weight] /= cumulativeWeightsSum;
        mParticles[weight].weight = mWeights[weight];
        std::cout << mParticles[weight].weight << "\n";
    }*/
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
    discrete_distribution<double> index(mWeights.begin(), mWeights.end());

    // I assume the weights add up to 1
    for (int p = 0 ; p < mNbParticles ; p++){
        Particle particle = mParticles[index(rng)];
        newSetParticles.push_back(Particle(p,
                                           particle.x,
                                           particle.y,
                                           particle.theta,
                                           1.0F));
    }
    // Replace the former set of particles by the newly created.
    mParticles = newSetParticles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
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
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
