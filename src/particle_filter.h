/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;

    Particle(int id_ = -1,
             double x_ = 0.0F,
             double y_ = 0.0F,
             double theta_ = 0.0F,
             double weight_ = 1.0F)
    : id(id_)
    , x(x_)
    , y(y_)
    , theta(theta_)
    , weight(weight_)
    {}

    void print(){
        std::cout << "Particle " << id << " " << x << " " << y << " " << theta << "\n";
    }
};


class ParticleFilter {
	
public:

    struct Config {
        
        double x;     ///< Initial x position [m] (simulated estimate from GPS)
        double y;     ///< Initial y position [m]
        double theta; ///< Initial orientation [rad]

        double std_dev_x;      ///< standard deviation of x [m]
        double std_dev_y;      ///< standard deviation of y [m]
        double std_dev_theta;  ///< standard deviation of theta [rad]

        /** Default constructor 
         *
         */
        Config(double x = 0.0,
               double y = 0.0,
               double theta = 0.0,
               double std_dev_x = 0.0,
               double std_dev_y = 0.0,
               double std_dev_theta = 0.0)
        : x(x)
        , y(y)
        , theta(theta)
        , std_dev_x(std_dev_x)
        , std_dev_y(std_dev_y)
        , std_dev_theta(std_dev_theta)
        {}
    };


	// Constructor
	// @param num_particles Number of particles
    ParticleFilter(int nbParticles = 100);

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
     *
     * @param[in] config  Configuration structure holding x, y, theta and their
                          respective standard deviations
     *
     * @return bool True if the initialization was successful
	 */
	bool init(Config& config);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations,
			const Map &map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

    // TODO : 4 next methods have nothing to do here
    // Should be moved into Particle structure

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y);
	
	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	inline const bool initialized() const {
		return mIsInitialized;
	}

    inline const std::vector<Particle>& particles(){
        return mParticles;
    }

private:

    // Number of particles to draw
    const int mNbParticles;

    // Flag, if filter is initialized
    // TODO: delete that
    bool mIsInitialized;

    // Vector of weights of all particles
    std::vector<double> mWeights;

    // Set of current particles
    std::vector<Particle> mParticles;

    Config mSavedConfig;
};



#endif /* PARTICLE_FILTER_H_ */
