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

	int id;     ///< Identification number of the particle
	float x;    ///< x position in map/global coordinates [meters]
	float y;    ///< y position in map/global coordinates [meters]
	float theta; ///< Angle of the car relative to the horizontal [rad]
	float weight; ///< Probability that this particle is at the car location
	std::vector<int> associations;
	std::vector<float> sense_x;
	std::vector<float> sense_y;

    /** Default constructor
     *
     */
    Particle(int id_ = -1,
             float x_ = 0.0F,
             float y_ = 0.0F,
             float theta_ = 0.0F,
             float weight_ = 1.0F)
    : id(id_)
    , x(x_)
    , y(y_)
    , theta(theta_)
    , weight(weight_)
    {}
};


class ParticleFilter {
	
public:

    struct Config {

        float x;     ///< Initial x position [m] (simulated estimate from GPS)
        float y;     ///< Initial y position [m]
        float theta; ///< Initial orientation [rad]

        float std_dev_x;      ///< standard deviation of x [m]
        float std_dev_y;      ///< standard deviation of y [m]
        float std_dev_theta;  ///< standard deviation of theta [rad]

        /** Default constructor 
         *
         */
        Config(float x = 0.0,
               float y = 0.0,
               float theta = 0.0,
               float std_dev_x = 0.0,
               float std_dev_y = 0.0,
               float std_dev_theta = 0.0)
        : x(x)
        , y(y)
        , theta(theta)
        , std_dev_x(std_dev_x)
        , std_dev_y(std_dev_y)
        , std_dev_theta(std_dev_theta)
        {}
    };


	/** Constructor
     *
	 * @param num_particles Number of particles
     *
     */
    ParticleFilter(int nbParticles = 40);

	/** Destructor
     *
     */
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
	void prediction(float delta_t, float std_pos[], float velocity, float yaw_rate);
	
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
	void updateWeights(float sensor_range, float std_landmark[], const std::vector<LandmarkObs> &observations,
			const Map &map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/**
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<float> sense_x, std::vector<float> sense_y);
	
	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * Returns whether particle filter is initialized yet or not.
	 */
	inline const bool initialized() const {
		return mIsInitialized;
	}


    /**
     Returns the set of particles

     @return mParticules Vector of Particle that are tracking the car
     */
    inline const std::vector<Particle>& particles(){
        return mParticles;
    }

private:

    const int mNbParticles; ///< Number of particles to draw

    bool mIsInitialized;    ///< Flag, if filter is initialized

    std::vector<float> mWeights; ///< Vector of weights of all particles

    std::vector<Particle> mParticles; ///< Set of current particles

    Config mSavedConfig;  ///< Configuration for initialization
};



#endif /* PARTICLE_FILTER_H_ */
