#include <iostream>
#include "gtest/gtest.h"
#include <vector>

#include "particle_filter.h"

#define MAX_ABSOLUTE_ERROR 1e-4

class ParticleFilterTest : public ::testing::Test
{
public:

    ParticleFilterTest()
        : mNbParticles(5)
        , mParticleFilter(mNbParticles)
        {}

    virtual void SetUp() {

    }

    virtual void TearDown() {
    
    }

    int             mNbParticles;

    ParticleFilter  mParticleFilter;

};


TEST_F(ParticleFilterTest, DataAssociationTest) {

    // Case 1
    std::vector<LandmarkObs> predictedLandmarks;
    std::vector<LandmarkObs> observedLandmarks;

    predictedLandmarks.push_back(LandmarkObs(0, 0.0F,  0.0F));
    predictedLandmarks.push_back(LandmarkObs(1, 1.0F,  1.0F));
    predictedLandmarks.push_back(LandmarkObs(2, 2.0F, -2.0F));

    observedLandmarks.push_back(LandmarkObs(0,  1.1F,  1.0F));
    observedLandmarks.push_back(LandmarkObs(1,  2.0F, -2.5F));
    observedLandmarks.push_back(LandmarkObs(2, -0.1F,  0.1F));

    mParticleFilter.dataAssociation(predictedLandmarks, observedLandmarks);

    ASSERT_EQ(observedLandmarks[0].id, 1);
    ASSERT_EQ(observedLandmarks[1].id, 2);
    ASSERT_EQ(observedLandmarks[2].id, 0);

    // Case 2
    observedLandmarks.push_back(LandmarkObs(3, 3.0F,  -2.0F));
    mParticleFilter.dataAssociation(predictedLandmarks, observedLandmarks);

    ASSERT_EQ(observedLandmarks[0].id, 1);
    ASSERT_EQ(observedLandmarks[1].id, 2);
    ASSERT_EQ(observedLandmarks[2].id, 0);
    ASSERT_EQ(observedLandmarks[3].id, 2);

    // Case 3 - Wrong ids - Shouldn't matter
    predictedLandmarks.push_back(LandmarkObs(0, 1.09F,  0.99F)); // 4th - indice 3
    mParticleFilter.dataAssociation(predictedLandmarks, observedLandmarks);

    ASSERT_EQ(observedLandmarks[0].id, 3);
    ASSERT_EQ(observedLandmarks[1].id, 2);
    ASSERT_EQ(observedLandmarks[2].id, 0);
    ASSERT_EQ(observedLandmarks[3].id, 2);
}


TEST_F(ParticleFilterTest, InitTest) {

    // Case 1 : Standard deviations set to 0
    //float std_devs[3] = {0.0F, 0.0F, 0.0F};
    ParticleFilter::Config config(5.0F, 5.0F, 10.0F, 0.0F, 0.0F, 0.0F);
    mParticleFilter.init(config);
    std::vector<Particle> particles = mParticleFilter.particles();

    ASSERT_EQ(particles.size(), mNbParticles);

    for (int particle = 0 ; particle < mNbParticles ; particle++){
        ASSERT_EQ(particles[particle].x, 5.0F);
        ASSERT_EQ(particles[particle].y, 5.0F);
        ASSERT_EQ(particles[particle].theta, 10.0F);
    }

    // Case 2 : ?
}
