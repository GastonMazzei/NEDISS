//
// Created by m4zz31 on 30/10/21.
//

#include "reproductibility.h"

void reproductibility_lock(unsigned int SEED) {
    boost::mt19937 gener(1);
    boost::normal_distribution<> normal(0,1);
    boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal);
    rng.engine().seed(SEED);
    rng.distribution().reset();
}