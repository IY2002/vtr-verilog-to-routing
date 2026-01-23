#ifndef VPR_WEIGHTED_CENTROID_MOVE_GEN_PROB_H
#define VPR_WEIGHTED_CENTROID_MOVE_GEN_PROB_H

#include "centroid_move_generator_probabilistic.h"

/**
 * @brief Weighted Centroid move generator
 *
 * This move generator is inspired by analytical placers: model net connections as springs and 
 * calculate the force equilibrium location.
 *
 * @details This class inherits from CentroidMoveGenerator to avoid code duplication.
 *
 * For more details, please refer to:
 * "Learn to Place: FPGA Placement using Reinforcement Learning and Directed Moves", ICFPT2020
 */
class WeightedCentroidProbMoveGenerator : public CentroidProbMoveGenerator {
  public:
    WeightedCentroidProbMoveGenerator() = delete;
    WeightedCentroidProbMoveGenerator(PlacerState& placer_state,
                                  e_reward_function reward_function,
                                  vtr::RngContainer& rng);
};

#endif
