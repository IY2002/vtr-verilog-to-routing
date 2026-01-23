#ifndef VPR_WEIGHTED_CENTROID_MOVE_GEN_LAYER_STUCK_H
#define VPR_WEIGHTED_CENTROID_MOVE_GEN_LAYER_STUCK_H

#include "centroid_move_generator_no_layer_change.h"

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
class WeightedCentroidLayerStuckMoveGenerator : public CentroidLayerStuckMoveGenerator {
  public:
    WeightedCentroidLayerStuckMoveGenerator() = delete;
    WeightedCentroidLayerStuckMoveGenerator(PlacerState& placer_state,
                                  e_reward_function reward_function,
                                  vtr::RngContainer& rng);
};

#endif
