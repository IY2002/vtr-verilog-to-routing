#include "weighted_centroid_move_generator_probabilistic.h"

WeightedCentroidProbMoveGenerator::WeightedCentroidProbMoveGenerator(PlacerState& placer_state,
                                                             e_reward_function reward_function,
                                                             vtr::RngContainer& rng)
    : CentroidProbMoveGenerator(placer_state, reward_function, rng) {
    weighted_ = true;
}
