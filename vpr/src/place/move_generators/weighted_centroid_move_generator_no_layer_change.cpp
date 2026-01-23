#include "weighted_centroid_move_generator_no_layer_change.h"

WeightedCentroidLayerStuckMoveGenerator::WeightedCentroidLayerStuckMoveGenerator(PlacerState& placer_state,
                                                             e_reward_function reward_function,
                                                             vtr::RngContainer& rng)
    : CentroidLayerStuckMoveGenerator(placer_state, reward_function, rng) {
    weighted_ = true;
}
