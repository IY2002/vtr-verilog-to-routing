#ifndef VPR_LAYER_SWAP_MOVE_GENERATOR_H
#define VPR_LAYER_SWAP_MOVE_GENERATOR_H

#include "move_generator.h"

/**
 * @brief Layer swap move generator
 * 
 * This move generator proposes moves that swap a block to a different layer
 * at the same X,Y location. This is useful for 3D architectures where blocks
 * can be placed on multiple layers.
 */
class LayerSwapMoveGenerator : public MoveGenerator {
  public:
    LayerSwapMoveGenerator() = delete;
    explicit LayerSwapMoveGenerator(PlacerState& placer_state,
                                   e_reward_function reward_function,
                                   vtr::RngContainer& rng);

  private:
    e_create_move propose_move(t_pl_blocks_to_be_moved& blocks_affected,
                              t_propose_action& proposed_action,
                              float rlim,
                              const t_placer_opts& placer_opts,
                              const PlacerCriticalities* criticalities) override;
};

#endif /* VPR_LAYER_SWAP_MOVE_GENERATOR_H */