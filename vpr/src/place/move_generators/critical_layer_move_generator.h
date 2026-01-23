// critical_layer_move_generator.h
#ifndef VPR_CRITICAL_LAYER_MOVE_GEN_H
#define VPR_CRITICAL_LAYER_MOVE_GEN_H

#include "move_generator.h"
#include "timing_place.h"

/**
 * @brief Critical Path Layer Optimization Move Generator
 * 
 * This move generator identifies blocks on critical paths and attempts to move them
 * to layers that minimize inter-layer hops along critical connections. It prioritizes
 * placing critically connected blocks on the same layer or adjacent layers.
 */
class CriticalLayerMoveGenerator : public MoveGenerator {
  public:
    CriticalLayerMoveGenerator() = delete;
    CriticalLayerMoveGenerator(PlacerState& placer_state,
                               e_reward_function reward_function,
                               vtr::RngContainer& rng);

  private:
    e_create_move propose_move(t_pl_blocks_to_be_moved& blocks_affected,
                              t_propose_action& proposed_action,
                              float rlim,
                              const t_placer_opts& placer_opts,
                              const PlacerCriticalities* criticalities) override;
    
    /**
     * @brief Calculate the optimal layer for a critical block
     * 
     * Determines the best layer by minimizing total weighted vertical distance
     * to connected blocks, where weights are based on connection criticality
     */
    int find_optimal_layer_for_critical_block(ClusterBlockId blk_id,
                                            ClusterNetId critical_net,
                                            int critical_pin,
                                            const PlacerCriticalities* criticalities);
};

#endif