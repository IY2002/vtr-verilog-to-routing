#include "layer_swap_ranged_move_generator.h"

#include "globals.h"
#include "place_constraints.h"
#include "placer_state.h"
#include "move_utils.h"

LayerSwapRangedMoveGenerator::LayerSwapRangedMoveGenerator(PlacerState& placer_state,
                                               e_reward_function reward_function,
                                               vtr::RngContainer& rng)
    : MoveGenerator(placer_state, reward_function, rng) {}

e_create_move LayerSwapRangedMoveGenerator::propose_move(t_pl_blocks_to_be_moved& blocks_affected,
                                                   t_propose_action& proposed_action,
                                                   float rlim,
                                                   const t_placer_opts& placer_opts,
                                                   const PlacerCriticalities* /*criticalities*/) {
    const auto& cluster_ctx = g_vpr_ctx.clustering();
    const auto& device_ctx = g_vpr_ctx.device();
    auto& placer_state = placer_state_.get();
    const auto& block_locs = placer_state.block_locs();
    const auto& blk_loc_registry = placer_state.blk_loc_registry();

    // Find a movable block based on blk_type
    ClusterBlockId b_from = propose_block_to_move(placer_opts,
                                                  proposed_action.logical_blk_type_index,
                                                  /*highly_crit_block=*/false,
                                                  /*net_from=*/nullptr,
                                                  /*pin_from=*/nullptr,
                                                  placer_state,
                                                  rng_);

    VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                   "Layer Swap Move Choose Block %d - rlim %f\n", 
                   size_t(b_from), 
                   rlim);

    if (!b_from) { // No movable block found
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                       "\tNo movable block found\n");
        return e_create_move::ABORT;
    }

    // Get current location of the block
    t_pl_loc from = block_locs[b_from].loc;
    auto cluster_from_type = cluster_ctx.clb_nlist.block_type(b_from);
    auto grid_from_type = device_ctx.grid.get_physical_type({from.x, from.y, from.layer});
    VTR_ASSERT(is_tile_compatible(grid_from_type, cluster_from_type));

    // Check if we have multiple layers
    int num_layers = device_ctx.grid.get_num_layers();
    if (num_layers <= 1) {
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                       "\tLayer swap aborted - only one layer available\n");
        return e_create_move::ABORT;
    }

    // Build list of candidate layers (all layers except current)
    std::vector<int> candidate_layers;
    for (int layer = 0; layer < num_layers; layer++) {
        if (layer != from.layer) {
            candidate_layers.push_back(layer);
        }
    }

    // Try random layers until we find a valid move
    bool found_valid_move = false;
    t_pl_loc to;
    
    while (!candidate_layers.empty() && !found_valid_move) {
        // Pick a random layer from candidates
        int random_index = rng_.irand(candidate_layers.size() - 1);
        int target_layer = candidate_layers[random_index];
        
        // Remove this layer from candidates
        candidate_layers.erase(candidate_layers.begin() + random_index);
        
        // Create a temporary "from" location on the target layer
        // This is used as the center point for finding a location within rlim
        t_pl_loc from_on_target_layer = from;
        from_on_target_layer.layer = target_layer;
        float search_range = std::max(3.0f, rlim); // Ensure a minimum search range of 3
        // Try to find a valid location within rlim on the target layer
        if (find_to_loc_uniform(cluster_from_type, search_range, from_on_target_layer, to, b_from, blk_loc_registry, rng_, true)) {
            // Make sure we actually changed layers (find_to_loc_uniform might return same layer)
            if (to.layer == target_layer) {
                found_valid_move = true;
            }
        }
    }

    if (!found_valid_move) {
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                       "\tNo valid location found for layer swap\n");
        return e_create_move::ABORT;
    }

    VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                   "\tSwapping block from layer %d at (%d,%d) to layer %d at (%d,%d)\n", 
                   from.layer, from.x, from.y, to.layer, to.x, to.y);

    // Create the move
    e_create_move create_move = ::create_move(blocks_affected, b_from, to, blk_loc_registry);

    // Check that all the blocks affected by the move would still be in a legal floorplan region
    if (!floorplan_legal(blocks_affected)) {
        return e_create_move::ABORT;
    }

    return create_move;
}