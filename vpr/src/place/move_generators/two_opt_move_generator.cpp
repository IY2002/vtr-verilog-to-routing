#include "two_opt_move_generator.h"

#include "globals.h"
#include "place_constraints.h"
#include "placer_state.h"
#include "move_utils.h"

#include <algorithm>
#include <vector>

TwoOptMoveGenerator::TwoOptMoveGenerator(PlacerState& placer_state,
                                         e_reward_function reward_function,
                                         vtr::RngContainer& rng)
    : MoveGenerator(placer_state, reward_function, rng) {}

e_create_move TwoOptMoveGenerator::propose_move(t_pl_blocks_to_be_moved& blocks_affected,
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

    VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, "2-opt Move Choose Block %d - rlim %f\n", size_t(b_from), rlim);

    if (!b_from) { // No movable block found
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, "\tNo movable block found\n");
        return e_create_move::ABORT;
    }

    // Get the type and location of the from block
    t_pl_loc from = block_locs[b_from].loc;
    auto cluster_from_type = cluster_ctx.clb_nlist.block_type(b_from);
    auto grid_from_type = device_ctx.grid.get_physical_type({from.x, from.y, from.layer});
    VTR_ASSERT(is_tile_compatible(grid_from_type, cluster_from_type));

    // Collect all compatible connected blocks
    std::vector<ClusterBlockId> compatible_connected_blocks;
    
    // Iterate over all pins of the from block
    for (ClusterPinId pin_id : cluster_ctx.clb_nlist.block_pins(b_from)) {
        ClusterNetId net_id = cluster_ctx.clb_nlist.pin_net(pin_id);
        
        if (cluster_ctx.clb_nlist.net_is_ignored(net_id))
            continue;
            
        // Skip high fanout nets
        if (int(cluster_ctx.clb_nlist.net_pins(net_id).size()) > placer_opts.place_high_fanout_net)
            continue;
        
        // Check driver block
        ClusterBlockId driver_block = cluster_ctx.clb_nlist.net_driver_block(net_id);
        if (driver_block && driver_block != b_from) {
            auto connected_type = cluster_ctx.clb_nlist.block_type(driver_block);
            
            // Check if the blocks have compatible types (can be swapped)
            if (is_tile_compatible(grid_from_type, connected_type)) {
                t_pl_loc connected_loc = block_locs[driver_block].loc;
                auto grid_connected_type = device_ctx.grid.get_physical_type({connected_loc.x, connected_loc.y, connected_loc.layer});
                
                if (is_tile_compatible(grid_connected_type, cluster_from_type)) {
                    // Check if we haven't already added this block
                    if (std::find(compatible_connected_blocks.begin(), compatible_connected_blocks.end(), driver_block) == compatible_connected_blocks.end()) {
                        compatible_connected_blocks.push_back(driver_block);
                    }
                }
            }
        }
        
        // Check all sink blocks
        for (ClusterPinId sink_pin : cluster_ctx.clb_nlist.net_sinks(net_id)) {
            ClusterBlockId sink_block = cluster_ctx.clb_nlist.pin_block(sink_pin);
            if (sink_block != b_from) {
                auto connected_type = cluster_ctx.clb_nlist.block_type(sink_block);
                
                // Check if the blocks have compatible types (can be swapped)
                if (is_tile_compatible(grid_from_type, connected_type)) {
                    t_pl_loc connected_loc = block_locs[sink_block].loc;
                    auto grid_connected_type = device_ctx.grid.get_physical_type({connected_loc.x, connected_loc.y, connected_loc.layer});
                    
                    if (is_tile_compatible(grid_connected_type, cluster_from_type)) {
                        // Check if we haven't already added this block
                        if (std::find(compatible_connected_blocks.begin(), compatible_connected_blocks.end(), sink_block) == compatible_connected_blocks.end()) {
                            compatible_connected_blocks.push_back(sink_block);
                        }
                    }
                }
            }
        }
    }
    
    // If no compatible connected blocks found, abort
    if (compatible_connected_blocks.empty()) {
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, "\tNo compatible connected blocks found\n");
        return e_create_move::ABORT;
    }
    
    // Randomly select one of the connected blocks
    int random_index = rng_.irand(compatible_connected_blocks.size() - 1);
    ClusterBlockId b_to = compatible_connected_blocks[random_index];
    
    // Get the location of the selected connected block
    t_pl_loc to = block_locs[b_to].loc;
    
    VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, "\tSwapping block %d at (%d,%d,%d) with connected block %d at (%d,%d,%d)\n", 
                   size_t(b_from), from.x, from.y, from.layer,
                   size_t(b_to), to.x, to.y, to.layer);
    
    // Create the move (swap b_from with b_to's location)
    e_create_move create_move = ::create_move(blocks_affected, b_from, to, blk_loc_registry);
    
    // Check that all the blocks affected by the move would still be in a legal floorplan region after the swap
    if (!floorplan_legal(blocks_affected)) {
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, "\tMove aborted - floorplan illegal\n");
        return e_create_move::ABORT;
    }
    
    return create_move;
}