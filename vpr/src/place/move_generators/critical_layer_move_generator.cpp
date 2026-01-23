// critical_layer_move_generator.cpp
#include "critical_layer_move_generator.h"
#include "globals.h"
#include "place_constraints.h"
#include "placer_state.h"
#include "move_utils.h"

#include <vector>
#include <limits>

CriticalLayerMoveGenerator::CriticalLayerMoveGenerator(PlacerState& placer_state,
                                                     e_reward_function reward_function,
                                                     vtr::RngContainer& rng)
    : MoveGenerator(placer_state, reward_function, rng) {}

e_create_move CriticalLayerMoveGenerator::propose_move(t_pl_blocks_to_be_moved& blocks_affected,
                                                       t_propose_action& proposed_action,
                                                       float rlim,
                                                       const t_placer_opts& placer_opts,
                                                       const PlacerCriticalities* criticalities) {
    const auto& cluster_ctx = g_vpr_ctx.clustering();
    const auto& device_ctx = g_vpr_ctx.device();
    auto& placer_state = placer_state_.get();
    const auto& block_locs = placer_state.block_locs();
    const auto& blk_loc_registry = placer_state.blk_loc_registry();

    // Find a critical block to move
    ClusterNetId net_from;
    int pin_from;
    ClusterBlockId b_from = propose_block_to_move(placer_opts,
                                                 proposed_action.logical_blk_type_index,
                                                 /*highly_crit_block=*/true,
                                                 &net_from,
                                                 &pin_from,
                                                 placer_state,
                                                 rng_);

    VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                   "Critical Layer Move Choose Block %d on critical net %d pin %d - rlim %f\n", 
                   size_t(b_from), size_t(net_from), pin_from, rlim);

    if (!b_from) { // No movable block found
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                       "\tNo movable critical block found\n");
        return e_create_move::ABORT;
    }

    // Get current location
    t_pl_loc from = block_locs[b_from].loc;
    auto cluster_from_type = cluster_ctx.clb_nlist.block_type(b_from);
    auto grid_from_type = device_ctx.grid.get_physical_type({from.x, from.y, from.layer});
    VTR_ASSERT(is_tile_compatible(grid_from_type, cluster_from_type));

    // Check if we have multiple layers
    int num_layers = device_ctx.grid.get_num_layers();
    if (num_layers <= 1) {
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                       "\tCritical layer move aborted - only one layer available\n");
        return e_create_move::ABORT;
    }

    // Find the optimal layer for this critical block
    int target_layer = find_optimal_layer_for_critical_block(b_from, net_from, pin_from, criticalities);
    
    if (target_layer == from.layer) {
        // Already on optimal layer, try second best
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                       "\tBlock already on optimal layer %d\n", target_layer);
        return e_create_move::ABORT;
    }

    // Create a temporary "from" location on the target layer
    t_pl_loc from_on_target_layer = from;
    from_on_target_layer.layer = target_layer;
    
    // Find a location within rlim on the target layer
    t_pl_loc to;
    if (!find_to_loc_uniform(cluster_from_type, rlim, from_on_target_layer, to, b_from, blk_loc_registry, rng_)) {
        VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                       "\tNo valid location found on target layer %d\n", target_layer);
        return e_create_move::ABORT;
    }

    // Ensure we're actually moving to the target layer
    if (to.layer != target_layer) {
        to.layer = target_layer;
        // Re-validate the location
        auto grid_to_type = device_ctx.grid.get_physical_type({to.x, to.y, to.layer});
        if (!is_tile_compatible(grid_to_type, cluster_from_type)) {
            return e_create_move::ABORT;
        }
    }

    VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                   "\tMoving critical block from layer %d to layer %d at (%d,%d)\n", 
                   from.layer, to.layer, to.x, to.y);

    // Create the move
    e_create_move create_move = ::create_move(blocks_affected, b_from, to, blk_loc_registry);

    // Check floorplan legality
    if (!floorplan_legal(blocks_affected)) {
        return e_create_move::ABORT;
    }

    return create_move;
}

int CriticalLayerMoveGenerator::find_optimal_layer_for_critical_block(
    ClusterBlockId blk_id,
    ClusterNetId critical_net,
    int critical_pin,
    const PlacerCriticalities* criticalities) {
    
    const auto& cluster_ctx = g_vpr_ctx.clustering();
    const auto& device_ctx = g_vpr_ctx.device();
    const auto& block_locs = placer_state_.get().block_locs();
    const auto& blk_loc_registry = placer_state_.get().blk_loc_registry();
    
    int num_layers = device_ctx.grid.get_num_layers();
    std::vector<float> layer_scores(num_layers, 0.0f);
    
    // For each layer, calculate the weighted sum of vertical distances to critical connections
    for (int target_layer = 0; target_layer < num_layers; target_layer++) {
        float score = 0.0f;
        float total_weight = 0.0f;
        
        // Check all pins of the block
        for (ClusterPinId pin_id : cluster_ctx.clb_nlist.block_pins(blk_id)) {
            ClusterNetId net_id = cluster_ctx.clb_nlist.pin_net(pin_id);
            
            if (cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
                continue;
            }
            
            // // Skip high fanout nets
            // if (cluster_ctx.clb_nlist.net_sinks(net_id).size() > HIGH_FANOUT_NET_LIM) {
            //     continue;
            // }
            
            float criticality_weight = 1.0f;
            
            if (cluster_ctx.clb_nlist.pin_type(pin_id) == PinType::DRIVER) {
                // For driver pins, consider all sinks
                for (auto sink_pin_id : cluster_ctx.clb_nlist.net_sinks(net_id)) {
                    int ipin = cluster_ctx.clb_nlist.pin_net_index(sink_pin_id);
                    
                    // Get criticality if available
                    if (criticalities != nullptr) {
                        criticality_weight = criticalities->criticality(net_id, ipin);
                    }
                    
                    // Extra weight for the critical connection that triggered this move
                    if (net_id == critical_net && ipin == critical_pin) {
                        criticality_weight *= 2.0f;
                    }
                    
                    // Get sink location
                    t_physical_tile_loc sink_loc = blk_loc_registry.get_coordinate_of_pin(sink_pin_id);
                    
                    // Calculate vertical distance
                    int layer_distance = std::abs(target_layer - sink_loc.layer_num);
                    
                    // Add to score (lower is better)
                    score += layer_distance * criticality_weight;
                    total_weight += criticality_weight;
                }
            } else {
                // For sink pins, consider the driver
                int ipin = cluster_ctx.clb_nlist.pin_net_index(pin_id);
                
                // Get criticality if available
                if (criticalities != nullptr) {
                    criticality_weight = criticalities->criticality(net_id, ipin);
                }
                
                // Extra weight for the critical connection
                if (net_id == critical_net && ipin == critical_pin) {
                    criticality_weight *= 2.0f;
                }
                
                ClusterPinId driver_pin = cluster_ctx.clb_nlist.net_driver(net_id);
                t_physical_tile_loc driver_loc = blk_loc_registry.get_coordinate_of_pin(driver_pin);
                
                // Calculate vertical distance
                int layer_distance = std::abs(target_layer - driver_loc.layer_num);
                
                // Add to score
                score += layer_distance * criticality_weight;
                total_weight += criticality_weight;
            }
        }
        
        // Normalize score
        if (total_weight > 0) {
            layer_scores[target_layer] = score / total_weight;
        } else {
            layer_scores[target_layer] = std::numeric_limits<float>::max();
        }
    }
    
    // Find the layer with minimum score (least total critical vertical distance)
    int best_layer = 0;
    float min_score = layer_scores[0];
    
    for (int layer = 1; layer < num_layers; layer++) {
        if (layer_scores[layer] < min_score) {
            min_score = layer_scores[layer];
            best_layer = layer;
        }
    }
    
    VTR_LOGV_DEBUG(g_vpr_ctx.placement().f_placer_debug, 
                   "\tOptimal layer for critical block: %d (score: %.3f)\n", 
                   best_layer, min_score);
    
    return best_layer;
}