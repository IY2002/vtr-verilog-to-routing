
#include "placer.h"

#include <utility>

#include "vtr_time.h"
#include "draw.h"
#include "read_place.h"
#include "analytic_placer.h"
#include "initial_placement.h"
#include "concrete_timing_info.h"
#include "verify_placement.h"
#include "place_timing_update.h"
#include "annealer.h"
#include "RL_agent_util.h"
#include "place_checkpoint.h"
#include "tatum/echo_writer.hpp"
#include "VprTimingGraphResolver.h"
#include "tatum/TimingReporter.hpp"

Placer::Placer(const Netlist<>& net_list,
               const t_placer_opts& placer_opts,
               const t_analysis_opts& analysis_opts,
               const t_noc_opts& noc_opts,
               const IntraLbPbPinLookup& pb_gpin_lookup,
               const ClusteredPinAtomPinsLookup& netlist_pin_lookup,
               const std::vector<t_direct_inf>& directs,
               std::shared_ptr<PlaceDelayModel> place_delay_model,
               bool cube_bb,
               bool is_flat,
               bool quiet)
    : placer_opts_(placer_opts)
    , analysis_opts_(analysis_opts)
    , noc_opts_(noc_opts)
    , pb_gpin_lookup_(pb_gpin_lookup)
    , netlist_pin_lookup_(netlist_pin_lookup)
    , costs_(placer_opts.place_algorithm, noc_opts.noc)
    , placer_state_(placer_opts.place_algorithm.is_timing_driven(), cube_bb)
    , rng_(placer_opts.seed)
    , net_cost_handler_(placer_opts, placer_state_, cube_bb)
    , place_delay_model_(std::move(place_delay_model))
    , log_printer_(*this, quiet)
    , is_flat_(is_flat) {
    const auto& cluster_ctx = g_vpr_ctx.clustering();

    pre_place_timing_stats_ = g_vpr_ctx.timing().stats;

    init_placement_context(placer_state_.mutable_blk_loc_registry(), directs);

    // create a NoC cost handler if NoC optimization is enabled
    if (noc_opts.noc) {
        noc_cost_handler_.emplace(placer_state_.block_locs());
    }

    /* To make sure the importance of NoC-related cost terms compared to
     * BB and timing cost is determine only through NoC placement weighting factor,
     * we normalize NoC-related cost weighting factors so that they add up to 1.
     * With this normalization, NoC-related cost weighting factors only determine
     * the relative importance of NoC cost terms with respect to each other, while
     * the importance of total NoC cost to conventional placement cost is determined
     * by NoC placement weighting factor.
     */
    if (noc_opts.noc) {
        normalize_noc_cost_weighting_factor(const_cast<t_noc_opts&>(noc_opts));
    }

    BlkLocRegistry& blk_loc_registry = placer_state_.mutable_blk_loc_registry();
    initial_placement(placer_opts, placer_opts.constraints_file.c_str(),
                      noc_opts, blk_loc_registry, noc_cost_handler_, rng_);

    const int move_lim = (int)(placer_opts.anneal_sched.inner_num * pow(net_list.blocks().size(), 1.3333));
    //create the move generator based on the chosen placement strategy
    auto [move_generator, move_generator2] = create_move_generators(placer_state_, placer_opts, move_lim, noc_opts.noc_centroid_weight, rng_);

    if (!placer_opts.write_initial_place_file.empty()) {
        print_place(nullptr, nullptr, placer_opts.write_initial_place_file.c_str(), placer_state_.block_locs());
    }

#ifdef ENABLE_ANALYTIC_PLACE
    /*
     * Cluster-level Analytic Placer:
     *  Passes in the initial_placement via vpr_context, and passes its placement back via locations marked on
     *  both the clb_netlist and the gird.
     *  Most of anneal is disabled later by setting initial temperature to 0 and only further optimizes in quench
     */
    if (placer_opts.enable_analytic_placer) {
        AnalyticPlacer{blk_loc_registry}.ap_place();
    }

#endif /* ENABLE_ANALYTIC_PLACE */

    // Update physical pin values
   for (const ClusterBlockId block_id : cluster_ctx.clb_nlist.blocks()) {
       blk_loc_registry.place_sync_external_block_connections(block_id);
   }

   if (!quiet) {
#ifndef NO_GRAPHICS
       if (noc_cost_handler_.has_value()) {
           get_draw_state_vars()->set_noc_link_bandwidth_usages_ref(noc_cost_handler_->get_link_bandwidth_usages());
       }
#endif

       // width_fac gives the width of the widest channel
       const int width_fac = placer_opts.place_chan_width;
       init_draw_coords((float)width_fac, placer_state_.blk_loc_registry());
   }

   // Gets initial cost and loads bounding boxes.
   costs_.bb_cost = net_cost_handler_.comp_bb_cost(e_cost_methods::NORMAL);
   costs_.bb_cost_norm = 1 / costs_.bb_cost;

   if (placer_opts.place_algorithm.is_timing_driven()) {
       alloc_and_init_timing_objects_(net_list, analysis_opts);
   } else {
       VTR_ASSERT(placer_opts.place_algorithm == e_place_algorithm::BOUNDING_BOX_PLACE);
       // Timing cost and normalization factors are not used
       constexpr double INVALID_COST = std::numeric_limits<double>::quiet_NaN();
       costs_.timing_cost = INVALID_COST;
       costs_.timing_cost_norm = INVALID_COST;
   }

   if (noc_opts.noc) {
       VTR_ASSERT(noc_cost_handler_.has_value());

       // get the costs associated with the NoC
       costs_.noc_cost_terms.aggregate_bandwidth = noc_cost_handler_->comp_noc_aggregate_bandwidth_cost();
       std::tie(costs_.noc_cost_terms.latency, costs_.noc_cost_terms.latency_overrun) = noc_cost_handler_->comp_noc_latency_cost();
       costs_.noc_cost_terms.congestion = noc_cost_handler_->comp_noc_congestion_cost();

       // initialize all the noc normalization factors
       noc_cost_handler_->update_noc_normalization_factors(costs_);
   }

   // set the starting total placement cost
   costs_.cost = costs_.get_total_cost(placer_opts, noc_opts, placer_opts_.timing_tradeoff);

   // Sanity check that initial placement is legal
   check_place_();

   log_printer_.print_initial_placement_stats();

   annealer_ = std::make_unique<PlacementAnnealer>(placer_opts_, placer_state_, costs_, net_cost_handler_, noc_cost_handler_,
                                                   noc_opts_, rng_, std::move(move_generator), std::move(move_generator2), place_delay_model_.get(),
                                                   placer_criticalities_.get(), placer_setup_slacks_.get(), timing_info_.get(), pin_timing_invalidator_.get(),
                                                   move_lim);
}

void Placer::alloc_and_init_timing_objects_(const Netlist<>& net_list,
                                            const t_analysis_opts& analysis_opts) {
   const auto& atom_ctx = g_vpr_ctx.atom();
   const auto& cluster_ctx = g_vpr_ctx.clustering();
   const auto& timing_ctx = g_vpr_ctx.timing();
   const auto& p_timing_ctx = placer_state_.timing();

   // Update the point-to-point delays from the initial placement
   comp_td_connection_delays(place_delay_model_.get(), placer_state_);

   // Initialize timing analysis
   placement_delay_calc_ = std::make_shared<PlacementDelayCalculator>(atom_ctx.nlist,
                                                                      atom_ctx.lookup,
                                                                      p_timing_ctx.connection_delay,
                                                                      is_flat_);
   placement_delay_calc_->set_tsu_margin_relative(placer_opts_.tsu_rel_margin);
   placement_delay_calc_->set_tsu_margin_absolute(placer_opts_.tsu_abs_margin);

   timing_info_ = make_setup_timing_info(placement_delay_calc_, placer_opts_.timing_update_type);

   placer_setup_slacks_ = std::make_unique<PlacerSetupSlacks>(cluster_ctx.clb_nlist,
                                                              netlist_pin_lookup_,
                                                              timing_info_);

   placer_criticalities_ = std::make_unique<PlacerCriticalities>(cluster_ctx.clb_nlist,
                                                                 netlist_pin_lookup_,
                                                                 timing_info_);

   pin_timing_invalidator_ = make_net_pin_timing_invalidator(placer_opts_.timing_update_type,
                                                             net_list,
                                                             netlist_pin_lookup_,
                                                             atom_ctx.nlist,
                                                             atom_ctx.lookup,
                                                             timing_info_,
                                                             is_flat_);

   // First time compute timing and costs, compute from scratch
   PlaceCritParams crit_params;
   crit_params.crit_exponent = placer_opts_.td_place_exp_first;
   crit_params.crit_limit = placer_opts_.place_crit_limit;

   initialize_timing_info(crit_params, place_delay_model_.get(), placer_criticalities_.get(),
                          placer_setup_slacks_.get(), pin_timing_invalidator_.get(),
                          timing_info_.get(), &costs_, placer_state_);

   critical_path_ = timing_info_->least_slack_critical_path();

   // Write out the initial timing echo file
   if (isEchoFileEnabled(E_ECHO_INITIAL_PLACEMENT_TIMING_GRAPH)) {
       tatum::write_echo(getEchoFileName(E_ECHO_INITIAL_PLACEMENT_TIMING_GRAPH),
                         *timing_ctx.graph, *timing_ctx.constraints,
                         *placement_delay_calc_, timing_info_->analyzer());

       tatum::NodeId debug_tnode = id_or_pin_name_to_tnode(analysis_opts.echo_dot_timing_graph_node);

       write_setup_timing_graph_dot(getEchoFileName(E_ECHO_INITIAL_PLACEMENT_TIMING_GRAPH) + std::string(".dot"),
                                    *timing_info_, debug_tnode);
   }

   costs_.timing_cost_norm = 1 / costs_.timing_cost;
}

void Placer::check_place_() {
    const ClusteredNetlist& clb_nlist = g_vpr_ctx.clustering().clb_nlist;
    const DeviceGrid& device_grid = g_vpr_ctx.device().grid;
    const auto& cluster_constraints = g_vpr_ctx.floorplanning().cluster_constraints;
    
    int error = 0;
    
   // Verify the placement invariants independent to the placement flow.
   error += verify_placement(placer_state_.blk_loc_registry(),
                             clb_nlist,
                             device_grid,
                             cluster_constraints);

   error += check_placement_costs_();

   if (noc_opts_.noc) {
       // check the NoC costs during placement if the user is using the NoC supported flow
       error += noc_cost_handler_->check_noc_placement_costs(costs_, PL_INCREMENTAL_COST_TOLERANCE, noc_opts_);
       // make sure NoC routing configuration does not create any cycles in CDG
       error += (int)noc_cost_handler_->noc_routing_has_cycle();
   }

   if (error == 0) {
       VTR_LOG("\n");
       VTR_LOG("Completed placement consistency check successfully.\n");

   } else {
       VPR_ERROR(VPR_ERROR_PLACE,
                 "\nCompleted placement consistency check, %d errors found.\n"
                 "Aborting program.\n",
                 error);
   }
}

int Placer::check_placement_costs_() {
   int error = 0;
   double timing_cost_check;

   double bb_cost_check = net_cost_handler_.comp_bb_cost(e_cost_methods::CHECK);

   if (fabs(bb_cost_check - costs_.bb_cost) > costs_.bb_cost * PL_INCREMENTAL_COST_TOLERANCE) {
       VTR_LOG_ERROR(
           "bb_cost_check: %g and bb_cost: %g differ in check_place.\n",
           bb_cost_check, costs_.bb_cost);
       error++;
   }

   if (placer_opts_.place_algorithm.is_timing_driven()) {
       comp_td_costs(place_delay_model_.get(), *placer_criticalities_, placer_state_, &timing_cost_check);
       //VTR_LOG("timing_cost recomputed from scratch: %g\n", timing_cost_check);
       if (fabs(timing_cost_check - costs_.timing_cost) > costs_.timing_cost * PL_INCREMENTAL_COST_TOLERANCE) {
           VTR_LOG_ERROR(
               "timing_cost_check: %g and timing_cost: %g differ in check_place.\n",
               timing_cost_check, costs_.timing_cost);
           error++;
       }
   }
   return error;
}

void Placer::place() {
   const auto& timing_ctx = g_vpr_ctx.timing();
   const auto& cluster_ctx = g_vpr_ctx.clustering();


   bool skip_anneal = false;
#ifdef ENABLE_ANALYTIC_PLACE
   // Cluster-level analytic placer: when enabled, skip most of the annealing and go straight to quench
   if (placer_opts_.enable_analytic_placer) {
       skip_anneal = true;
   }
#endif

   // soft partitioning init_place
   if (placer_opts_.soft_partitioning == e_soft_partitioning::INIT_PLACE) {
        VTR_LOG("Soft partitioning: INIT_PLACE\n");
        // clear the partition constraints
        g_vpr_ctx.mutable_floorplanning().constraints = UserPlaceConstraints();

        for (ClusterBlockId cluster_id : cluster_ctx.clb_nlist.blocks()) {
            g_vpr_ctx.mutable_floorplanning().cluster_constraints[cluster_id].get_mutable_regions().clear();
        }
            
        g_vpr_ctx.mutable_floorplanning().compressed_cluster_constraints = std::vector<vtr::vector<ClusterBlockId, PartitionRegion>>();
   }


   bool cluster_cleared = false;

   // check place every 10 iterations
   int check_place_counter = 0;

   float lsmc_activation_threshold = placer_opts_.kick_move_activation_percentage / 100.0;

   if (!skip_anneal) {
       // Table header
       log_printer_.print_place_status_header();

       // Track consecutive kick moves without checkpoint saves
       int consecutive_kick_moves_without_save = 0;
       
       // How many outer loop iterations since last kick move (to avoid placement running forever)
       int iterations_since_last_kick_move = 0;

        const int success_rate_avg_window = 5;
        float success_rates[success_rate_avg_window] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        int success_rate_idx = 0;
        int success_rate_count = 0;
        float success_rate_sum = 5.0f; // since initialized to 1.0 each
        float success_rate_moving_avg = 1.0f;

        const float sr_start = placer_opts_.timing_layer_weight_start_sr;  // success rate where 3D influence begins
        const float sr_end   = placer_opts_.timing_layer_weight_end_sr;  // success rate where 3D influence is fully on
        const float layer_weight_start = placer_opts_.timing_layer_weight_start; // layer weight when 3D influence begins
        const float layer_weight_end   = placer_opts_.timing_layer_weight_end; // layer weight when 3D influence is fully on

        float layer_weight = layer_weight_start;
        // if (layer_weight_start > layer_weight_end){ layer_weight = }
        float timing_tradeoff = placer_opts_.timing_tradeoff_start;
       // Outer loop of the simulated annealing begins
       do {
           vtr::Timer temperature_timer;

           const auto& [swap_stats_pre_iter, move_type_stats_pre_iter, placer_stats_pre_iter] = annealer_->get_stats();
           float success_rate = placer_stats_pre_iter.success_rate;

           if (success_rate == 0.0) { //This means this is the first iteration
               success_rate = 1.0;
            }

            

            // update sliding window average
            success_rate_sum -= success_rates[success_rate_idx];
            success_rate_sum += success_rate;
            success_rates[success_rate_idx] = success_rate;

            success_rate_idx = (success_rate_idx + 1) % success_rate_avg_window;
            if (success_rate_count < success_rate_avg_window)
                success_rate_count++;

            success_rate_moving_avg = success_rate_sum / success_rate_count;


            


            float new_layer_weight = layer_weight;

            if (placer_opts_.timing_layer_weight_adjustor == e_timing_layer_weight_adjustor::NONE) {
                new_layer_weight = 1.0;
            } else if (placer_opts_.timing_layer_weight_adjustor == e_timing_layer_weight_adjustor::STEP){
                if (success_rate_moving_avg >= sr_start) new_layer_weight = 0.0;  // still flat 2D
                else new_layer_weight = 1.0;  // fully 3D
                new_layer_weight = layer_weight_start + new_layer_weight * (layer_weight_end - layer_weight_start);
            } else if (placer_opts_.timing_layer_weight_adjustor == e_timing_layer_weight_adjustor::LINEAR){
                if (success_rate_moving_avg >= sr_start) new_layer_weight = 0.0;  // still flat 2D
                else if (success_rate_moving_avg <= sr_end) new_layer_weight = 1.0;  // fully 3D
                else new_layer_weight = ((sr_start - success_rate_moving_avg) / (sr_start - sr_end));
                new_layer_weight = layer_weight_start + new_layer_weight * (layer_weight_end - layer_weight_start);
            } else if (placer_opts_.timing_layer_weight_adjustor == e_timing_layer_weight_adjustor::QUADRATIC){
                if (success_rate_moving_avg >= sr_start) new_layer_weight = 0.0;  // still flat 2D
                else if (success_rate_moving_avg <= sr_end) new_layer_weight = 1.0;  // fully 3D
                else new_layer_weight = ((sr_start - success_rate_moving_avg) / (sr_start - sr_end)) * ((sr_start - success_rate_moving_avg) / (sr_start - sr_end));
                new_layer_weight = layer_weight_start + new_layer_weight * (layer_weight_end - layer_weight_start);
            }

            // if (new_layer_weight > layer_weight) { // only update the weight if it is increasing, reduces the sensitivity to noise in success rate
            //     layer_weight = new_layer_weight;
            // }

            

            if (layer_weight_start < layer_weight_end) {
                // increasing layer weight
                layer_weight = std::max(layer_weight, new_layer_weight);
            } else {
                // decreasing layer weight
                layer_weight = std::min(layer_weight, new_layer_weight);
            }

            if (placer_opts_.timing_layer_weight_adjustor == e_timing_layer_weight_adjustor::NONE) {
                layer_weight = 1.0; // sanity check
            }

            // VTR_LOG("New layer weight: %.2f, layer_weight: %.2f\n", new_layer_weight, layer_weight);

            float tt_start = placer_opts_.timing_tradeoff_start; // WL-focused
            float tt_end   = placer_opts_.timing_tradeoff_end; // timing-focused
            float tt_start_sr = placer_opts_.timing_tradeoff_start_sr; // WL-focused
            float tt_end_sr   = placer_opts_.timing_tradeoff_end_sr; // timing-focused
            float new_timing_tradeoff = timing_tradeoff;
            if (placer_opts_.timing_tradeoff_adjustor == e_timing_tradeoff_adjustor::NONE) {
                // do nothing, use the user-specified constant timing tradeoff
                new_timing_tradeoff = placer_opts_.timing_tradeoff;
            } else if (placer_opts_.timing_tradeoff_adjustor == e_timing_tradeoff_adjustor::LINEAR) {
                // Linear ramp based on success rate
                if (success_rate_moving_avg >= tt_start_sr) new_timing_tradeoff = tt_start;
                else if (success_rate_moving_avg <= tt_end_sr) new_timing_tradeoff = tt_end;
                else new_timing_tradeoff = tt_end - (success_rate_moving_avg * (tt_end - tt_start));

            } else if (placer_opts_.timing_tradeoff_adjustor == e_timing_tradeoff_adjustor::QUADRATIC) {
                // Quadratic ramp based on success rate
                if (success_rate_moving_avg >= tt_start_sr) new_timing_tradeoff = tt_start;
                else if (success_rate_moving_avg <= tt_end_sr) new_timing_tradeoff = tt_end;
                else new_timing_tradeoff = tt_end - (success_rate_moving_avg * success_rate_moving_avg * (tt_end - tt_start));
            } else {
                VPR_FATAL_ERROR(VPR_ERROR_PLACE,
                                "Unknown timing tradeoff adjustor type %d.\n",
                                static_cast<int>(placer_opts_.timing_tradeoff_adjustor));
            }
            // if (new_timing_tradeoff > timing_tradeoff) { // only update the timing tradeoff if it is decreasing, reduces the sensitivity to noise in success rate
            //     timing_tradeoff = new_timing_tradeoff;
            // }

            if (tt_start < tt_end) {
                // increasing timing tradeoff
                timing_tradeoff = std::max(timing_tradeoff, new_timing_tradeoff);
            } else {
                // decreasing timing tradeoff
                timing_tradeoff = std::min(timing_tradeoff, new_timing_tradeoff);
            }

           annealer_->outer_loop_update_timing_info(layer_weight, timing_tradeoff);

           

           if (placer_opts_.place_algorithm.is_timing_driven()) {
               critical_path_ = timing_info_->least_slack_critical_path();

               // see if we should save the current placement solution as a checkpoint
               if (placer_opts_.place_checkpointing && annealer_->get_agent_state() == e_agent_state::LATE_IN_THE_ANNEAL && success_rate < 0.15) { // only save checkpoint in the late anneal stage and when success rate is low
                   
                   save_placement_checkpoint_if_needed(placer_state_.mutable_block_locs(),
                                                       placement_checkpoint_,
                                                       timing_info_, costs_, critical_path_.delay());
               }
           }

        //    check_place_counter++;
        //     if (check_place_counter >= 10) {
        //         check_place_counter = 0;
        //         check_place_();
        //     }
        //    if (annealer_->get_agent_state() == e_agent_state::LATE_IN_THE_ANNEAL && !cluster_cleared) {
        //        // clear the cluster constraints
        //        g_vpr_ctx.mutable_floorplanning().constraints = UserPlaceConstraints();
        //        for (ClusterBlockId cluster_id : cluster_ctx.clb_nlist.blocks()) {
        //            g_vpr_ctx.mutable_floorplanning().cluster_constraints[cluster_id].get_mutable_regions().clear();
        //        }
        //        g_vpr_ctx.mutable_floorplanning().compressed_cluster_constraints.clear();
        //        cluster_cleared = true;
        //    }

           

        //    VTR_LOG("Placement success rate: %.2f%% (moving avg: %.2f%%), Layer weight: %.2f, Timing tradeoff: %.2f\n",
        //            success_rate * 100, success_rate_moving_avg * 100, layer_weight, timing_tradeoff);

            

           // do a complete inner loop iteration
           bool kick_move = annealer_->placement_inner_loop(layer_weight, timing_tradeoff);

           bool should_check_for_kick = kick_move;
           const auto& [swap_stats, move_type_stats, placer_stats] = annealer_->get_stats();
           success_rate = placer_stats.success_rate;
           
           if (placer_opts_.soft_partitioning == e_soft_partitioning::MID_PLACE && !cluster_cleared && success_rate < placer_opts_.mid_soft_partitioning_enable_percent){
            VTR_LOG("Soft partitioning: MID_PLACE\n");
            // clear the partition constraints
            g_vpr_ctx.mutable_floorplanning().constraints = UserPlaceConstraints();

            for (ClusterBlockId cluster_id : cluster_ctx.clb_nlist.blocks()) {
                g_vpr_ctx.mutable_floorplanning().cluster_constraints[cluster_id].get_mutable_regions().clear();
            }
            g_vpr_ctx.mutable_floorplanning().compressed_cluster_constraints = std::vector<vtr::vector<ClusterBlockId, PartitionRegion>>();
            cluster_cleared = true;
           }

           if (placer_opts_.enable_kick_move && success_rate < lsmc_activation_threshold && consecutive_kick_moves_without_save < placer_opts_.kick_move_num) {
                should_check_for_kick = true;
                VTR_LOG("LSMC kick move: Activation threshold met (Success rate: %.1f%% < %.1f%%).\n", success_rate*100, lsmc_activation_threshold*100);
            } else if (placer_opts_.enable_kick_move && success_rate < lsmc_activation_threshold && consecutive_kick_moves_without_save == placer_opts_.kick_move_num) {
                // Restore the best placement solution so far
                const t_annealing_state& annealing_state = annealer_->get_annealing_state();
                PlaceCritParams crit_params;
                crit_params.crit_exponent = annealing_state.crit_exponent;
                crit_params.crit_limit = placer_opts_.place_crit_limit;
                save_placement_checkpoint_if_needed(placer_state_.mutable_block_locs(),
                                                       placement_checkpoint_,
                                                       timing_info_, costs_, critical_path_.delay());
                restore_best_placement(placer_state_,
                                   placement_checkpoint_, timing_info_, costs_,
                                   placer_criticalities_, placer_setup_slacks_, place_delay_model_,
                                   pin_timing_invalidator_, crit_params, noc_cost_handler_);

                // do a complete inner loop cost update
                annealer_->update_costs_after_kick_move(layer_weight, timing_tradeoff);

                consecutive_kick_moves_without_save++;
                VTR_LOG("LSMC kick move: restoring best placement solution found during search.\n");
            }   

           iterations_since_last_kick_move++;
           
           if (should_check_for_kick && annealer_->get_agent_state() != e_agent_state::LATE_IN_THE_ANNEAL && placer_opts_.enable_kick_move && consecutive_kick_moves_without_save < placer_opts_.kick_move_num) {
               // do a kick move to help the annealer escape local minima
               //first checkpoint the placement
               bool checkpoint_saved = false;
               
               if (placer_opts_.kick_move_checkpointing != e_kick_move_checkpointing::NONE) {
                   auto critical_path_ = timing_info_->least_slack_critical_path();
                   checkpoint_saved = save_placement_checkpoint_if_needed(placer_state_.mutable_block_locs(),
                                                       placement_checkpoint_,
                                                       timing_info_, costs_, critical_path_.delay());
                   VTR_LOG("LSMC kick move: %s placement checkpoint before the kick move.\n", 
                           checkpoint_saved ? "saved" : "did not save");
               }
               
               iterations_since_last_kick_move = 0;

               // Update the counter based on whether checkpoint was saved
               if (checkpoint_saved) {
                   if (consecutive_kick_moves_without_save > 0) {
                       consecutive_kick_moves_without_save--;  // Reset counter
                   }
               } else {
                   consecutive_kick_moves_without_save++;
               }
               
               VTR_LOG("LSMC kick move: Consecutive kick moves without checkpoint save: %d\n", 
                       consecutive_kick_moves_without_save);

               if (placer_opts_.kick_move_checkpointing == e_kick_move_checkpointing::PROGRESSIVE) {
                   // if we are in the progressive checkpointing mode, we should restore the best placement solution before the kick move
                   // restore the best placement solution
                   const t_annealing_state& annealing_state = annealer_->get_annealing_state();
                   PlaceCritParams crit_params;
                   crit_params.crit_exponent = annealing_state.crit_exponent;
                   crit_params.crit_limit = placer_opts_.place_crit_limit;
                   restore_best_placement(placer_state_,
                                   placement_checkpoint_, timing_info_, costs_,
                                   placer_criticalities_, placer_setup_slacks_, place_delay_model_,
                                   pin_timing_invalidator_, crit_params, noc_cost_handler_);
                   // do a complete inner loop cost update
                   annealer_->update_costs_after_kick_move(layer_weight, timing_tradeoff);
                   VTR_LOG("LSMC kick move: restoring best placement solution before the kick move.\n");
               }

               annealer_->lsmc_kick_move(placer_opts_.kick_move_percent_to_swap, layer_weight);
           } else if (success_rate < placer_opts_.rl_second_state_activation_percent && annealer_->get_agent_state() != e_agent_state::LATE_IN_THE_ANNEAL) { // means annealer should entered LATE_IN_THE_ANNEAL state
                   // if we are in the late annealing state, we should restore the best placement solution
                   // restore the best placement solution
                //    VTR_LOG("LSMC kick move: restoring best placement solution before the late annealing state.\n");
                   const t_annealing_state& annealing_state = annealer_->get_annealing_state();
                   PlaceCritParams crit_params;
                   crit_params.crit_exponent = annealing_state.crit_exponent;
                   crit_params.crit_limit = placer_opts_.place_crit_limit;
                //    restore_best_placement(placer_state_,
                //                  placement_checkpoint_, timing_info_, costs_,
                //                  placer_criticalities_, placer_setup_slacks_, place_delay_model_,
                //                  pin_timing_invalidator_, crit_params, noc_cost_handler_);

                    iterations_since_last_kick_move = 0;

                   // do a complete inner loop cost update
                   annealer_->update_costs_after_kick_move(layer_weight, timing_tradeoff);

                   annealer_->set_agent_state(e_agent_state::LATE_IN_THE_ANNEAL);

                   if (placer_opts_.soft_partitioning == e_soft_partitioning::LATE_PLACE) {
                       VTR_LOG("Soft partitioning: LATE_PLACE\n");
                       // clear the partition constraints
                       g_vpr_ctx.mutable_floorplanning().constraints = UserPlaceConstraints();

                       for (ClusterBlockId cluster_id : cluster_ctx.clb_nlist.blocks()) {
                           g_vpr_ctx.mutable_floorplanning().cluster_constraints[cluster_id].get_mutable_regions().clear();
                       }
                           
                       g_vpr_ctx.mutable_floorplanning().compressed_cluster_constraints = std::vector<vtr::vector<ClusterBlockId, PartitionRegion>>();
                   }
            }
           

           log_printer_.print_place_status(temperature_timer.elapsed_sec());

           // Outer loop of the simulated annealing ends
       } while (annealer_->outer_loop_update_state() || (placer_opts_.enable_kick_move && annealer_->get_agent_state() != e_agent_state::LATE_IN_THE_ANNEAL && iterations_since_last_kick_move < 1000)); // exit only if the annealing state is LATE_IN_THE_ANNEAL or 1000 iterations have passed when kick move is enabled
   } //skip_anneal ends

    // Start Quench
    annealer_->start_quench();

    pre_quench_timing_stats_ = timing_ctx.stats;
    { // Quench
       vtr::ScopedFinishTimer temperature_timer("Placement Quench");

       annealer_->outer_loop_update_timing_info(1.0, placer_opts_.timing_tradeoff);

       /* Run inner loop again with temperature = 0 so as to accept only swaps
        * which reduce the cost of the placement */
       annealer_->placement_inner_loop();

       if (placer_opts_.place_quench_algorithm.is_timing_driven()) {
           critical_path_ = timing_info_->least_slack_critical_path();

            
            

       }

       log_printer_.print_place_status(temperature_timer.elapsed_sec());
    }
    post_quench_timing_stats_ = timing_ctx.stats;

    // Final timing analysis
    const t_annealing_state& annealing_state = annealer_->get_annealing_state();
    PlaceCritParams crit_params;
    crit_params.crit_exponent = annealing_state.crit_exponent;
    crit_params.crit_limit = placer_opts_.place_crit_limit;

    if (placer_opts_.place_algorithm.is_timing_driven()) {
       perform_full_timing_update(crit_params, place_delay_model_.get(), placer_criticalities_.get(),
                                  placer_setup_slacks_.get(), pin_timing_invalidator_.get(),
                                  timing_info_.get(), &costs_, placer_state_);

       critical_path_ = timing_info_->least_slack_critical_path();

       VTR_LOG("post-quench CPD = %g (ns) \n",
               1e9 * critical_path_.delay());
    }

    // See if our latest checkpoint is better than the current placement solution
    if (placer_opts_.place_checkpointing) {
       restore_best_placement(placer_state_,
                              placement_checkpoint_, timing_info_, costs_,
                              placer_criticalities_, placer_setup_slacks_, place_delay_model_,
                              pin_timing_invalidator_, crit_params, noc_cost_handler_);
    }

    if (placer_opts_.placement_saves_per_temperature >= 1) {
       std::string filename = vtr::string_fmt("placement_%03d_%03d.place",
                                              annealing_state.num_temps + 1, 0);
       VTR_LOG("Saving final placement to file: %s\n", filename.c_str());
       print_place(nullptr, nullptr, filename.c_str(), placer_state_.mutable_block_locs());
    }

    
    // Update physical pin values
    for (const ClusterBlockId block_id : cluster_ctx.clb_nlist.blocks()) {
        placer_state_.mutable_blk_loc_registry().place_sync_external_block_connections(block_id);
    }
    
    check_place_();
    
    const auto& atom_ctx = g_vpr_ctx.atom();

    //setup timing resolver to resolve the timing graph
    VprTimingGraphResolver resolver(atom_ctx.nlist,
                                        atom_ctx.lookup, *timing_ctx.graph, *placement_delay_calc_, is_flat_, placer_state_.blk_loc_registry());
    resolver.set_detail_level(e_timing_report_detail::DETAILED_ROUTING);

    tatum::TimingReporter timing_reporter(resolver, *timing_ctx.graph,
                                        *timing_ctx.constraints);

    timing_reporter.report_timing_setup(std::string("post_quench_setup_timing.txt"),
                                        *timing_info_->setup_analyzer(),
                                        (size_t) 1000);
                                        
    log_printer_.print_post_placement_stats();
}

void Placer::copy_locs_to_global_state(PlacementContext& place_ctx) {
    // the placement location variables should be unlocked before being accessed
    place_ctx.unlock_loc_vars();

    // copy the local location variables into the global state
    auto& global_blk_loc_registry = place_ctx.mutable_blk_loc_registry();
    global_blk_loc_registry = placer_state_.blk_loc_registry();

#ifndef NO_GRAPHICS
    // update the graphics' reference to placement location variables
    get_draw_state_vars()->set_graphics_blk_loc_registry_ref(global_blk_loc_registry);
#endif
}
