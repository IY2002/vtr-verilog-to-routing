#ifndef VPR_PARTITION_CREATOR_H
#define VPR_PARTITION_CREATOR_H

#include <string>
#include "atom_netlist.h"
#include "clustered_netlist.h"
#include "hyper_graph.h"
#include "vpr_context.h"
#include "prepack.h"
#include "tatum/TimingReporter.hpp"
#include "tatum/echo_writer.hpp"
#include "cluster_legalizer.h"
#include "pack_types.h"
#include "cluster_util.h"

class PartitionCreator {
  public:
    /**
     * @brief Create placement constraints from partitioning results
     * @param atom_netlist The atom netlist
     * @param num_partitions Number of partitions
     * @param imbalance_rate Maximum allowed imbalance
     * @param config_file KaHyPar configuration file (optional)
     */
    void create_partition_constraints(
                                    FloorplanningContext& floorplanning_ctx,
                                    DeviceContext& device_ctx,
                                    const AtomNetlist& atom_netlist,
                                    int num_partitions,
                                    float imbalance_rate,
                                    float cost_alpha = 0.0f,
                                    float cost_beta = 0.0f,
                                    const std::string& config_file = "");

    /**
     * @brief Create placement constraints from partitioning results after packing
     * @param cluster_netlist The clustered netlist
     * @param num_partitions Number of partitions
     * @param imbalance_rate Maximum allowed imbalance
     * @param cost_alpha Cost factor for hyperedge weights (default 0.0)
     * @param cost_beta Cost factor for vertex weights (default 0.0)
     * @param config_file KaHyPar configuration file (optional)
     */
    void create_partition_constraints_post_pack(
                                    FloorplanningContext& floorplanning_ctx,
                                    DeviceContext& device_ctx,
                                    const ClusteredNetlist& cluster_netlist,
                                    int num_partitions,
                                    float imbalance_rate,
                                    float cost_alpha,
                                    float cost_beta,
                                    const std::string& config_file, vtr::vector<AtomBlockId, float>& atom_criticality,vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup);

    /**
     * @brief Create placement constraints from partitioning results using Prepacker
     * @param atom_netlist The atom netlist
     * @param prepacker The Prepacker object containing preprocessed data
     * @param num_partitions Number of partitions
     * @param imbalance_rate Maximum allowed imbalance
     * @param cost_alpha Cost factor for hyperedge weights (default 0.0)
     * @param cost_beta Cost factor for vertex weights (default 0.0)
     * @param config_file KaHyPar configuration file (optional)
     * */
    void create_partition_constraints(
                                    FloorplanningContext& floorplanning_ctx,
                                    DeviceContext& device_ctx,
                                    const AtomNetlist& atom_netlist,
                                    const Prepacker& prepacker,
                                    int num_partitions,
                                    float imbalance_rate,
                                    t_packer_opts* packer_opts,
                                    const t_analysis_opts* analysis_opts,
                                    float cost_alpha = 0.0f,
                                    float cost_beta = 0.0f,
                                    const std::string& config_file = "");
};

#endif // VPR_PARTITION_CREATOR_H