#include "partition_creator.h"
#include "hyper_graph.h"
#include "partitioning_engine.h"
#include "user_place_constraints.h"
#include "partition.h"
#include "partition_region.h"
#include "region.h"
#include "vtr_log.h"
#include "vpr_types.h"
#include "PreClusterTimingGraphResolver.h"
#include "PreClusterDelayCalculator.h"
#include "concrete_timing_info.h"
#include "vtr_time.h"



void PartitionCreator::create_partition_constraints(FloorplanningContext& floorplanning_ctx,
                                                    DeviceContext& device_ctx,
                                                  const AtomNetlist& atom_netlist,
                                                  int num_partitions,
                                                  float imbalance_rate,
                                                  float cost_alpha,
                                                  float cost_beta,
                                                  const std::string& config_file) {
    
    // Create hypergraph from atom netlist
    AtomHypergraph hypergraph(atom_netlist);
    
    // Run partitioning
    AtomPartitioningEngine engine(hypergraph, num_partitions, imbalance_rate, cost_alpha, cost_beta, config_file);
    // engine.partition_with_triton();
    // engine.partition();
    
    // Create UserPlaceConstraints object
    UserPlaceConstraints constraints;
    
    // Create partitions
    for (int part_id = 0; part_id < num_partitions; part_id++) {
        // Create Partition object
        Partition partition;
        
        // Create PartitionRegion object
        PartitionRegion part_region;
        
        // Create vector of Regions (using magic numbers for now)
        std::vector<Region> regions;

        //Get Device dimensions
        int device_width = std::numeric_limits<int>::max();
        int device_height = std::numeric_limits<int>::max();

        VTR_LOG("Creating partition %d with dimensions %dx%d\n", part_id, device_width, device_height);

        // Using partition number as layer, magic numbers for coordinates
        Region region(0, 0, device_width, device_height, part_id);  // x_min, y_min, x_max, y_max, layer
        regions.push_back(region);
        
        // Set the regions in PartitionRegion
        part_region.set_partition_region(regions);
        
        // Set the PartitionRegion in Partition
        partition.set_part_region(part_region);
        
        // Add Partition to UserPlaceConstraints
        constraints.add_partition(partition);
    }
    
    // Add atoms to their respective partitions
    for (auto atom_id : atom_netlist.blocks()) {
        // Get partition assignment for this atom
        int assigned_partition = engine.get_block_partition(atom_id);
        
        if (assigned_partition >= 0) {
            // The partition index in UserPlaceConstraints is the same as assigned_partition
            // since we added them in order (0 to num_partitions-1)
            constraints.add_constrained_atom(atom_id, (PartitionId) assigned_partition);
        }
    }
    
    // Set the constraints in VPR context
    floorplanning_ctx.constraints = constraints;
}

void PartitionCreator::create_partition_constraints_post_pack(
                                    FloorplanningContext& floorplanning_ctx,
                                    DeviceContext& device_ctx,
                                    const ClusteredNetlist& cluster_netlist,
                                    int num_partitions,
                                    float imbalance_rate,
                                    float cost_alpha,
                                    float cost_beta,
                                    const std::string& config_file, vtr::vector<AtomBlockId, float>& atom_criticality,vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup){
                                        // Create Partition object
        // Create hypergraph from atom netlist
    vtr::ScopedStartFinishTimer timer("Post-Packing Partitioning");
    ClusterHypergraph hypergraph(cluster_netlist);
    
    // Run partitioning
    ClusterPartitioningEngine engine(hypergraph, num_partitions, imbalance_rate, cost_alpha, cost_beta, config_file);
    std::string triton_path = "/path/to/tritonpart_run_script.sh";
    engine.partition_with_triton(triton_path,1, atom_criticality, atoms_lookup);
    // engine.partition(atom_criticality, atoms_lookup);
    
    //Create Vector of Partition objects
    std::vector<Partition> partitions;
    std::vector<std::vector<Region>> partition_regions;
    // Create partitions
    for (int part_id = 0; part_id < num_partitions; part_id++) {
        Partition partition;
        
        // Create PartitionRegion object
        PartitionRegion part_region;
        
        // Create vector of Regions (using magic numbers for now)
        std::vector<Region> regions;

        //Set to max int, want to cover the whole layer
        int device_width = std::numeric_limits<int>::max();
        int device_height = std::numeric_limits<int>::max();

        VTR_LOG("Creating partition %d with dimensions %dx%d\n", part_id, device_width, device_height);

        // Using partition number as layer, magic numbers for coordinates
        Region region(0, 0, device_width, device_height, part_id);  // x_min, y_min, x_max, y_max, layer
        regions.push_back(region);

        // Add the region to the partition regions vector
        partition_regions.push_back(regions);
        
        // Set the regions in PartitionRegion
        part_region.set_partition_region(regions);
        
        // Set the PartitionRegion in Partition
        partition.set_part_region(part_region);
        
        // Add Partition to partitions vector
        partitions.push_back(partition);

        floorplanning_ctx.constraints.add_partition(partition);
    }
    

    for (auto& block : cluster_netlist.blocks()) {
        int assigned_partition = engine.get_block_partition(block);
        if (assigned_partition >= 0) {
            // Map the block to its partition
            auto& constrain_region = floorplanning_ctx.cluster_constraints[block].get_mutable_regions();
            constrain_region = partition_regions[assigned_partition];
        }
    }
    
    // // Set the constraints in VPR context
    // floorplanning_ctx.constraints = constraints;
                                    }

void create_timing_file(std::string filename,
                        const t_packer_opts* packer_opts,
                        const t_analysis_opts* analysis_opts,
                        const Prepacker& prepacker,
                        std::shared_ptr<PreClusterDelayCalculator>& clustering_delay_calc,
                        std::shared_ptr<SetupTimingInfo>& timing_info,
                        const std::map<t_pack_molecule*, int>& molecule_to_vertex) {
    


    const AtomContext& atom_ctx = g_vpr_ctx.atom();

    /*
     * Initialize the timing analyzer
     */
    clustering_delay_calc = std::make_shared<PreClusterDelayCalculator>(atom_ctx.nlist, atom_ctx.lookup, packer_opts->inter_cluster_net_delay, prepacker);
    timing_info = make_setup_timing_info(clustering_delay_calc, packer_opts->timing_update_type);

    //Calculate the initial timing
    timing_info->update();

    {
        auto& timing_ctx = g_vpr_ctx.timing();
        PreClusterTimingGraphResolver resolver(atom_ctx.nlist,
                                               atom_ctx.lookup, *timing_ctx.graph, *clustering_delay_calc);
        resolver.set_detail_level(analysis_opts->timing_report_detail);

        tatum::TimingReporter timing_reporter(resolver, *timing_ctx.graph,
                                              *timing_ctx.constraints);

        timing_reporter.create_timing_file("timing_file.txt", *timing_info->setup_analyzer(), 9999999, molecule_to_vertex, prepacker);
    }
}

void PartitionCreator::create_partition_constraints(FloorplanningContext& floorplanning_ctx,
                                                    DeviceContext& device_ctx,
                                                  const AtomNetlist& atom_netlist,
                                                  const Prepacker& prepacker,
                                                  int num_partitions,
                                                  float imbalance_rate,
                                                  t_packer_opts* packer_opts,
                                                  const t_analysis_opts* analysis_opts,
                                                  float cost_alpha,
                                                  float cost_beta,
                                                  const std::string& config_file) {
    
    // Create molecule-based hypergraph from atom netlist and prepacker
    AtomHypergraph hypergraph(atom_netlist);

    hypergraph.write_hmetis_format("atom_hypergraph.hgr", std::vector<int>(), std::vector<int>());
    // Rebuild hypergraph to include molecules
    std::map<t_pack_molecule*, int> molecule_to_vertex = build_molecule_hypergraph(hypergraph, atom_netlist, prepacker);

    // Create timing file
    std::shared_ptr<PreClusterDelayCalculator> clustering_delay_calc;
    std::shared_ptr<SetupTimingInfo> timing_info;
    create_timing_file("timing_file.txt", packer_opts, analysis_opts,
                        prepacker, clustering_delay_calc, timing_info, molecule_to_vertex);

    // Run partitioning with molecule support
    AtomPartitioningEngine engine(hypergraph, num_partitions, imbalance_rate, cost_alpha, cost_beta, config_file);
    // engine.partition_with_triton();
    // engine.partition();
    
    // Create UserPlaceConstraints object
    UserPlaceConstraints constraints;
    
    // Create partitions
    for (int part_id = 0; part_id < num_partitions; part_id++) {
        // Create Partition object
        Partition partition;
        
        // Create PartitionRegion object
        PartitionRegion part_region;
        
        // Create vector of Regions (using magic numbers for now)
        std::vector<Region> regions;

        //Get Device dimensions
        int device_width = std::numeric_limits<int>::max();
        int device_height = std::numeric_limits<int>::max();

        VTR_LOG("Creating partition %d with dimensions %dx%d\n", part_id, device_width, device_height);

        // Using partition number as layer, magic numbers for coordinates
        Region region(0, 0, device_width, device_height, part_id);  // x_min, y_min, x_max, y_max, layer
        regions.push_back(region);
        
        // Set the regions in PartitionRegion
        part_region.set_partition_region(regions);
        
        // Set the PartitionRegion in Partition
        partition.set_part_region(part_region);
        
        // Add Partition to UserPlaceConstraints
        constraints.add_partition(partition);
    }
    
    // Add atoms to their respective partitions
    // The partitioning engine handles the atom-to-molecule mapping internally
    for (auto atom_id : atom_netlist.blocks()) {
        // Get partition assignment for this atom
        int assigned_partition = engine.get_block_partition(atom_id);
        
        if (assigned_partition >= 0) {
            // The partition index in UserPlaceConstraints is the same as assigned_partition
            // since we added them in order (0 to num_partitions-1)
            constraints.add_constrained_atom(atom_id, (PartitionId) assigned_partition);
        }
    }
    
    // Set the constraints in VPR context
    floorplanning_ctx.constraints = constraints;
    
    VTR_LOG("Created molecule-based partition constraints\n");
}

