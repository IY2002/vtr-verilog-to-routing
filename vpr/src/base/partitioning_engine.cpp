#include "partitioning_engine.h"
#include "vtr_log.h"
#include "vtr_assert.h"
#include <fstream>
#include <algorithm>
#include <map>
#include <chrono>

// Constructor
template<typename HypergraphType>
PartitioningEngine<HypergraphType>::PartitioningEngine(
    const HypergraphType& hypergraph,
    int num_partitions,
    float imbalance_rate,
    float cost_alpha,
    float cost_beta,
    const std::string& config_file)
    : hypergraph_(hypergraph)
    , num_partitions_(num_partitions)
    , imbalance_rate_(imbalance_rate)
    , cost_alpha_(cost_alpha)
    , cost_beta_(cost_beta)
    , config_file_(config_file)
    , context_(nullptr)
    , objective_(0)
    , is_partitioned_(false),
      is_molecule_based_(false) {
    
    VTR_ASSERT(num_partitions >= 2);
    VTR_LOG("Initializing partitioning engine with %d partitions, imbalance rate %.2f, alpha %.2f, beta %.2f\n", num_partitions, imbalance_rate_, cost_alpha_, cost_beta_);
    // VTR_ASSERT(imbalance_rate >= 0.0 && imbalance_rate <= 1.0);
    
    // Initialize KaHyPar context
    // context_ = kahypar_context_new();
    
    // if (!config_file_.empty()) {
    //     kahypar_configure_context_from_file(context_, config_file_.c_str());
    // }
    
    // kahypar_set_seed(context_, 42);  // For reproducibility
    
    // Initialize partition vector
    partition_.resize(hypergraph_.num_vertices(), -1);
}

// Destructor
template<typename HypergraphType>
PartitioningEngine<HypergraphType>::~PartitioningEngine() {
    // if (context_) {
    //     kahypar_context_free(context_);
    // }
}

// Run partitioning
template<typename HypergraphType>
bool PartitioningEngine<HypergraphType>::partition(vtr::vector<AtomBlockId, float>& atom_criticality,vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup) {
    VTR_LOG("\nRunning KaHyPar partitioning...\n");
    VTR_LOG("  Vertices: %d\n", hypergraph_.num_vertices());
    VTR_LOG("  Hyperedges: %d\n", hypergraph_.num_hyperedges());
    VTR_LOG("  Partitions: %d\n", num_partitions_);
    VTR_LOG("  Max imbalance: %.1f%%\n", imbalance_rate_ * 100);
    VTR_LOG("  Cost alpha: %.1f\n", cost_alpha_);
    VTR_LOG("  Cost beta: %.1f\n", cost_beta_);

    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Prepare data for KaHyPar
    const auto& indices = hypergraph_.hyperedge_indices();
    const auto& edges = hypergraph_.hyperedges();

    const auto hyperedge_weights = hypergraph_.get_edge_weights(cost_alpha_, atom_criticality, atoms_lookup);
    // const auto hyperedge_weights = nullptr;  // Use default weights if not provided

    // std::shared_ptr<PreClusterDelayCalculator> clustering_delay_calc;
    // std::shared_ptr<SetupTimingInfo> timing_info;



    const auto vertex_weights = hypergraph_.get_vertices_weights(cost_beta_, atom_criticality, atoms_lookup);

    // TODO: Add weighted edges based on criticality

    // Call KaHyPar
    // kahypar_partition(hypergraph_.num_vertices(),
    //                  hypergraph_.num_hyperedges(),
    //                  imbalance_rate_,
    //                  num_partitions_,
    //                  vertex_weights.data(), // Vertex weights
    //                  hyperedge_weights.data(), // Hyperedge weights
    //                  indices.data(),
    //                  edges.data(),
    //                  &objective_,
    //                  context_,
    //                  partition_.data());
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    is_partitioned_ = true;
    
    VTR_LOG("Partitioning completed in %.2f seconds\n", duration.count() / 1000.0);
    VTR_LOG("  Objective (cut): %d\n", objective_);

    // save intermediate partition file
    save_partition_file("kahypar_intermediate.part");

    return true;
}

// Get partition for vertex
template<typename HypergraphType>
int PartitioningEngine<HypergraphType>::get_partition(int vertex_id) const {
    VTR_ASSERT(is_partitioned_);
    VTR_ASSERT(vertex_id >= 0 && vertex_id < hypergraph_.num_vertices());
    return partition_[vertex_id];
}

// Get partition for block
template<typename HypergraphType>
int PartitioningEngine<HypergraphType>::get_block_partition(BlockIdType block_id) const {
    VTR_ASSERT(is_partitioned_);
    
    if (is_molecule_based_) {
        // This should only happen for AtomBlockId
        auto it = atom_to_molecule_.find(block_id);
        if (it != atom_to_molecule_.end()) {
            int molecule_vertex_id = it->second;
            if (molecule_vertex_id >= 0 && molecule_vertex_id < hypergraph_.num_vertices()) {
                return partition_[molecule_vertex_id];
            }
        }
        return -1;
    } else {
        int vertex_id = hypergraph_.block_to_vertex(block_id);
        if (vertex_id >= 0) {
            return get_partition(vertex_id);
        }
        return -1;
    }
}

// Print statistics
template<typename HypergraphType>
void PartitioningEngine<HypergraphType>::print_stats() const {
    if (!is_partitioned_) {
        VTR_LOG("No partition computed yet\n");
        return;
    }
    
    VTR_LOG("\nPartitioning Statistics:\n");
    VTR_LOG("=======================\n");
    VTR_LOG("Objective (hyperedge cut): %d\n", objective_);
    
    // Calculate partition sizes
    std::vector<int> partition_sizes(num_partitions_, 0);
    for (int p : partition_) {
        partition_sizes[p]++;
    }
    
    // Print partition sizes and balance
    VTR_LOG("\nPartition sizes:\n");
    double avg_size = static_cast<double>(hypergraph_.num_vertices()) / num_partitions_;
    for (int i = 0; i < num_partitions_; i++) {
        double imbalance = (partition_sizes[i] - avg_size) / avg_size;
        VTR_LOG("  Partition %d: %d vertices (%.1f%% from average)\n",
                i, partition_sizes[i], imbalance * 100);
    }
    
    // Calculate cut hyperedges
    int cut_hyperedges = 0;
    for (int he = 0; he < hypergraph_.num_hyperedges(); he++) {
        size_t start = hypergraph_.hyperedge_indices()[he];
        size_t end = hypergraph_.hyperedge_indices()[he + 1];
        
        // Check if hyperedge is cut
        int first_partition = partition_[hypergraph_.hyperedges()[start]];
        bool is_cut = false;
        for (size_t i = start + 1; i < end; i++) {
            if (partition_[hypergraph_.hyperedges()[i]] != first_partition) {
                is_cut = true;
                break;
            }
        }
        if (is_cut) cut_hyperedges++;
    }
    
    VTR_LOG("\nCut hyperedges: %d / %d (%.1f%%)\n",
            cut_hyperedges, hypergraph_.num_hyperedges(),
            100.0 * cut_hyperedges / hypergraph_.num_hyperedges());
}

// Save partition assignment
template<typename HypergraphType>
void PartitioningEngine<HypergraphType>::save_partition_file(const std::string& filename) const {
    if (!is_partitioned_) {
        VTR_LOG_ERROR("Cannot save partition - no partition computed\n");
        return;
    }
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        VTR_LOG_ERROR("Failed to open file '%s' for writing\n", filename.c_str());
        return;
    }
    
    // Write partition assignment for each vertex
    for (int v = 0; v < hypergraph_.num_vertices(); v++) {
        out << partition_[v] << "\n";
    }
    
    out.close();
    VTR_LOG("Saved partition assignment to '%s'\n", filename.c_str());
}

// Save detailed report
template<typename HypergraphType>
void PartitioningEngine<HypergraphType>::save_partition_report(const std::string& filename) const {
    if (!is_partitioned_) {
        VTR_LOG_ERROR("Cannot save report - no partition computed\n");
        return;
    }
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        VTR_LOG_ERROR("Failed to open file '%s' for writing\n", filename.c_str());
        return;
    }
    
    out << "KaHyPar Partitioning Report\n";
    out << "===========================\n\n";
    
    out << "Input:\n";
    out << "  Vertices: " << hypergraph_.num_vertices() << "\n";
    out << "  Hyperedges: " << hypergraph_.num_hyperedges() << "\n";
    out << "  Partitions: " << num_partitions_ << "\n";
    out << "  Max imbalance: " << (imbalance_rate_ * 100) << "%\n\n";
    
    out << "Results:\n";
    out << "  Objective (cut): " << objective_ << "\n\n";
    
    // Partition sizes
    std::vector<int> partition_sizes(num_partitions_, 0);
    for (int p : partition_) {
        partition_sizes[p]++;
    }
    
    out << "Partition sizes:\n";
    for (int i = 0; i < num_partitions_; i++) {
        out << "  Partition " << i << ": " << partition_sizes[i] << " vertices\n";
    }
    
    out.close();
    VTR_LOG("Saved partitioning report to '%s'\n", filename.c_str());
}

// Add Triton partitioning method
template<typename HypergraphType>
bool PartitioningEngine<HypergraphType>::partition_with_triton(
    const std::string& triton_script_path, 
    int seed, vtr::vector<AtomBlockId, float>& atom_criticality ,vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup ) {
    
    VTR_LOG("\nRunning Triton partitioning...\n");
    VTR_LOG("  Vertices: %d\n", hypergraph_.num_vertices());
    VTR_LOG("  Hyperedges: %d\n", hypergraph_.num_hyperedges());
    VTR_LOG("  Partitions: %d\n", num_partitions_);
    VTR_LOG("  Balance constraint: %.2f\n", imbalance_rate_);
    VTR_LOG("  Seed: %d\n", seed);

    auto start_time = std::chrono::high_resolution_clock::now();
    
    auto edge_weights = hypergraph_.get_edge_weights(cost_alpha_, atom_criticality, atoms_lookup);
    auto vertex_weights = hypergraph_.get_vertices_weights(cost_beta_, atom_criticality, atoms_lookup);


    // Write hypergraph to file in hMETIS format
    std::string hypergraph_file = "triton_hypergraph.hgr";
    hypergraph_.write_hmetis_format(hypergraph_file, vertex_weights, edge_weights);
    
    // Run Triton partitioner
    bool success = run_triton_script(hypergraph_file, triton_script_path, seed);
    
    if (success) {
        // Read partition assignment from output file
        // Triton typically outputs to filename.part.k format
        std::string partition_file = hypergraph_file + ".part." + std::to_string(num_partitions_);
        success = parse_triton_partition_output(partition_file);
        
        if (success) {
            // Calculate objective (cut size)
            objective_ = calculate_cut_objective();
            is_partitioned_ = true;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    if (success) {
        VTR_LOG("Triton partitioning completed in %.2f seconds\n", duration.count() / 1000.0);
        VTR_LOG("  Objective (cut): %d\n", objective_);
    } else {
        VTR_LOG_ERROR("Triton partitioning failed\n");
    }
    
    // Clean up temporary files
    std::remove(hypergraph_file.c_str());
    std::remove("run_triton_part.tcl");
    std::remove("timing_file.txt");
    
    return success;
}

// Run Triton script
template<typename HypergraphType>
bool PartitioningEngine<HypergraphType>::run_triton_script(
    const std::string& hypergraph_file,
    const std::string& triton_script_path,
    int seed) {
    
    // Build command
    std::stringstream cmd;
    cmd << triton_script_path << " ";
    cmd << hypergraph_file << " ";
    cmd << num_partitions_ << " ";
    cmd << imbalance_rate_ << " ";  // Using imbalance_rate as balance_constraint
    cmd << seed << " ";
    // cmd << "timing_file.txt";  // Assuming timing file is required
    
    VTR_LOG("Executing: %s\n", cmd.str().c_str());
    
    // Execute command
    int result = std::system(cmd.str().c_str());
    
    if (result != 0) {
        VTR_LOG_ERROR("Triton partitioner execution failed with code %d\n", result);
        return false;
    }
    
    
    return true;
}

// Parse Triton partition output
template<typename HypergraphType>
bool PartitioningEngine<HypergraphType>::parse_triton_partition_output(
    const std::string& partition_file) {
    
    std::ifstream in(partition_file);
    if (!in.is_open()) {
        VTR_LOG_ERROR("Failed to open Triton partition file '%s'\n", partition_file.c_str());
        return false;
    }
    
    // Read partition assignment for each vertex
    int vertex_id = 0;
    std::string line;
    while (std::getline(in, line) && vertex_id < hypergraph_.num_vertices()) {
        try {
            int part = std::stoi(line);
            if (part < 0 || part >= num_partitions_) {
                VTR_LOG_ERROR("Invalid partition ID %d for vertex %d\n", part, vertex_id);
                return false;
            }
            partition_[vertex_id] = part;
            vertex_id++;
        } catch (const std::exception& e) {
            VTR_LOG_ERROR("Error parsing partition file at line %d: %s\n", vertex_id + 1, e.what());
            return false;
        }
    }
    
    if (vertex_id != hypergraph_.num_vertices()) {
        VTR_LOG_ERROR("Partition file has %d entries but hypergraph has %d vertices\n", 
                      vertex_id, hypergraph_.num_vertices());
        return false;
    }
    
    in.close();
    return true;
}

// Calculate cut objective
template<typename HypergraphType>
int PartitioningEngine<HypergraphType>::calculate_cut_objective() const {
    int cut_hyperedges = 0;
    
    for (int he = 0; he < hypergraph_.num_hyperedges(); he++) {
        size_t start = hypergraph_.hyperedge_indices()[he];
        size_t end = hypergraph_.hyperedge_indices()[he + 1];
        
        // Check if hyperedge is cut
        int first_partition = partition_[hypergraph_.hyperedges()[start]];
        bool is_cut = false;
        
        for (size_t i = start + 1; i < end; i++) {
            if (partition_[hypergraph_.hyperedges()[i]] != first_partition) {
                is_cut = true;
                break;
            }
        }
        
        if (is_cut) {
            cut_hyperedges++;
        }
    }
    
    return cut_hyperedges;
}

// Explicit instantiations
template class PartitioningEngine<AtomHypergraph>;
template class PartitioningEngine<ClusterHypergraph>;
