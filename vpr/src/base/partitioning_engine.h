#ifndef VPR_PARTITIONING_ENGINE_H
#define VPR_PARTITIONING_ENGINE_H

#include <vector>
#include <string>
#include "hyper_graph.h"
#include "prepack.h"
// extern "C" {
//     #include "libkahypar.h"
// }

/**
 * @brief Engine for partitioning VPR hypergraphs using KaHyPar
 * 
 * This class interfaces with the KaHyPar library to partition
 * hypergraphs generated from VPR netlists.
 */
template<typename HypergraphType>
class PartitioningEngine {
  public:
    using BlockIdType = typename HypergraphType::BlockIdType;
    
    /**
     * @brief Construct partitioning engine
     * @param hypergraph The hypergraph to partition
     * @param num_partitions Number of partitions (k)
     * @param imbalance_rate Maximum allowed imbalance (e.g., 0.03 for 3%)
     * @param config_file KaHyPar configuration file (optional)
     */
    PartitioningEngine(const HypergraphType& hypergraph,
                      int num_partitions,
                      float imbalance_rate,
                      float cost_alpha = 0.0f,
                      float cost_beta = 0.0f,
                      const std::string& config_file = "");
    
    /**
     * @brief Destructor - cleans up KaHyPar context
     */
    ~PartitioningEngine();
    
    /**
     * @brief Run the partitioning algorithm
     * @return true if successful
     */
    bool partition(vtr::vector<AtomBlockId, float>& atom_criticality = nullptr,vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup = nullptr);

    /**
     * @brief Run partitioning using Triton script
     * @param triton_script_path Path to Triton script
     * @param seed Random seed for reproducibility
     * @return true if successful
     */
    bool partition_with_triton(const std::string& triton_script_path, 
                              int seed, vtr::vector<AtomBlockId, float>& atom_criticality,vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup);
    
    /**
     * @brief Get partition ID for a vertex
     * @param vertex_id Vertex ID in hypergraph
     * @return Partition ID (0 to k-1)
     */
    int get_partition(int vertex_id) const;
    
    /**
     * @brief Get partition ID for a block
     * @param block_id Block ID from netlist
     * @return Partition ID (0 to k-1), or -1 if not found
     */
    int get_block_partition(BlockIdType block_id) const;
    
    /**
     * @brief Get the objective value (cut size)
     * @return The hyperedge cut of the partition
     */
    int get_objective() const { return objective_; }
    
    /**
     * @brief Print partitioning statistics
     */
    void print_stats() const;
    
    /**
     * @brief Save partition assignment to file
     * @param filename Output filename
     */
    void save_partition_file(const std::string& filename) const;
    
    /**
     * @brief Save detailed partitioning report
     * @param filename Output filename
     */
    void save_partition_report(const std::string& filename) const;
    
  private:
    const HypergraphType& hypergraph_;
    const int num_partitions_;
    const float imbalance_rate_;
    std::string config_file_;
    const float cost_alpha_;
    const float cost_beta_;

    // KaHyPar context
    int* context_;
    
    // Partitioning results
    std::vector<int> partition_;
    int objective_;
    bool is_partitioned_;

    // For molecule-based partitioning
    bool is_molecule_based_;
    std::map<BlockIdType, int> atom_to_molecule_; // Maps atom ID to molecule vertex ID

    
    bool run_triton_script(const std::string& hypergraph_file, 
                          const std::string& triton_script_path, 
                          int seed);
    bool parse_triton_partition_output(const std::string& partition_file);
    int calculate_cut_objective() const;
};



// Type aliases
using AtomPartitioningEngine = PartitioningEngine<AtomHypergraph>;
using ClusterPartitioningEngine = PartitioningEngine<ClusterHypergraph>;

#endif // VPR_PARTITIONING_ENGINE_H