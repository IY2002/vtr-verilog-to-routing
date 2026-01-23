#ifndef VPR_HYPERGRAPH_H
#define VPR_HYPERGRAPH_H

#include <vector>
#include <unordered_map>
#include <string>
#include "netlist.h"
#include "clustered_netlist.h"
#include "prepack.h"

// Type traits to map netlist types to their ID types
template<typename NetlistType>
struct NetlistTraits;

// Specialization for AtomNetlist
template<>
struct NetlistTraits<AtomNetlist> {
    using block_id_type = AtomBlockId;
    using net_id_type = AtomNetId;
    using pin_id_type = AtomPinId;
};

// Specialization for ClusteredNetlist
template<>
struct NetlistTraits<ClusteredNetlist> {
    using block_id_type = ClusterBlockId;
    using net_id_type = ClusterNetId;
    using pin_id_type = ClusterPinId;
};


/**
 * @brief A hypergraph representation of VPR netlists
 * 
 * This class converts VPR netlists into a hypergraph format where:
 * - Vertices represent blocks
 * - Hyperedges represent nets
 * 
 * The hypergraph data is stored in a compressed format suitable
 * for hypergraph partitioners.
 */
template<typename NetlistType>
class VprHypergraph {
  public:
    using BlockIdType = typename NetlistTraits<NetlistType>::block_id_type;
    using NetIdType = typename NetlistTraits<NetlistType>::net_id_type;
    using PinIdType = typename NetlistTraits<NetlistType>::pin_id_type;
    
    /**
     * @brief Construct hypergraph from a netlist
     * @param netlist The netlist to convert
     * @param min_net_size Minimum net size to include (default 2)
     */
    explicit VprHypergraph(const NetlistType& netlist, size_t min_net_size = 2);

    
    // Accessors
    int num_vertices() const { return num_vertices_; }
    int num_hyperedges() const { return num_hyperedges_; }
    
    /**
     * @brief Get the hyperedge indices array
     * 
     * hyperedge_indices[i] points to the start of hyperedge i in the hyperedges array
     * hyperedge_indices[i+1] - hyperedge_indices[i] gives the size of hyperedge i
     */
    const std::vector<size_t>& hyperedge_indices() const { 
        return hyperedge_indices_; 
    }

    /**
     * @brief Get the hyperedge sizes
     * 
     * Returns a vector containing the size of each hyperedge (number of vertices in each hyperedge)
     * The size of the vector is num_hyperedges_
     */
    std::vector<int> get_edge_weights(float alpha, vtr::vector<AtomBlockId, float>& atom_criticality = nullptr, vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup = nullptr) const;

    /**
     * @brief Get the vertex weights
     * 
     * Returns a vector containing the number of hyperedges each vertex is part of
     * The size of the vector is num_vertices_
     */
    std::vector<int> get_vertices_weights(float beta, vtr::vector<AtomBlockId, float>& atom_criticality = nullptr, vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup = nullptr, bool timing = false) const;
    
    /**
     * @brief Get the hyperedges array
     * 
     * Flattened array containing all vertex IDs in all hyperedges
     * Use hyperedge_indices to determine which vertices belong to which hyperedge
     */
    const std::vector<unsigned int>& hyperedges() const { 
        return hyperedges_; 
    }
    
    /**
     * @brief Convert vertex ID to block ID
     * @param vertex_id The vertex ID in hypergraph (0 to num_vertices-1)
     * @return The corresponding block ID in netlist
     */
    BlockIdType vertex_to_block(int vertex_id) const;
    
    /**
     * @brief Convert block ID to vertex ID
     * @param block_id The block ID in netlist
     * @return The corresponding vertex ID in hypergraph (-1 if not found)
     */
    int block_to_vertex(BlockIdType block_id) const;
    
    /**
     * @brief Get the net ID for a hyperedge
     * @param hyperedge_id The hyperedge index (0 to num_hyperedges-1)
     * @return The corresponding net ID in netlist
     */
    NetIdType hyperedge_to_net(int hyperedge_id) const;
    
    /**
     * @brief Print hypergraph statistics to log
     */
    void print_stats() const;
    
    /**
     * @brief Write hypergraph to file in hMETIS format
     * @param filename Output filename
     */
    void write_hmetis_format(const std::string& filename, const std::vector<int>& vertex_weights, const std::vector<int>& edge_weights) const;
    
    /**
     * @brief Build the hypergraph from the netlist
     */
    void build_hypergraph();

    // Add this friend declaration to allow access to private members
    friend std::map<t_pack_molecule*, int> build_molecule_hypergraph(VprHypergraph<AtomNetlist>& hypergraph, 
                                         const AtomNetlist& netlist,
                                         const Prepacker& prepacker);

    
  private:
    const NetlistType& netlist_;
    const size_t min_net_size_;
    
    // Hypergraph data
    int num_vertices_;
    int num_hyperedges_;
    std::vector<size_t> hyperedge_indices_;
    std::vector<unsigned int> hyperedges_;
    
    // Mappings between hypergraph and netlist
    std::unordered_map<BlockIdType, int> block_to_vertex_;
    std::unordered_map<int, BlockIdType> vertex_to_block_;
    std::vector<NetIdType> hyperedge_to_net_;
};
 



// Explicit instantiations
extern template class VprHypergraph<AtomNetlist>;
extern template class VprHypergraph<ClusteredNetlist>;

// Type aliases for convenience
using AtomHypergraph = VprHypergraph<AtomNetlist>;
using ClusterHypergraph = VprHypergraph<ClusteredNetlist>;


#endif // VPR_HYPERGRAPH_H