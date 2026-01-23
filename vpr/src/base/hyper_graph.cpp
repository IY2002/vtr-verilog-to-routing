#include "hyper_graph.h"
#include "atom_netlist.h"
#include "clustered_netlist.h"
#include "vtr_log.h"
#include "vtr_assert.h"
#include <fstream>
#include <set>
#include <map>

// Constructor
template<typename NetlistType>
VprHypergraph<NetlistType>::VprHypergraph(const NetlistType& netlist, size_t min_net_size)
    : netlist_(netlist)
    , min_net_size_(min_net_size)
    , num_vertices_(0)
    , num_hyperedges_(0) {
    
    VTR_ASSERT(min_net_size >= 2);
    build_hypergraph();
}


// Build hypergraph from netlist
template<typename NetlistType>
void VprHypergraph<NetlistType>::build_hypergraph() {
    // Clear any existing data
    block_to_vertex_.clear();
    vertex_to_block_.clear();
    hyperedge_indices_.clear();
    hyperedges_.clear();
    hyperedge_to_net_.clear();
    
    // Step 1: Create vertex for each block
    num_vertices_ = 0;
    for (auto blk_id : netlist_.blocks()) {
        block_to_vertex_[blk_id] = num_vertices_;
        vertex_to_block_[num_vertices_] = blk_id;
        num_vertices_++;
    }
    
    // Step 2: Create hyperedges from nets
    hyperedge_indices_.push_back(0);  // First hyperedge starts at index 0
    num_hyperedges_ = 0;
    
    for (auto net_id : netlist_.nets()) {
        // Collect unique blocks connected to this net
        std::set<unsigned int> connected_vertices;  // Changed from int to unsigned int
        
        for (auto pin_id : netlist_.net_pins(net_id)) {
            auto blk_id = netlist_.pin_block(pin_id);
            auto it = block_to_vertex_.find(blk_id);
            if (it != block_to_vertex_.end()) {
                connected_vertices.insert(static_cast<unsigned int>(it->second));  // Cast to unsigned int
            }
        }
        
        // Only create hyperedge if it meets minimum size requirement
        if (connected_vertices.size() >= min_net_size_) {
            // Add vertices to hyperedge
            for (unsigned int vertex_id : connected_vertices) {  // Changed from int to unsigned int
                hyperedges_.push_back(vertex_id);
            }
            
            // Record where next hyperedge starts
            hyperedge_indices_.push_back(hyperedges_.size());
            
            // Track net ID for this hyperedge
            hyperedge_to_net_.push_back(net_id);
            
            num_hyperedges_++;
        }
    }
    
    VTR_LOG("Built hypergraph: %d vertices, %d hyperedges\n", 
            num_vertices_, num_hyperedges_);

    // save intermediate hypergraph as hMETIS format
    write_hmetis_format("vpr_hypergraph.hgr", std::vector<int>(), std::vector<int>());
}



// // Get edge weights (number of vertices in each hyperedge)
// template<typename NetlistType>
// std::vector<int> VprHypergraph<NetlistType>::get_edge_weights(float alpha) const {
//     std::vector<int> edge_weights;
//     edge_weights.reserve(num_hyperedges_);
    
//     for (int i = 0; i < num_hyperedges_; i++) {
//         int edge_size = hyperedge_indices_[i + 1] - hyperedge_indices_[i];
//         edge_weights.push_back(edge_size);
//     }

//     // Apply alpha scaling
//     for (auto& weight : edge_weights) {
//         if (alpha != 0.0)
//             weight = static_cast<int>(std::pow(weight, alpha) * 1000);  // Scale by alpha and convert to int (1000 for precision)
//         else
//             weight = 1;  // Uniform weight if alpha is 0
//     }

//     return edge_weights;
// }

// Get edge weights based on atom-level timing info
template<typename NetlistType>
std::vector<int> VprHypergraph<NetlistType>::get_edge_weights(
        float alpha,
        vtr::vector<AtomBlockId, float>& atom_criticality,
        vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup) const {

    std::vector<int> edge_weights(num_hyperedges_, 0);

    constexpr double kScale = 1e6;
    constexpr int kMin = 1;

    if constexpr (std::is_same_v<NetlistType, ClusteredNetlist>) {
        // --- Only compiled for ClusteredNetlist ---
        for (int e = 0; e < num_hyperedges_; ++e) {
            NetIdType net = hyperedge_to_net_[e];
            float max_crit = 0.f;

            // For each pin/block connected to this net
            for (auto pin_id : netlist_.net_pins(net)) {
                ClusterBlockId clb = netlist_.pin_block(pin_id);

                if (clb.is_valid() && size_t(clb) < atoms_lookup.size()) {
                    // Find the most critical atom inside this cluster block
                    for (AtomBlockId atom : atoms_lookup[clb]) {
                        if (atom.is_valid() && size_t(atom) < atom_criticality.size()) {
                            max_crit = std::max(max_crit, atom_criticality[atom]);
                        }
                    }
                }
            }

            // Exponentiate by alpha (same logic as vertex weights use beta)
            float pow_exp = (alpha != 0.0f) ? alpha : 1.0f;
            double val = std::pow(max_crit, pow_exp);
            int weight = std::max(kMin, static_cast<int>(std::round(val * kScale)));

            edge_weights[e] = weight;
        }

    } else {
        // --- For AtomNetlist or other cases ---
        for (int e = 0; e < num_hyperedges_; ++e) {
            // fallback: weight proportional to net size (fanout)
            int edge_size = hyperedge_indices_[e + 1] - hyperedge_indices_[e];
            int weight = static_cast<int>(std::pow(edge_size, alpha) * 1000);
            edge_weights[e] = std::max(kMin, weight);
        }
    }

    // save edge weights for debugging
    std::ofstream ew_out("vpr_edge_weights.txt");
    for (int e = 0; e < num_hyperedges_; ++e) {
        ew_out << edge_weights[e] << "\n";
    }
    ew_out.close();

    return edge_weights;
}



// Get vertex weights (number of hyperedges each vertex is part of)
template<typename NetlistType>
std::vector<int> VprHypergraph<NetlistType>::get_vertices_weights(
        float beta,
        vtr::vector<AtomBlockId, float>& atom_criticality, vtr::vector<ClusterBlockId, std::unordered_set<AtomBlockId>>& atoms_lookup, bool timing) const {

    std::vector<int> vertex_weights(num_vertices_, 0);

    constexpr double kScale = 1e6;
    constexpr int kMin = 1;
    std::fill(vertex_weights.begin(), vertex_weights.end(), kMin);
    // if constexpr (std::is_same_v<NetlistType, ClusteredNetlist>) {
    //     // --- Only compiled for ClusteredNetlist ---
    //     if (!timing) {
    //         std::fill(vertex_weights.begin(), vertex_weights.end(), kMin);
    //     }
    //     else {
    //         for (int v = 0; v < num_vertices_; ++v) {
    //             ClusterBlockId clb = vertex_to_block(v);

    //             float max_crit = 0.f;
    //             if (clb.is_valid() && size_t(clb) < atoms_lookup.size()) {
    //                 for (AtomBlockId atom : atoms_lookup[clb]) {
    //                     if (atom.is_valid() && size_t(atom) < atom_criticality.size())
    //                         max_crit = std::max(max_crit, atom_criticality[atom]);
    //                 }
    //             }
    //             float pow_exp = 1.0f;
    //             if (beta != 0.0f) {
    //                 pow_exp = beta;
    //             }
    //             double val = std::pow(max_crit, pow_exp);
    //             int weight = std::max(kMin, static_cast<int>(std::round(val * kScale)));
    //             vertex_weights[v] = weight;
    //         }
        
    //     }

    // } else {
    //     // --- For AtomNetlist or anything else ---
    //     std::fill(vertex_weights.begin(), vertex_weights.end(), kMin);
    // }

    // save vertex weights for debugging
    std::ofstream vw_out("vpr_vertex_weights.txt");
    for (int v = 0; v < num_vertices_; ++v) {
        vw_out << vertex_weights[v] << "\n";
    }
    vw_out.close();

    return vertex_weights;
}

// Convert vertex ID to block ID
template<typename NetlistType>
typename VprHypergraph<NetlistType>::BlockIdType 
VprHypergraph<NetlistType>::vertex_to_block(int vertex_id) const {
    VTR_ASSERT(vertex_id >= 0 && vertex_id < num_vertices_);
    auto it = vertex_to_block_.find(vertex_id);
    VTR_ASSERT(it != vertex_to_block_.end());
    return it->second;
}

// Convert block ID to vertex ID
template<typename NetlistType>
int VprHypergraph<NetlistType>::block_to_vertex(BlockIdType block_id) const {
    auto it = block_to_vertex_.find(block_id);
    if (it != block_to_vertex_.end()) {
        return it->second;
    }
    return -1;  // Block not found
}

// Get net ID for hyperedge
template<typename NetlistType>
typename VprHypergraph<NetlistType>::NetIdType 
VprHypergraph<NetlistType>::hyperedge_to_net(int hyperedge_id) const {
    VTR_ASSERT(hyperedge_id >= 0 && hyperedge_id < num_hyperedges_);
    return hyperedge_to_net_[hyperedge_id];
}

// Print statistics
template<typename NetlistType>
void VprHypergraph<NetlistType>::print_stats() const {
    VTR_LOG("\n");
    VTR_LOG("Hypergraph Statistics:\n");
    VTR_LOG("---------------------\n");
    VTR_LOG("  Vertices (blocks): %d\n", num_vertices_);
    VTR_LOG("  Hyperedges (nets): %d\n", num_hyperedges_);
    VTR_LOG("  Total pins: %zu\n", hyperedges_.size());
    
    if (num_hyperedges_ > 0) {
        double avg_size = static_cast<double>(hyperedges_.size()) / num_hyperedges_;
        VTR_LOG("  Average hyperedge size: %.2f\n", avg_size);
        
        // Calculate hyperedge size distribution
        std::map<int, int> size_distribution;
        int max_size = 0;
        int min_size = INT_MAX;
        
        for (int i = 0; i < num_hyperedges_; i++) {
            int size = hyperedge_indices_[i+1] - hyperedge_indices_[i];
            size_distribution[size]++;
            max_size = std::max(max_size, size);
            min_size = std::min(min_size, size);
        }
        
        VTR_LOG("  Min hyperedge size: %d\n", min_size);
        VTR_LOG("  Max hyperedge size: %d\n", max_size);
        
        // Show distribution for small sizes
        VTR_LOG("\n  Hyperedge size distribution:\n");
        for (const auto& [size, count] : size_distribution) {
            VTR_LOG("    Size %4d: %6d nets (%.1f%%)\n", 
                    size, count, 100.0 * count / num_hyperedges_);
            
        }
    }
    VTR_LOG("\n");
}

template<typename NetlistType>
void VprHypergraph<NetlistType>::write_hmetis_format(
    const std::string& filename,
    const std::vector<int>& vertex_weights,
    const std::vector<int>& edge_weights
) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        VTR_LOG_ERROR("Failed to open file '%s' for writing\n", filename.c_str());
        return;
    }

    // Determine format flag (fmt)
    int fmt = 0;
    if (!edge_weights.empty()) fmt += 1;   // bit 0
    if (!vertex_weights.empty()) fmt += 10; // bit 1

    // Header: num_hyperedges num_vertices fmt
    out << num_hyperedges_ << " " << num_vertices_;
    if (fmt > 0)
        out << " " << fmt;
    out << "\n";

    // --- Write hyperedges (1-based vertex indices) ---
    for (int i = 0; i < num_hyperedges_; ++i) {
        // Optional edge weight
        if (!edge_weights.empty()) {
            out << edge_weights[i] << " ";
        }

        size_t start = hyperedge_indices_[i];
        size_t end = hyperedge_indices_[i + 1];

        for (size_t j = start; j < end; ++j) {
            out << (hyperedges_[j] + 1);
            if (j < end - 1)
                out << " ";
        }
        out << "\n";
    }

    // --- Write vertex weights ---
    if (!vertex_weights.empty()) {
        for (size_t i = 0; i < vertex_weights.size(); ++i) {
            out << vertex_weights[i] << "\n";
        }
    }

    out.close();
    VTR_LOG("Wrote hypergraph to '%s' in hMETIS format\n", filename.c_str());
}

// Non-member function to build molecule-based hypergraph
std::map<t_pack_molecule*, int> build_molecule_hypergraph(VprHypergraph<AtomNetlist>& hypergraph, 
                              const AtomNetlist& netlist,
                              const Prepacker& prepacker) {
    // Clear any existing data
    hypergraph.block_to_vertex_.clear();
    hypergraph.vertex_to_block_.clear();
    hypergraph.hyperedge_indices_.clear();
    hypergraph.hyperedges_.clear();
    hypergraph.hyperedge_to_net_.clear();
    
    // Step 1: Create vertex for each molecule
    std::map<t_pack_molecule*, int> molecule_to_vertex;
    hypergraph.num_vertices_ = 0;
    
    auto molecules = prepacker.get_molecules_vector();
    for (size_t mol_idx = 0; mol_idx < molecules.size(); mol_idx++) {
        t_pack_molecule* mol_ptr = molecules[mol_idx];
        molecule_to_vertex[mol_ptr] = hypergraph.num_vertices_;
        
        // Store the first atom of the molecule as representative
        if (!mol_ptr->atom_block_ids.empty()) {
            auto repr_atom = mol_ptr->atom_block_ids[0];
            hypergraph.vertex_to_block_[hypergraph.num_vertices_] = repr_atom;
            hypergraph.block_to_vertex_[repr_atom] = hypergraph.num_vertices_;
        }
        hypergraph.num_vertices_++;
    }
    
    // Step 2: Create hyperedges from nets
    hypergraph.hyperedge_indices_.push_back(0);
    hypergraph.num_hyperedges_ = 0;
    
    for (auto net_id : netlist.nets()) {
        // Collect unique molecules connected to this net
        std::set<unsigned int> connected_molecules;
        
        for (auto pin_id : netlist.net_pins(net_id)) {
            auto blk_id = netlist.pin_block(pin_id);
            t_pack_molecule* molecule_ptr = prepacker.get_atom_molecule(blk_id);
            
            if (molecule_ptr != nullptr) {
                auto it = molecule_to_vertex.find(molecule_ptr);
                if (it != molecule_to_vertex.end()) {
                    connected_molecules.insert(static_cast<unsigned int>(it->second));
                }
            }
        }
        
        // Only create hyperedge if it meets minimum size requirement
        if (connected_molecules.size() >= hypergraph.min_net_size_) {
            // Add vertices to hyperedge
            for (unsigned int vertex_id : connected_molecules) {
                hypergraph.hyperedges_.push_back(vertex_id);
            }
            
            // Record where next hyperedge starts
            hypergraph.hyperedge_indices_.push_back(hypergraph.hyperedges_.size());
            
            // Track net ID for this hyperedge
            hypergraph.hyperedge_to_net_.push_back(net_id);
            
            hypergraph.num_hyperedges_++;
        }
    }
    
    VTR_LOG("Built molecule hypergraph: %d vertices (molecules), %d hyperedges\n", 
            hypergraph.num_vertices_, hypergraph.num_hyperedges_);

    return molecule_to_vertex;  // Return mapping of molecules to vertex IDs
}

// Explicit instantiations
template class VprHypergraph<AtomNetlist>;
template class VprHypergraph<ClusteredNetlist>;

