#ifndef OPENFPGA_RR_GRAPH_UTILS_H
#define OPENFPGA_RR_GRAPH_UTILS_H

/********************************************************************
 * Include header files that are required by function declaration
 *******************************************************************/
/* Headers from vtrutil library */
#include "vtr_geometry.h"

/* Headers from vpr library */
#include "rr_graph_obj.h"
#include "rr_graph_view.h"

/********************************************************************
 * Function declaration
 *******************************************************************/

vtr::Point<size_t> get_track_rr_node_start_coordinate(const RRGraphView& rr_graph,
                                                      const RRNodeId& track_rr_node);

vtr::Point<size_t> get_track_rr_node_end_coordinate(const RRGraphView& rr_graph,
                                                    const RRNodeId& track_rr_node);

std::vector<RRSwitchId> get_rr_graph_driver_switches(const RRGraphView& rr_graph,
                                                     const RRNodeId& node);

std::vector<RRNodeId> get_rr_graph_driver_nodes(const RRGraphView& rr_graph,
                                                const RRNodeId& node);

std::vector<RRNodeId> get_rr_graph_configurable_driver_nodes(const RRGraphView& rr_graph,
                                                             const RRNodeId& node);

std::vector<RRNodeId> get_rr_graph_non_configurable_driver_nodes(const RRGraphView& rr_graph,
                                                                 const RRNodeId& node);

bool is_opin_direct_connected_ipin(const RRGraphView& rr_graph,
                                   const RRNodeId& node);

bool is_ipin_direct_connected_opin(const RRGraphView& rr_graph,
                                   const RRNodeId& node);

e_side get_rr_graph_single_node_side(const RRGraphView& rr_graph,
                                     const RRNodeId& node);

#endif
