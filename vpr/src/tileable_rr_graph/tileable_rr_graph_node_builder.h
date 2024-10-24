#ifndef TILEABLE_RR_GRAPH_NODE_BUILDER_H
#define TILEABLE_RR_GRAPH_NODE_BUILDER_H

/********************************************************************
 * Include header files that are required by function declaration
 *******************************************************************/
/* Headers from vtrutil library */
#include "vtr_geometry.h"

/* Headers from readarch library */
#include "physical_types.h"

/* Headers from vpr library */
#include "device_grid.h"
#include "device_grid_annotation.h"
#include "rr_node_types.h"
#include "rr_graph_type.h"
#include "rr_graph_view.h"
#include "rr_graph_builder.h"

/********************************************************************
 * Function declaration
 *******************************************************************/

void alloc_tileable_rr_graph_nodes(RRGraphBuilder& rr_graph_builder,
                                   vtr::vector<RRNodeId, RRSwitchId>& driver_switches,
                                   const DeviceGrid& grids,
                                   const size_t& layer,
                                   const vtr::Point<size_t>& chan_width,
                                   const std::vector<t_segment_inf>& segment_inf_x,
                                   const std::vector<t_segment_inf>& segment_inf_y,
                                   const DeviceGridAnnotation& device_grid_annotation,
                                   const bool& shrink_boundary,
                                   const bool& perimeter_cb,
                                   const bool& through_channel);

void create_tileable_rr_graph_nodes(const RRGraphView& rr_graph,
                                    RRGraphBuilder& rr_graph_builder,
                                    vtr::vector<RRNodeId, RRSwitchId>& rr_node_driver_switches,
                                    std::map<RRNodeId, std::vector<size_t>>& rr_node_track_ids,
                                    std::vector<t_rr_rc_data>& rr_rc_data,
                                    const DeviceGrid& grids,
                                    const size_t& layer,
                                    const vtr::Point<size_t>& chan_width,
                                    const std::vector<t_segment_inf>& segment_inf_x,
                                    const std::vector<t_segment_inf>& segment_inf_y,
                                    const t_unified_to_parallel_seg_index& segment_index_map,
                                    const RRSwitchId& wire_to_ipin_switch,
                                    const RRSwitchId& delayless_switch,
                                    const DeviceGridAnnotation& device_grid_annotation,
                                    const bool& shrink_boundary,
                                    const bool& perimeter_cb,
                                    const bool& through_channel);

#endif
