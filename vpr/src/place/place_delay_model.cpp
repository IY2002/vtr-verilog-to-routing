/**
 * @file place_delay_model.cpp
 * @brief This file implements all the class methods and individual
 *        routines related to the placer delay model.
 */

#include <queue>
#include "place_delay_model.h"
#include "globals.h"
#include "router_lookahead_map.h"
#include "rr_graph2.h"

#include "timing_place_lookup.h"
#include "placer_state.h"

#include "vtr_log.h"
#include "vtr_math.h"
#include "vpr_error.h"

#ifdef VTR_ENABLE_CAPNPROTO
#    include "capnp/serialize.h"
#    include "place_delay_model.capnp.h"
#    include "ndmatrix_serdes.h"
#    include "mmap_file.h"
#    include "serdes_utils.h"
#endif /* VTR_ENABLE_CAPNPROTO */

///@brief DeltaDelayModel methods.
float DeltaDelayModel::delay(const t_physical_tile_loc& from_loc, int /*from_pin*/, const t_physical_tile_loc& to_loc, int /*to_pin*/,  e_side from_side, e_side to_side) const {
    int delta_x = std::abs(from_loc.x - to_loc.x);
    int delta_y = std::abs(from_loc.y - to_loc.y);

    return delays_[from_loc.layer_num][to_loc.layer_num][delta_x][delta_y];
}

void DeltaDelayModel::dump_echo(std::string filepath) const {
    FILE* f = vtr::fopen(filepath.c_str(), "w");
    fprintf(f, "         ");
    for (size_t from_layer_num = 0; from_layer_num < delays_.dim_size(0); ++from_layer_num) {
        for (size_t to_layer_num = 0; to_layer_num < delays_.dim_size(1); ++to_layer_num) {
            fprintf(f, " %9zu", from_layer_num);
            fprintf(f, "\n");
            for (size_t dx = 0; dx < delays_.dim_size(2); ++dx) {
                fprintf(f, " %9zu", dx);
            }
            fprintf(f, "\n");
            for (size_t dy = 0; dy < delays_.dim_size(3); ++dy) {
                fprintf(f, "%9zu", dy);
                for (size_t dx = 0; dx < delays_.dim_size(2); ++dx) {
                    fprintf(f, " %9.2e", delays_[from_layer_num][to_layer_num][dx][dy]);
                }
                fprintf(f, "\n");
            }
        }
    }
    vtr::fclose(f);
}

const DeltaDelayModel* OverrideDelayModel::base_delay_model() const {
    return base_delay_model_.get();
}

///@brief OverrideDelayModel methods.
float OverrideDelayModel::delay(const t_physical_tile_loc& from_loc, int from_pin, const t_physical_tile_loc& to_loc, int to_pin,  e_side from_side, e_side to_side) const {
    //First check to if there is an override delay value
    auto& device_ctx = g_vpr_ctx.device();
    auto& grid = device_ctx.grid;

    t_physical_tile_type_ptr from_type_ptr = grid.get_physical_type(from_loc);
    t_physical_tile_type_ptr to_type_ptr = grid.get_physical_type(to_loc);

    t_override override_key;
    override_key.from_type = from_type_ptr->index;
    override_key.from_class = from_type_ptr->pin_class[from_pin];
    override_key.to_type = to_type_ptr->index;
    override_key.to_class = to_type_ptr->pin_class[to_pin];

    //Delay overrides may be different for +/- delta so do not use
    //an absolute delta for the look-up
    override_key.delta_x = to_loc.x - from_loc.x;
    override_key.delta_y = to_loc.y - from_loc.y;

    float delay_val = std::numeric_limits<float>::quiet_NaN();
    auto override_iter = delay_overrides_.find(override_key);
    if (override_iter != delay_overrides_.end()) {
        //Found an override
        delay_val = override_iter->second;
    } else {
        //Fall back to the base delay model if no override was found
        delay_val = base_delay_model_->delay(from_loc, from_pin, to_loc, to_pin, NUM_2D_SIDES, NUM_2D_SIDES);
    }

    return delay_val;
}

void OverrideDelayModel::set_delay_override(int from_type, int from_class, int to_type, int to_class, int delta_x, int delta_y, float delay_val) {
    t_override override_key;
    override_key.from_type = from_type;
    override_key.from_class = from_class;
    override_key.to_type = to_type;
    override_key.to_class = to_class;
    override_key.delta_x = delta_x;
    override_key.delta_y = delta_y;

    auto res = delay_overrides_.insert(std::make_pair(override_key, delay_val));
    if (!res.second) {                 //Key already exists
        res.first->second = delay_val; //Overwrite existing delay
    }
}

void OverrideDelayModel::dump_echo(std::string filepath) const {
    base_delay_model_->dump_echo(filepath);

    FILE* f = vtr::fopen(filepath.c_str(), "a");

    fprintf(f, "\n");
    fprintf(f, "# Delay Overrides\n");
    auto& device_ctx = g_vpr_ctx.device();
    for (auto kv : delay_overrides_) {
        auto override_key = kv.first;
        float delay_val = kv.second;
        fprintf(f, "from_type: %s to_type: %s from_pin_class: %d to_pin_class: %d delta_x: %d delta_y: %d -> delay: %g\n",
                device_ctx.physical_tile_types[override_key.from_type].name.c_str(),
                device_ctx.physical_tile_types[override_key.to_type].name.c_str(),
                override_key.from_class,
                override_key.to_class,
                override_key.delta_x,
                override_key.delta_y,
                delay_val);
    }

    vtr::fclose(f);
}

float OverrideDelayModel::get_delay_override(int from_type, int from_class, int to_type, int to_class, int delta_x, int delta_y) const {
    t_override key;
    key.from_type = from_type;
    key.from_class = from_class;
    key.to_type = to_type;
    key.to_class = to_class;
    key.delta_x = delta_x;
    key.delta_y = delta_y;

    auto iter = delay_overrides_.find(key);
    if (iter == delay_overrides_.end()) {
        VPR_THROW(VPR_ERROR_PLACE, "Key not found.");
    }
    return iter->second;
}

void OverrideDelayModel::set_base_delay_model(std::unique_ptr<DeltaDelayModel> base_delay_model_obj) {
    base_delay_model_ = std::move(base_delay_model_obj);
}

float SimpleDelayModel::delay(const t_physical_tile_loc& from_loc, int /*from_pin*/, const t_physical_tile_loc& to_loc, int /*to_pin*/,  e_side from_side, e_side to_side) const {
    int delta_x = std::abs(from_loc.x - to_loc.x);
    int delta_y = std::abs(from_loc.y - to_loc.y);

    int from_tile_idx = g_vpr_ctx.device().grid.get_physical_type(from_loc)->index;
    return delays_[from_tile_idx][from_loc.layer_num][to_loc.layer_num][delta_x][delta_y];
}

// Helper functions
e_side get_opposite_side(e_side side) {
    switch(side) {
        case TOP: return BOTTOM;
        case BOTTOM: return TOP;
        case LEFT: return RIGHT;
        case RIGHT: return LEFT;
        default: return NUM_2D_SIDES;
    }
}

e_side get_perpendicular_side(e_side side) {
    switch(side) {
        case TOP: return LEFT;
        case BOTTOM: return RIGHT;
        case LEFT: return TOP;
        case RIGHT: return BOTTOM;
        default: return NUM_2D_SIDES;
    }
}

// float SimpleDelayModel::delay(const t_physical_tile_loc& from_loc, int from_pin, const t_physical_tile_loc& to_loc, int to_pin,  e_side from_side, e_side to_side) const {
//     // Get original deltas (signed)
//     int original_delta_x = to_loc.x - from_loc.x;
//     int original_delta_y = to_loc.y - from_loc.y;
    
//     // Transform coordinates based on source side (matching Python logic)
//     int delta_x, delta_y;
//     delta_x = original_delta_x;
//     delta_y = original_delta_y;
    
//     if (from_side == TOP) {
//         // Scenario 1
//         if (delta_x == 0 && delta_y == 1){
//             if (to_side == TOP) // + 1 penalty
//                 delta_y += 1;

//             if (to_side == BOTTOM) // - 1 savings
//                 delta_y -= 1;
//         } 
//         // Scenario 2
//         else if (delta_x == 0 && delta_y > 1){
//             if (to_side == TOP) // + 1 penalty
//                 delta_y += 1;
//         }
//         // Scenario 3
//         else if (delta_x == 0 && delta_y < 0){
//             if (to_side == BOTTOM) // + 2 penalty
//                 delta_y -= 2;
//             else // + 1 penalty
//                 delta_y -= 1;
//         }
//         // Scenario 4
//         else if (delta_y > 0 && delta_x != 0){
//             if (delta_x > 0){
//                 if (to_side == BOTTOM) // - 1 saving
//                     delta_y -= 1;
//                 else if (to_side == LEFT) // - 1 saving
//                     delta_x -= 1;
//             }
//             else {
//                 if (to_side == BOTTOM) // - 1 saving
//                     delta_y -= 1;
//                 else if (to_side == RIGHT) // - 1 saving
//                     delta_x += 1;
//             }
//         }
//         // Scenario 5
//         else if (delta_y <= 0 && delta_x != 0){
//             if (delta_x > 0){
//                 if (to_side == BOTTOM) // + 1 penalty
//                     delta_y -= 1;
//                 else if (to_side == RIGHT) // + 1 penalty
//                     delta_x += 1;
//             }
//             else {
//                 if (to_side == BOTTOM) // + 1 penalty
//                     delta_y -= 1;

//                 else if (to_side == LEFT) // + 1 penalty
//                     delta_x -= 1;
//             }
//         }
//     } else if (from_side == BOTTOM) {
//         // Scenario 1
//         if (delta_x == 0 && delta_y == -1){
//             if (to_side == BOTTOM) // + 1 penalty
//                 delta_y -= 1;

//             if (to_side == TOP) // - 1 savings
//                 delta_y += 1;
//         } 
//         // Scenario 2
//         else if (delta_x == 0 && delta_y < -1){
//             if (to_side == BOTTOM) // + 1 penalty
//                 delta_y -= 1;
//         }
//         // Scenario 3
//         else if (delta_x == 0 && delta_y > 0){
//             if (to_side == TOP) // + 2 penalty
//                 delta_y += 2;
//             else // + 1 penalty
//                 delta_y += 1;
//         }

//         // Scenario 4
//         else if (delta_y < 0 && delta_x != 0){
//             if (delta_x > 0){
//                 if (to_side == TOP) // - 1 saving
//                     delta_y += 1;
                
//                 else if (to_side == LEFT) // - 1 saving
//                     delta_x -= 1;

//             }
//             else {
//                 if (to_side == TOP) // - 1 saving
//                     delta_y += 1;

//                 else if (to_side == RIGHT) // - 1 saving
//                     delta_x += 1;
//             }
//         }
//         // Scenario 5
//         else if (delta_y >= 0 && delta_x != 0){
//             if (delta_x > 0){
//                 if (to_side == TOP) // + 1 penalty
//                     delta_y += 1;
//                 else if (to_side == RIGHT) // + 1 penalty
//                     delta_x += 1;
//             }
//             else {
//                 if (to_side == TOP) // + 1 penalty
//                     delta_y += 1;

//                 else if (to_side == LEFT) // + 1 penalty
//                     delta_x -= 1;
//             }
//         }
//     } else if (from_side == LEFT) {
//         // Scenario 1
//         if (delta_y == 0 && delta_x == -1){
//             if (to_side == LEFT) // + 1 penalty
//                 delta_x -= 1;

//             if (to_side == RIGHT) // - 1 savings
//                 delta_x += 1;
//         } 
//         // Scenario 2
//         else if (delta_y == 0 && delta_x < -1){
//             if (to_side == LEFT) // + 1 penalty
//                 delta_x -= 1;
//         }
//         // Scenario 3
//         else if (delta_y == 0 && delta_x > 0){
//             if (to_side == RIGHT) // + 2 penalty
//                 delta_x += 2;
//             else // + 1 penalty
//                 delta_x += 1;
//         }

//         // Scenario 4
//         else if (delta_x < 0 && delta_y != 0){
//             if (delta_y > 0){
//                 if (to_side == RIGHT) // - 1 saving
//                     delta_x += 1;
//                 else if (to_side == BOTTOM) // - 1 saving
//                     delta_y -= 1;
//             }
//             else {
//                 if (to_side == RIGHT) // - 1 saving
//                     delta_x += 1;

//                 else if (to_side == TOP) // - 1 saving
//                     delta_y += 1;
//             }
//         }
//         // Scenario 5
//         else if (delta_x >= 0 && delta_y != 0){
//             if (delta_y > 0){
//                 if (to_side == RIGHT) // + 1 penalty
//                     delta_x += 1;
                
//                 else if (to_side == TOP) // + 1 penalty
//                     delta_y += 1;
//             }
//             else {
//                 if (to_side == RIGHT) // + 1 penalty
//                     delta_x += 1;

//                 else if (to_side == BOTTOM) // + 1 penalty
//                     delta_y -= 1;
//             }
//         }
//     } else {
//         VTR_ASSERT(from_side == RIGHT);
//         // Scenario 1
//         if (delta_y == 0 && delta_x == 1){
//             if (to_side == RIGHT) // + 1 penalty
//                 delta_x += 1;

//             if (to_side == LEFT) // - 1 savings
//                 delta_x -= 1;
//         } 
//         // Scenario 2
//         else if (delta_y == 0 && delta_x > 1){
//             if (to_side == RIGHT) // + 1 penalty
//                 delta_x += 1;
//         }
//         // Scenario 3
//         else if (delta_y == 0 && delta_x < 0){
//             if (to_side == LEFT) // + 2 penalty
//                 delta_x -= 2;
//             else // + 1 penalty
//                 delta_x -= 1;
//         }

//         // Scenario 4
//         else if (delta_x > 0 && delta_y != 0){
//             if (delta_y > 0){
//                 if (to_side == LEFT) // - 1 saving
//                     delta_x -= 1;

//                 else if (to_side == BOTTOM) // - 1 saving
//                     delta_y -= 1;
//             }
//             else {
//                 if (to_side == LEFT) // - 1 saving
//                     delta_x -= 1;
                
//                 else if (to_side == TOP) // - 1 saving
//                     delta_y += 1;
//             }
//         }
//         // Scenario 5
//         else if (delta_x <= 0 && delta_y != 0){
//             if (delta_y > 0){
//                 if (to_side == LEFT) // + 1 penalty
//                     delta_x -= 1;
                
//                 else if (to_side == TOP) // + 1 penalty
//                     delta_y += 1;
                
//             }
//             else {
//                 if (to_side == LEFT) // + 1 penalty
//                     delta_x -= 1;

//                 else if (to_side == BOTTOM) // + 1 penalty
//                     delta_y -= 1;
//             }
//         }
//     }
    
    

//     // // Initialize return values
//     // int ret_x = delta_x;
//     // int ret_y = delta_y;
    
//     // // Implement the scenario logic
//     // if (delta_x == 0 && delta_y == 1) { // Scenario 1: Directly adjacent side 
//     //     if (to_side == from_side) {
//     //         ret_y += 1;  // Penalty
//     //     }
//     //     if (to_side == get_opposite_side(from_side)) {
//     //         ret_y -= 1;  // Savings
//     //     }
//     // } else if (delta_x == 0 && delta_y > 1) { // Scenario 2: Adjacent side farther in source direction
//     //     if (to_side == from_side) {
//     //         ret_y += 1;  // Penalty
//     //     }
//     // } else if (delta_x > 0 && delta_y > 0) { // Scenario 3: Positive quadrant.  
//     //     if (to_side == get_perpendicular_side(from_side)) {
//     //         ret_x -= 1;  // Savings
//     //     }
//     //     if (to_side == get_opposite_side(from_side)) {
//     //         ret_y -= 1;  // Savings
//     //     }
//     // } else if (delta_x > 0 && delta_y <= 0) { // Scenario 4: Mixed quadrant
//     //     if (to_side == get_opposite_side(get_perpendicular_side(from_side))) {
//     //         ret_x += 1;  // Penalty
//     //     }
//     //     if (to_side == get_opposite_side(from_side)) {
//     //         ret_y -= 1;  // Penalty
//     //     }
//     // } else if (delta_x == 0 && delta_y < 0) { // Scenario 5: Opposite direction
//     //     if (to_side == get_opposite_side(from_side)) {
//     //         ret_y -= 2;  // Big penalty
//     //     } else if (to_side == from_side) {
//     //         ret_y -= 1;  // Standard penalty
//     //     } else {
//     //         ret_x -= 1;  // Other sides penalty
//     //     }
//     // } else if (delta_x < 0 && delta_y <= 0) { // Scenario 6: Negative quadrant
//     //     if (to_side == get_perpendicular_side(from_side)) {
//     //         ret_x -= 1;  // Penalty
//     //     }
//     //     if (to_side == get_opposite_side(from_side)) {
//     //         ret_y -= 1;  // Penalty
//     //     }
//     // } else if (delta_x < 0 && delta_y > 0) { // Scenario 7: Mixed negative quadrant
//     //     if (to_side == get_opposite_side(get_perpendicular_side(from_side))) {
//     //         ret_x += 1;  // Savings
//     //     }
//     //     if (to_side == get_opposite_side(from_side)) {
//     //         ret_y -= 1;  // Savings
//     //     }
//     // }
    
//     // Use absolute values for delay lookup (since delay table expects positive indices)
//     int lookup_delta_x = std::abs(delta_x);
//     int lookup_delta_y = std::abs(delta_y);

    
//     int from_tile_idx = g_vpr_ctx.device().grid.get_physical_type(from_loc)->index;
//     return delays_[from_tile_idx][from_loc.layer_num][to_loc.layer_num][lookup_delta_x][lookup_delta_y];
// }

void SimpleDelayModel::write(const std::string& file) const {
    // Use text-based output since you're not interested in Cap'n Proto
    FILE* f = vtr::fopen(file.c_str(), "w");
    
    // Write header information
    fprintf(f, "# SimpleDelayModel Delay Matrix\n");
    fprintf(f, "# Format: [physical_type][from_layer][to_layer][delta_x][delta_y] = delay\n");
    fprintf(f, "# Dimensions: %zu physical_types, %zu layers, %zu max_dx, %zu max_dy\n\n",
            delays_.dim_size(0), delays_.dim_size(1), delays_.dim_size(3), delays_.dim_size(4));
    
    auto& device_ctx = g_vpr_ctx.device();
    
    // Write the delay values
    for (size_t phys_type = 0; phys_type < delays_.dim_size(0); ++phys_type) {
        fprintf(f, "# Physical Type: %s (index: %zu)\n", 
                device_ctx.physical_tile_types[phys_type].name.c_str(), phys_type);
        
        for (size_t from_layer = 0; from_layer < delays_.dim_size(1); ++from_layer) {
            for (size_t to_layer = 0; to_layer < delays_.dim_size(2); ++to_layer) {
                fprintf(f, "## From Layer %zu to Layer %zu\n", from_layer, to_layer);
                
                // Print header row with delta_x values
                fprintf(f, "     ");
                for (size_t dx = 0; dx < delays_.dim_size(3); ++dx) {
                    fprintf(f, " %8zu", dx);
                }
                fprintf(f, "\n");
                
                // Print each row with delta_y and delay values
                for (size_t dy = 0; dy < delays_.dim_size(4); ++dy) {
                    fprintf(f, "%4zu:", dy);
                    for (size_t dx = 0; dx < delays_.dim_size(3); ++dx) {
                        fprintf(f, " %8.3e", delays_[phys_type][from_layer][to_layer][dx][dy]);
                    }
                    fprintf(f, "\n");
                }
                fprintf(f, "\n");
            }
        }
    }
    
    vtr::fclose(f);
}


#ifndef VTR_ENABLE_CAPNPROTO

#    define DISABLE_ERROR                              \
        "is disable because VTR_ENABLE_CAPNPROTO=OFF." \
        "Re-compile with CMake option VTR_ENABLE_CAPNPROTO=ON to enable."

void DeltaDelayModel::read(const std::string& /*file*/) {
    VPR_THROW(VPR_ERROR_PLACE, "DeltaDelayModel::read " DISABLE_ERROR);
}

void DeltaDelayModel::write(const std::string& /*file*/) const {
    VPR_THROW(VPR_ERROR_PLACE, "DeltaDelayModel::write " DISABLE_ERROR);
}

void OverrideDelayModel::read(const std::string& /*file*/) {
    VPR_THROW(VPR_ERROR_PLACE, "OverrideDelayModel::read " DISABLE_ERROR);
}

void OverrideDelayModel::write(const std::string& /*file*/) const {
    VPR_THROW(VPR_ERROR_PLACE, "OverrideDelayModel::write " DISABLE_ERROR);
}

#else /* VTR_ENABLE_CAPNPROTO */

static void ToFloat(float* out, const VprFloatEntry::Reader& in) {
    // Getting a scalar field is always "get<field name>()".
    *out = in.getValue();
}

static void FromFloat(VprFloatEntry::Builder* out, const float& in) {
    // Setting a scalar field is always "set<field name>(value)".
    out->setValue(in);
}

void DeltaDelayModel::read(const std::string& file) {
    // MmapFile object creates an mmap of the specified path, and will munmap
    // when the object leaves scope.
    MmapFile f(file);

    /* Increase reader limit to 1G words to allow for large files. */
    ::capnp::ReaderOptions opts = default_large_capnp_opts();

    // FlatArrayMessageReader is used to read the message from the data array
    // provided by MmapFile.
    ::capnp::FlatArrayMessageReader reader(f.getData(), opts);

    // When reading capnproto files the Reader object to use is named
    // <schema name>::Reader.
    //
    // Initially this object is an empty VprDeltaDelayModel.
    VprDeltaDelayModel::Reader model;

    // The reader.getRoot performs a cast from the generic capnproto to fit
    // with the specified schema.
    //
    // Note that capnproto does not validate that the incoming data matches the
    // schema.  If this property is required, some form of check would be
    // required.
    model = reader.getRoot<VprDeltaDelayModel>();

    // ToNdMatrix is a generic function for converting a Matrix capnproto
    // to a vtr::NdMatrix.
    //
    // The use must supply the matrix dimension (2 in this case), the source
    // capnproto type (VprFloatEntry),
    // target C++ type (flat), and a function to convert from the source capnproto
    // type to the target C++ type (ToFloat).
    //
    // The second argument should be of type Matrix<X>::Reader where X is the
    // capnproto element type.
    ToNdMatrix<4, VprFloatEntry, float>(&delays_, model.getDelays(), ToFloat);
}

void DeltaDelayModel::write(const std::string& file) const {
    // MallocMessageBuilder object is the generate capnproto message builder,
    // using malloc for buffer allocation.
    ::capnp::MallocMessageBuilder builder;

    // initRoot<X> returns a X::Builder object that can be used to set the
    // fields in the message.
    auto model = builder.initRoot<VprDeltaDelayModel>();

    // FromNdMatrix is a generic function for converting a vtr::NdMatrix to a
    // Matrix message.  It is the mirror function of ToNdMatrix described in
    // read above.
    auto delay_values = model.getDelays();
    FromNdMatrix<4, VprFloatEntry, float>(&delay_values, delays_, FromFloat);

    // writeMessageToFile writes message to the specified file.
    writeMessageToFile(file, &builder);
}

void OverrideDelayModel::read(const std::string& file) {
    MmapFile f(file);

    /* Increase reader limit to 1G words to allow for large files. */
    ::capnp::ReaderOptions opts = default_large_capnp_opts();
    ::capnp::FlatArrayMessageReader reader(f.getData(), opts);

    vtr::NdMatrix<float, 4> delays;
    auto model = reader.getRoot<VprOverrideDelayModel>();
    ToNdMatrix<4, VprFloatEntry, float>(&delays, model.getDelays(), ToFloat);

    base_delay_model_ = std::make_unique<DeltaDelayModel>(cross_layer_delay_, delays, is_flat_);

    // Reading non-scalar capnproto fields is roughly equivilant to using
    // a std::vector of the field type.  Actual type is capnp::List<X>::Reader.
    auto overrides = model.getDelayOverrides();
    std::vector<std::pair<t_override, float> > overrides_arr(overrides.size());
    for (size_t i = 0; i < overrides.size(); ++i) {
        const auto& elem = overrides[i];
        overrides_arr[i].first.from_type = elem.getFromType();
        overrides_arr[i].first.to_type = elem.getToType();
        overrides_arr[i].first.from_class = elem.getFromClass();
        overrides_arr[i].first.to_class = elem.getToClass();
        overrides_arr[i].first.delta_x = elem.getDeltaX();
        overrides_arr[i].first.delta_y = elem.getDeltaY();

        overrides_arr[i].second = elem.getDelay();
    }

    delay_overrides_ = vtr::make_flat_map2(std::move(overrides_arr));
}

void OverrideDelayModel::write(const std::string& file) const {
    ::capnp::MallocMessageBuilder builder;
    auto model = builder.initRoot<VprOverrideDelayModel>();

    auto delays = model.getDelays();
    FromNdMatrix<4, VprFloatEntry, float>(&delays, base_delay_model_->delays(), FromFloat);

    // Non-scalar capnproto fields should be first initialized with
    // init<field  name>(count), and then accessed from the returned
    // std::vector-like Builder object (specifically capnp::List<X>::Builder).
    auto overrides = model.initDelayOverrides(delay_overrides_.size());
    auto dst_iter = overrides.begin();
    for (const auto& src : delay_overrides_) {
        auto elem = *dst_iter++;
        elem.setFromType(src.first.from_type);
        elem.setToType(src.first.to_type);
        elem.setFromClass(src.first.from_class);
        elem.setToClass(src.first.to_class);
        elem.setDeltaX(src.first.delta_x);
        elem.setDeltaY(src.first.delta_y);

        elem.setDelay(src.second);
    }

    writeMessageToFile(file, &builder);
}

#endif

///@brief Initialize the placer delay model.
std::unique_ptr<PlaceDelayModel> alloc_lookups_and_delay_model(const Netlist<>& net_list,
                                                               t_chan_width_dist chan_width_dist,
                                                               const t_placer_opts& placer_opts,
                                                               const t_router_opts& router_opts,
                                                               t_det_routing_arch* det_routing_arch,
                                                               std::vector<t_segment_inf>& segment_inf,
                                                               const std::vector<t_direct_inf>& directs,
                                                               bool is_flat) {
    return compute_place_delay_model(placer_opts,
                                     router_opts,
                                     net_list,
                                     det_routing_arch,
                                     segment_inf,
                                     chan_width_dist,
                                     directs,
                                     is_flat);
}

bool written = false;

/**
 * @brief Returns the delay of one point to point connection.
 *
 * Only estimate delay for signals routed through the inter-block routing network.
 * TODO: Do how should we compute the delay for globals. "Global signals are assumed to have zero delay."
 */
float comp_td_single_connection_delay(const PlaceDelayModel* delay_model,
                                      const vtr::vector_map<ClusterBlockId, t_block_loc>& block_locs,
                                      ClusterNetId net_id,
                                      int ipin, float layer_weight) {
    auto& cluster_ctx = g_vpr_ctx.clustering();

    float delay_source_to_sink = 0.;

    if (!cluster_ctx.clb_nlist.net_is_ignored(net_id)) {
        ClusterPinId source_pin = cluster_ctx.clb_nlist.net_driver(net_id);
        ClusterPinId sink_pin = cluster_ctx.clb_nlist.net_pin(net_id, ipin);

        ClusterBlockId source_block = cluster_ctx.clb_nlist.pin_block(source_pin);
        ClusterBlockId sink_block = cluster_ctx.clb_nlist.pin_block(sink_pin);

        int source_block_ipin = cluster_ctx.clb_nlist.pin_logical_index(source_pin);
        int sink_block_ipin = cluster_ctx.clb_nlist.pin_logical_index(sink_pin);

        t_pl_loc source_block_loc = block_locs[source_block].loc;
        t_pl_loc sink_block_loc = block_locs[sink_block].loc;

        // auto source_physical_type = g_vpr_ctx.device().grid.get_physical_type(t_physical_tile_loc(source_block_loc.x, source_block_loc.y, source_block_loc.layer));
        // auto sink_physical_type = g_vpr_ctx.device().grid.get_physical_type(t_physical_tile_loc(sink_block_loc.x, sink_block_loc.y, sink_block_loc.layer));

        // auto source_pin_logical_index = cluster_ctx.clb_nlist.pin_logical_index(source_pin);
        // auto sink_pin_logical_index = cluster_ctx.clb_nlist.pin_logical_index(sink_pin);

        // bool is_pin_on_side = sink_physical_type->pinloc[0][0][1][sink_pin_logical_index];
        // bool is_pin_on_side_2 = source_physical_type->pinloc[0][0][1][source_pin_logical_index];

        e_side source_side = TOP; // doesn't matter for now, value is not used in delay model
        e_side sink_side = TOP; // doesn't matter for now, value is not used in delay model

        // if(sink_block_loc.x == 0) sink_side = RIGHT;
        // else if(sink_block_loc.x == g_vpr_ctx.device().grid.width() - 1) sink_side = LEFT;
        // else if(sink_block_loc.y == 0) sink_side = TOP;
        // else if(sink_block_loc.y == g_vpr_ctx.device().grid.height() - 1) sink_side = BOTTOM;
        // else {
        //     for (e_side side : {TOP, BOTTOM, LEFT, RIGHT}) {
        //         if (sink_physical_type->pinloc[0][0][(int)side][sink_pin_logical_index]) {
        //             sink_side = side;
        //         }
        //     }
        // }

        // if(source_block_loc.x == 0) source_side = RIGHT;
        // else if(source_block_loc.x == g_vpr_ctx.device().grid.width() - 1) source_side = LEFT;
        // else if(source_block_loc.y == 0) source_side = TOP;
        // else if(source_block_loc.y == g_vpr_ctx.device().grid.height() - 1) source_side = BOTTOM;
        // else {
        //     for (e_side side : {TOP, BOTTOM, LEFT, RIGHT}) {
        //         if (source_physical_type->pinloc[0][0][(int)side][source_block_ipin]) {
        //             source_side = side;
        //         }
        //     }
        // }

        /**
         * This heuristic only considers delta_x and delta_y, a much better
         * heuristic would be to to create a more comprehensive lookup table.
         *
         * In particular this approach does not accurately capture the effect
         * of fast carry-chain connections.
         */
        // delay_source_to_sink = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
        //                                           {sink_block_loc.x, sink_block_loc.y, sink_block_loc.layer}, sink_block_ipin, source_side, sink_side);

        if (source_block_loc.layer != sink_block_loc.layer) {
            // Account for the fact that the delay model is not exact when crossing layers.
            float delay_2d = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
                                                  {sink_block_loc.x, sink_block_loc.y, source_block_loc.layer}, sink_block_ipin, source_side, sink_side);
            float delay_3d = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
                                                    {sink_block_loc.x, sink_block_loc.y, sink_block_loc.layer}, sink_block_ipin, source_side, sink_side);
            float inter_layer_delay = delay_3d - delay_2d;
            
            if (delay_2d > delay_3d){
                inter_layer_delay = 0;
            }

            float blended_delay = delay_2d + (inter_layer_delay * layer_weight);
            
            delay_source_to_sink = blended_delay;
        } else{
        delay_source_to_sink = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
                                                  {sink_block_loc.x, sink_block_loc.y, sink_block_loc.layer}, sink_block_ipin, source_side, sink_side);
        }
        /*

            I imagine if I want to do a transitioning delay model, I would do something like this, where layer_weight is between 0 and 1, and increases from 0 to 1 as the annealing process is going along. Finally capping at 1.0 so the final delay_model is exact. <3:

        delay_2d = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
                                                  {sink_block_loc.x, sink_block_loc.y, source_block_loc.layer}, sink_block_ipin, source_side, sink_side);
        delay_3d = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
                                                  {sink_block_loc.x, sink_block_loc.y, sink_block_loc.layer}, sink_block_ipin, source_side, sink_side);
        float inter_layer_delay = delay_3d - delay_2d;    

        layer_weight = get_layer_weight();
        blended_delay = delay_2d + inter_layer_delay * layer_weight;
        
        delay_source_to_sink = blended_delay;
        */

        if (delay_source_to_sink < 0) {
            float delay_2d = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
                                                  {sink_block_loc.x, sink_block_loc.y, source_block_loc.layer}, sink_block_ipin, source_side, sink_side);
            float delay_3d = delay_model->delay({source_block_loc.x, source_block_loc.y, source_block_loc.layer}, source_block_ipin,
                                                    {sink_block_loc.x, sink_block_loc.y, sink_block_loc.layer}, sink_block_ipin, source_side, sink_side);
            VTR_LOG("delay 2d: %g, delay_3d: %g, inter_layer_delay: %g, layer_weight: %g, blended_delay: %g\n",
                    delay_2d, delay_3d, delay_3d - delay_2d, layer_weight, delay_source_to_sink);
            VPR_ERROR(VPR_ERROR_PLACE,
                      "in comp_td_single_connection_delay: Bad delay_source_to_sink value %g from %s (at %d,%d,%d) to %s (at %d,%d,%d)\n"
                      "in comp_td_single_connection_delay: Delay is less than 0\n",
                      block_type_pin_index_to_name(physical_tile_type(source_block_loc), source_block_ipin, false).c_str(),
                      source_block_loc.x, source_block_loc.y, source_block_loc.layer,
                      block_type_pin_index_to_name(physical_tile_type(sink_block_loc), sink_block_ipin, false).c_str(),
                      sink_block_loc.x, sink_block_loc.y, sink_block_loc.layer,
                      delay_source_to_sink);
        }
        // For debugging purposes, write the delay model to file once 
        if (!written){
            delay_model->write(std::string("delay_model.txt"));
            written = true;
        }
            
    }

    return (delay_source_to_sink);
}

///@brief Recompute all point to point delays, updating `connection_delay` matrix.
void comp_td_connection_delays(const PlaceDelayModel* delay_model,
                               PlacerState& placer_state, float layer_weight) {
    const auto& cluster_ctx = g_vpr_ctx.clustering();
    auto& p_timing_ctx = placer_state.mutable_timing();
    auto& block_locs = placer_state.block_locs();
    auto& connection_delay = p_timing_ctx.connection_delay;

    for (ClusterNetId net_id : cluster_ctx.clb_nlist.nets()) {
        for (size_t ipin = 1; ipin < cluster_ctx.clb_nlist.net_pins(net_id).size(); ++ipin) {
            connection_delay[net_id][ipin] = comp_td_single_connection_delay(delay_model, block_locs, net_id, ipin, layer_weight);
        }
    }
}
