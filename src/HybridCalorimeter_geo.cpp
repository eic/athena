//==========================================================================
//  A general implementation for homogeneous calorimeter
//  it supports three types of placements
//  1. Individual module placement with module dimensions and positions
//  2. Array placement with module dimensions and numbers of row and columns
//  3. Disk placement with module dimensions and (Rmin, Rmax), and (Phimin, Phimax)
//  4. Lines placement with module dimensions and (mirrorx, mirrory)
//     (NOTE: anchor point is the 0th block of the line instead of line center)
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 06/09/2021
//==========================================================================

#include "GeometryHelpers.h"
#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <math.h>
#include <fmt/core.h>

using namespace dd4hep;

// main
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{

  using namespace std;
  using namespace fmt;

    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    DetElement det(detName, detID);
    sens.setType("calorimeter");

    auto glass_material = desc.material("SciGlass");
    auto crystal_material = desc.material("PbWO4");
    auto air_material = desc.material("Air");

    double ROut = desc.constantAsDouble("EcalEndcapN_rmax");
    double RIn = desc.constantAsDouble("EcalEndcapN_rmin");
    double SizeZ = desc.constantAsDouble("EcalEndcapN_thickness");
    double thickness = desc.constantAsDouble("EcalEndcapN_thickness");
    double trans_radius = desc.constantAsDouble("EcalEndcapNCrystal_rmax");
    double Glass_z0 = desc.constantAsDouble("GlassModule_z0");
    double Glass_Width = desc.constantAsDouble("GlassModule_width");
    double Glass_thickness = desc.constantAsDouble("GlassModule_length");
    double Glass_Gap = desc.constantAsDouble("GlassModule_wrap");
    double glass_distance = desc.constantAsDouble("GlassModule_distance");

    double Crystal_Width = desc.constantAsDouble("CrystalModule_width");
    double Crystal_thickness = desc.constantAsDouble("CrystalModule_length");
    double Crystal_Gap = desc.constantAsDouble("CrystalModule_wrap");
    double crystal_distance = desc.constantAsDouble("CrystalModule_distance");
    double Crystal_z0 = desc.constantAsDouble("CrystalModule_z0");

    // RIn and ROut will define outer tube embedding the calorimeter
    // centers_rin/out define the maximum radius of module centers
    // so that modules are not overlapping with mother tube volume
    double hypotenuse = sqrt(0.5 * glass_distance * glass_distance);
    double centers_rin = RIn + hypotenuse + 1*mm;
    double centers_rout = ROut - hypotenuse - 1*mm;

    const double Crystal_offset = -0.5*(Crystal_thickness - thickness);
    const double Glass_offset = -0.5*(Glass_thickness - thickness);

    // envelope

    // Assembly assembly(detName);

    Tube outer_tube(RIn, ROut, SizeZ / 2.0, 0., 360.0 * deg);
    Volume ecal_vol("negative_ecal", outer_tube, air_material);
    ecal_vol.setVisAttributes(desc.visAttributes("HybridEcalOuterVis"));

    double Glass_OuterR = ROut - 1 * cm ;
    double Glass_InnerR = trans_radius;

    // Geometry of modules
    Box glass_box("glass_box", Glass_Width * 0.5, Glass_Width * 0.5, Glass_thickness * 0.5);
    Volume glass_module("glass_module", glass_box, glass_material);
    glass_module.setVisAttributes(desc.visAttributes("EcalEndcapNModuleVis"));
    glass_module.setSensitiveDetector(sens);
    
    Box crystal_box("crystal_box",  Crystal_Width* 0.5, Crystal_Width * 0.5, Crystal_thickness * 0.5);
    Volume crystal_module("crystal_module", crystal_box, crystal_material);
    crystal_module.setVisAttributes(desc.visAttributes("EcalEndcapNModuleVis"));
    crystal_module.setSensitiveDetector(sens);

    // GLASS
    double diameter = 2 * Glass_OuterR;

    // How many towers do we have per row/columnt?
    // Add a gap + diameter as if we have N towers, we have N-1 gaps;
    int towersInRow = std::ceil((diameter + Glass_Gap) /  (Glass_Width + Glass_Gap));

    // Is it odd or even number of towersInRow
    double leftTowerPos, topTowerPos;
    if(towersInRow%2) {
      //             |
      //      [ ][ ][ ][ ][ ]
      //       ^     |
      int towersInHalfRow = std::ceil(towersInRow/2.0);
      topTowerPos = leftTowerPos = -towersInHalfRow * (Glass_Width + Glass_Gap);

    } else {
      //               |
      //      [ ][ ][ ][ ][ ][ ]
      //       ^      |
      int towersInHalfRow = towersInRow/2;
      topTowerPos = leftTowerPos = -(towersInHalfRow - 0.5) * (Glass_Width + Glass_Gap);
    }

    int moduleIndex = 0;

//    fmt::print("\nCE EMCAL GLASS SQUARE START\n");
//    fmt::print("Glass_thickness = {} cm;\n", Glass_thickness / cm);
//    fmt::print("Glass_Width     = {} cm;\n", Glass_Width / cm);
//    fmt::print("Glass_Gap       = {} cm;\n", Glass_Gap / cm);
//    fmt::print("Glass_InnerR    = {} cm;\n", Glass_InnerR / cm);
//    fmt::print("Glass_OuterR    = {} cm;\n", Glass_OuterR / cm);
//    fmt::print("Glass_PosZ      = {} cm;\n", glass_shift_z / cm);
//    fmt::print("Towers in Row/Col   = {} cm;\n", glass_shift_z / cm);
//    fmt::print("Top left tower pos  = {:<10} {:<10} cm;\n", -leftTowerPos / cm, topTowerPos / cm);
// fmt::print("#Towers info:\n");
// fmt::print("#{:<5} {:<6} {:<3} {:<3} {:>10} {:>10}   {}\n", "idx",  "code", "col", "row", "x", "y", "name");


    // We first do a "dry run", not really placing modules,
    // but figuring out the ID scheme, number of modules, etc.
    int glassModuleCount = 0;
    int crystalModuleCount = 0;
    int firstCrystRow = 1000000;   // The first row, where crystals are started
    int firstCrystCol = 1000000;   // The fist column, where crystals are started
    for(int rowIndex=0; rowIndex < towersInRow; rowIndex++) {
      for(int colIndex=0; colIndex < towersInRow; colIndex++) {
        double glass_x = leftTowerPos + colIndex * glass_distance;
        double glass_y = topTowerPos + rowIndex * glass_distance;
        double r = sqrt(glass_x * glass_x + glass_y * glass_y);

        if (r < centers_rout && r > centers_rin) {
          // we are placing something
          if(r<trans_radius) {
            // 4 crystal modules will be placed
            crystalModuleCount+=4;

            // Finding the first col and row where crystals start
            // is the same algorithm as finding a minimum in array
            if(colIndex<firstCrystCol) {
              firstCrystCol = colIndex;
            }
            if(rowIndex<firstCrystRow) {
              firstCrystRow = rowIndex;
            }
          }
          else
          {
            // glass module will be places
            glassModuleCount++;
          }
        }
      }
    }
    // fmt::print("#Towers info:\n");
    // fmt::print("#{:<5} {:<6} {:<3} {:<3} {:>10} {:>10}   {}\n", "idx",  "code", "col", "row", "x", "y", "name");
    int glass_module_index = 0;
    int cryst_module_index = 0;
    for(int rowIndex=0; rowIndex < towersInRow; rowIndex++) {
      for(int colIndex=0; colIndex < towersInRow; colIndex++) {
        double glass_x = leftTowerPos + colIndex * (Glass_Width + Glass_Gap);
        double glass_y = topTowerPos + rowIndex * (Glass_Width + Glass_Gap);
        double r = sqrt(glass_x * glass_x + glass_y * glass_y);

        if (r < centers_rout && r > centers_rin) {
          int code = 1000 * rowIndex + colIndex;
          std::string name = fmt::format("ce_EMCAL_glass_phys_{}", code);

          if(r<trans_radius) {

            // first crystal module
            double crystal_x = glass_x - crystal_distance / 2;
            double crystal_y = glass_y - crystal_distance / 2;
            auto   placement = ecal_vol.placeVolume(crystal_module, Position(crystal_x, crystal_y, Crystal_z0 + Crystal_offset));
            placement.addPhysVolID("sector", 1);
            placement.addPhysVolID("module", cryst_module_index++);

            // second crystal module
            crystal_x = glass_x + crystal_distance / 2;
            crystal_y = glass_y - crystal_distance / 2;
            placement = ecal_vol.placeVolume(crystal_module, Position(crystal_x, crystal_y, Crystal_z0 + Crystal_offset));
            placement.addPhysVolID("sector", 1);
            placement.addPhysVolID("module", cryst_module_index++);

            // third crystal module
            crystal_x = glass_x - crystal_distance / 2;
            crystal_y = glass_y + crystal_distance / 2;
            placement = ecal_vol.placeVolume(crystal_module, Position(crystal_x, crystal_y, Crystal_z0 + Crystal_offset));
            placement.addPhysVolID("sector", 1);
            placement.addPhysVolID("module", cryst_module_index++);

            // forth crystal module
            crystal_x = glass_x + crystal_distance / 2;
            crystal_y = glass_y + crystal_distance / 2;
            placement = ecal_vol.placeVolume(crystal_module, Position(crystal_x, crystal_y, Crystal_z0 + Crystal_offset));
            placement.addPhysVolID("sector", 1);
            placement.addPhysVolID("module", cryst_module_index++);
          }
          else
          {
            // glass module
            auto placement = ecal_vol.placeVolume(glass_module, Position(glass_x, glass_y, Glass_z0 + Glass_offset));
            placement.addPhysVolID("sector", 2);
            placement.addPhysVolID("module", glass_module_index++);
          }


          // fmt::print(" {:<5} {:<6} {:<3} {:<3} {:>10.4f} {:>10.4f}   {}\n", towerIndex, code, colIndex, rowIndex, x / cm, y / cm, name);
          //glass_module_index++;
        }
      }
    }

    desc.add(Constant("EcalEndcapN_NModules_Sector1", std::to_string(cryst_module_index)));
    desc.add(Constant("EcalEndcapN_NModules_Sector2", std::to_string(glass_module_index)));
//    fmt::print("Total Glass modules: {}\n", towerIndex);
//    fmt::print("CE EMCAL GLASS END\n\n");

     // detector position and rotation
     auto pos = detElem.position();
     auto rot = detElem.rotation();
     Volume motherVol = desc.pickMotherVolume(det);
     Transform3D tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
     PlacedVolume envPV = motherVol.placeVolume(ecal_vol, tr);
     envPV.addPhysVolID("system", detID);
     det.setPlacement(envPV);
     return det;

}

//@}
DECLARE_DETELEMENT(HybridCalorimeter, create_detector)

