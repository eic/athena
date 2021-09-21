//==========================================================================
// Implementation for Gridpix detector
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

    xml_det_t xml_det   = handle;
    string    det_name = xml_det.nameStr();
    int       det_id   = xml_det.id();
    xml_dim_t det_dim   = xml_det.dimensions();
    xml_dim_t det_pos   = xml_det.position();
    double    det_rin   = det_dim.rmin();
    double    det_rout  = det_dim.rmax();
    double    det_length = det_dim.length();


    long verbose = desc.constantAsLong("Gridpix_verbose");

    // We could use Tube as a detector outer volume, this would put sub elements strictly
    // hierarchically in Geant4 (it is nice!). But then we can't overlap this outer tube
    // Assembly put elements just in mother volume in Geant4.
    // Tube     det_geo(det_rin, det_rout, det_length / 2., 0., 360.0 * deg);
    // Volume   det_volume(det_geo,

    Assembly   det_volume("Gridpix");
    det_volume.setVisAttributes(desc.visAttributes(xml_det.visStr()));

    DetElement   det(det_name, det_id);
    Volume       mother_vol = desc.pickMotherVolume(det);

    Transform3D  tr(Position(0.0, 0.0, det_pos.z()));
    PlacedVolume det_plvol = mother_vol.placeVolume(det_volume, tr);

    det_plvol.addPhysVolID("system", det_id);
    det.setPlacement(det_plvol);

    for (xml_coll_t layer_iter(xml_det, _U(layer)); layer_iter; ++layer_iter) {
        xml_comp_t layer_xml = layer_iter;
        string     layer_name = layer_xml.nameStr();

        double r_min = layer_xml.rmin();
        double r_max = layer_xml.rmax();
        double length = layer_xml.length();
        if(verbose) {
            print("{}: creating layer {:<15} r_in={} r_out={} len={}\n", det_name, layer_name, r_min, r_max, length );
        }

        bool is_gas_volume = layer_name == "gas_volume";

        if(!is_gas_volume) {
            // This is a normal tube shaped layer, whatever it is: kapton, endcap, whatever...
            Tube layer_tube(r_min, r_max,  length / 2.0, 0., 360.0 * deg);
            Volume layer_vol(layer_name, layer_tube, desc.material(layer_xml.materialStr()));
            layer_vol.setVisAttributes(desc.visAttributes(layer_xml.visStr()));
            det_volume.placeVolume(layer_vol, Position(0, 0, layer_xml.z_offset()));
        } else {
            // This! is a gas volume tube. In should be sensitive.
            // To reconstruct later the results first we need to split it to sub tubes.
            // The problem solved here: in geant4, without additional settings, a particle
            // propagated in gas volume will leave only hits on enter and exit. We need tracks
            // so we divide gas in small layers.
            if(!layer_xml.hasAttr(_Unicode(ndivide))) {
                auto message = format("{} layer must have 'ndivide' attribute", layer_name);
                throw runtime_error(message);
            }

            int ndivide = layer_xml.attr<int>(_Unicode(ndivide));

            double sublayer_thickness = (r_max - r_min)/ndivide;
            for(int i=0; i<ndivide; i++) {
                double sublay_rmin = r_min + sublayer_thickness*i;
                double sublay_rmax = r_min + sublayer_thickness*(i+1);
                string sublay_name = format("{}_{}", layer_name, i);
                Tube layer_tube(sublay_rmin, sublay_rmax,  length / 2.0, 0., 360.0 * deg);
                Volume layer_vol(layer_name, layer_tube, desc.material(layer_xml.materialStr()));
                layer_vol.setVisAttributes(desc.visAttributes(layer_xml.visStr()));
                auto pl = det_volume.placeVolume(layer_vol, Position(0, 0, layer_xml.z_offset()));
                pl.addPhysVolID("sector", i);
                layer_vol.setSensitiveDetector(sens);
            }
        }
    }

    // Cooling boxes
    for (xml_coll_t ci(xml_det, _Unicode(cbox)); ci; ++ci) {
        xml_comp_t cbox_xml = ci;
        string     cbox_name = cbox_xml.nameStr();
        auto cbox_material = desc.material(cbox_xml.materialStr());
        auto cbox_vis = desc.visAttributes(cbox_xml.visStr());
        xml_dim_t cbox_dim   = cbox_xml.dimensions();

        Box cbox_geom(cbox_dim.width()/2, cbox_dim.height()/2, cbox_dim.length()/2);
        Volume layer_vol(cbox_name, cbox_geom, cbox_material);

        if(verbose) {
            print("{}: creating cbox {:<15} width={} height={} length={}\n",
                  det_name, cbox_name, cbox_dim.width(), cbox_dim.height(), cbox_dim.length());
        }

        for (xml_coll_t module_iter(cbox_xml, _U(module)); module_iter; ++module_iter) {
            xml_comp_t module_xml = module_iter;

            if(verbose) {
                print("{}:    placing module x={} y={} z={}\n", det_name, module_xml.x(), module_xml.y(), module_xml.z());
            }
            det_volume.placeVolume(layer_vol, Position(module_xml.x(), module_xml.y(), module_xml.z()));
        }
    }







//
//    cyl = new PHG4CylinderSubsystem("membrane_2", 4);
//    cyl->set_double_param("radius", G4GRIDPIX::inner_radius);
//    cyl->set_double_param("length", G4GRIDPIX::memb_length);
//    cyl->set_double_param("thickness", G4GRIDPIX::outer_radius - G4GRIDPIX::inner_radius);
//    cyl->set_double_param("place_x", 0);
//    cyl->set_double_param("place_y", 0);
//    cyl->set_double_param("place_z", -endcap_z - G4GRIDPIX::memb_length*3/2 - G4GRIDPIX::cbox_length);
//    cyl->set_string_param("material", "G4_KAPTON");
//    cyl->SetActive(0);
//    cyl->SuperDetector("GRIDPIX");
//    cyl->OverlapCheck(Enable::GRIDPIX_OVERLAPCHECK);
//    g4Reco->registerSubsystem(cyl);
//
//    // C box part of membrane
//    const int n_cbox = 4;
//    double cbox_x[n_cbox] = {G4GRIDPIX::inner_radius + G4GRIDPIX::cbox_length / 2,
//                             G4GRIDPIX::outer_radius - G4GRIDPIX::cbox_length / 2,
//                             -G4GRIDPIX::inner_radius - G4GRIDPIX::cbox_length / 2,
//                             -G4GRIDPIX::outer_radius + G4GRIDPIX::cbox_length / 2};
//
//    // place c-box
//    for(int ibox = 0; ibox < n_cbox; ibox++)
//    {
//        double cbox_z = -endcap_z - G4GRIDPIX::memb_length - G4GRIDPIX::cbox_length/2;
//        auto *cbox = new PHG4BlockSubsystem(Form("CBOX_1_%d",ibox));
//        cbox->set_double_param("size_x", 2.0);
//        cbox->set_double_param("size_y", 2.0);
//        cbox->set_double_param("size_z", 2.0);
//        cbox->set_double_param("place_x", cbox_x[ibox]);
//        cbox->set_double_param("place_y", 0);
//        cbox->set_double_param("place_z", cbox_z);
//        cbox->set_string_param("material", "CF");
//        cbox->OverlapCheck(Enable::GRIDPIX_OVERLAPCHECK);
//        g4Reco->registerSubsystem(cbox);
//    }
//
//    endcap_z = G4GRIDPIX::length/2 - endcap_2_length;
//    // endcap ; gridpix side
//    cyl = new PHG4CylinderSubsystem("ECAP_LAYER_1", 5);
//    cyl->set_double_param("radius", G4GRIDPIX::inner_radius);
//    cyl->set_double_param("length", G4GRIDPIX::ecap_length);
//    cyl->set_double_param("thickness", G4GRIDPIX::outer_radius - G4GRIDPIX::inner_radius);
//    cyl->set_double_param("place_x", 0);
//    cyl->set_double_param("place_y", 0);
//    cyl->set_double_param("place_z", endcap_z + G4GRIDPIX::ecap_length/2);
//    cyl->set_string_param("material", "FR4");
//    cyl->SuperDetector("GRIDPIX");
//    cyl->OverlapCheck(Enable::GRIDPIX_OVERLAPCHECK);
//    g4Reco->registerSubsystem(cyl);
//
//    cyl = new PHG4CylinderSubsystem("ECAP_LAYER_2", 6);
//    cyl->set_double_param("radius", G4GRIDPIX::inner_radius);
//    cyl->set_double_param("length", G4GRIDPIX::ecap_length);
//    cyl->set_double_param("thickness", G4GRIDPIX::outer_radius - G4GRIDPIX::inner_radius);
//    cyl->set_double_param("place_x", 0);
//    cyl->set_double_param("place_y", 0);
//    cyl->set_double_param("place_z", endcap_z + G4GRIDPIX::cbox_length + G4GRIDPIX::ecap_length*3/2);
//    cyl->set_string_param("material", "FR4");
//    cyl->SuperDetector("GRIDPIX");
//    cyl->OverlapCheck(Enable::GRIDPIX_OVERLAPCHECK);
//    g4Reco->registerSubsystem(cyl);
//
//    for(int ibox = 0; ibox < n_cbox; ibox++)
//    {
//        double cbox_z = endcap_z + G4GRIDPIX::ecap_length + G4GRIDPIX::cbox_length/2;
//        auto *cbox = new PHG4BlockSubsystem(Form("CBOX_2_%d",ibox));
//        cbox->set_double_param("size_x", 2.0);
//        cbox->set_double_param("size_y", 2.0);
//        cbox->set_double_param("size_z", 2.0);
//        cbox->set_double_param("place_x", cbox_x[ibox]);
//        cbox->set_double_param("place_y", 0);
//        cbox->set_double_param("place_z", cbox_z);
//        cbox->set_string_param("material", "CF");
//        cbox->OverlapCheck(Enable::GRIDPIX_OVERLAPCHECK);
//        g4Reco->registerSubsystem(cbox);
//    }
//
//
//    return G4GRIDPIX::outer_radius;
//
//     // detector position and rotation
//     auto pos = detElem.position();
//     auto rot = detElem.rotation();
//     Volume motherVol = desc.pickMotherVolume(det);
//     Transform3D tr(RotationZYX(rot.z(), rot.y(), rot.x()), Position(pos.x(), pos.y(), pos.z()));
//     PlacedVolume envPV = motherVol.placeVolume(ecal_vol, tr);
//     envPV.addPhysVolID("system", detID);
//     det.setPlacement(envPV);
//     return det;
    return det;
}

//@}
DECLARE_DETELEMENT(athena_Gridpix, create_detector)

