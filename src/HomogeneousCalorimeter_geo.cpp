//==========================================================================
//  A general implementation for homogeneous calorimeter
//  it supports three types of placements
//  1. Module placement with module dimensions and positions
//  2. Array placement with module dimensions and numbers of row and columns
//  3. Disk placement with module dimensions and (Rmin, Rmax), and (Phimin, Phimax)
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 06/09/2021
//==========================================================================
#include "GeometryHelpers.h"
#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <iostream>
#include <algorithm>
#include <math.h>

using namespace dd4hep;
using namespace dd4hep::detail;

/** \addtogroup calorimeters Calorimeters
 */

/** \addtogroup Homogeneous Calorimeter
 * \brief Type: **HomogeneousCalorimeter**.
 * \author C. Peng
 * \ingroup calorimeters
 *
 *
 * \code
 *   <detector id="1" name="HyCal" type="HomogeneousCalorimeter" readout="EcalHits" vis="GreenVis">
 *     <dimensions shape="box" sizex="120*cm" sizey="120*cm" sizez="46*cm"/>
 *     <position x="0" y="0" z="0"/>
 *     <rotation x="0" y="0" z="0"/>
 *     <placements>
 *       <array nrow="34" ncol="34" sector="1">
 *         <position x="0" y="0" z="-9.73*cm"/>
 *         <module sizex="2.05*cm" sizey="2.05*cm" sizez="18*cm" vis="GreenVis" material="PbWO4"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *         <removal row="16" col="16"/>
 *         <removal row="16" col="17"/>
 *         <removal row="17" col="16"/>
 *         <removal row="17" col="17"/>
 *       </array>
 *       <array nrow="6" ncol="24" sector="2">
 *         <position x="-17*(2.05+0.015)*cm+12*(3.8+0.015)*cm" y="17*(2.05+0.015)*cm+3*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *       <array nrow="24" ncol="6" sector="3">
 *         <position x="17*(2.05+0.015)*cm+3*(3.8+0.015)*cm" y="17*(2.05+0.015)*cm-12*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *       <array nrow="6" ncol="24" sector="4">
 *         <position x="17*(2.05+0.015)*cm-12*(3.8+0.015)*cm" y="-17*(2.05+0.015)*cm-3*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *       <array nrow="24" ncol="6" sector="3">
 *         <position x="-17*(2.05+0.015)*cm-3*(3.8+0.015)*cm" y="-17*(2.05+0.015)*cm+12*(3.8+0.015)*cm" z="0"/>
 *         <module sizex="3.8*cm" sizey="3.8*cm" sizez="45*cm" vis="BlueVis" material="PbGlass"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *       </array>
 *     </placements>
 *   </detector>
 *
 *   <detector id="2" name="SomeBlocks" type="HomogeneousCalorimeter" readout="EcalHits" vis="GreenVis">
 *     <dimensions shape="box" sizex="100*cm" sizey="100*cm" sizez="20.5*cm"/>
 *     <position x="0" y="0" z="30*cm"/>
 *     <rotation x="0" y="0" z="0"/>
 *     <placements>
 *       <blocks sector="1"/>
 *         <module sizex="2.05*cm" sizey="2.05*cm" sizez="20*cm" vis="GreenVis" material="PbWO4"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *         <placement x="1*cm" y="1*cm" z="0" id="1"/>
 *         <placement x="-1*cm" y="1*cm" z="0" id="2"/>
 *         <placement x="1*cm" y="-1*cm" z="0" id="3"/>
 *         <placement x="-1*cm" y="-1*cm" z="0" id="4"/>
 *       </blocks>
 *     </placements>
 *   </detector>
 *
 *   <detector id="2" name="DiskShapeCalorimeter" type="HomogeneousCalorimeter" readout="EcalHits" vis="GreenVis">
 *     <dimensions shape="disk" rmin="25*cm" rmax="125*cm" length="20.5*cm" phimin="0" phimax="360*degree"/>
 *     <position x="0" y="0" z="-30*cm"/>
 *     <rotation x="0" y="0" z="0"/>
 *     <placements>
 *       <disk rmin="25*cm" rmax="125*cm" length="20.5*cm" phimin="0" phimax="360*degree" sector="1"/>
 *         <module sizex="2.05*cm" sizey="2.05*cm" sizez="20*cm" vis="GreenVis" material="PbWO4"/>
 *         <wrapper thickness="0.015*cm" material="Epoxy" vis="WhiteVis"/>
 *     </placements>
 *   </detector>
 * \endcode
 *
 * @{
 */

// headers
static void add_blocks(Detector& desc, Volume &env, xml::Collection_t &plm, SensitiveDetector &sens, int id);
static void add_array(Detector& desc, Volume &env, xml::Collection_t &plm, SensitiveDetector &sens, int id);
static void add_disk(Detector& desc, Volume &env, xml::Collection_t &plm, SensitiveDetector &sens, int id);

// helper function to get x, y, z if defined in a xml component
template<class XmlComp>
Position get_xml_xyz(XmlComp &comp, dd4hep::xml::Strng_t name)
{
    Position pos(0., 0., 0.);
    if (comp.hasChild(name)) {
        auto child = comp.child(name);
        pos.SetX(dd4hep::getAttrOrDefault<double>(child, _Unicode(x), 0.));
        pos.SetY(dd4hep::getAttrOrDefault<double>(child, _Unicode(y), 0.));
        pos.SetZ(dd4hep::getAttrOrDefault<double>(child, _Unicode(z), 0.));
    }
    return pos;
}

// main
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    static const std::string func = "HomogeneousCalorimeter";
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    DetElement det(detName, detID);
    sens.setType("calorimeter");

    // top-level children
    xml::Component dims = detElem.dimensions();

    // build envelop from dimensions
    std::string shape = dd4hep::getAttrOrDefault(dims, _Unicode(shape), "");
    // no shape input, try to determine shape from dimension variables
    if (shape.empty()) {
        if (dims.hasAttr(_Unicode(rmin)) && dims.hasAttr(_Unicode(rmax)) && dims.hasAttr(_Unicode(length))) {
            shape = "disk";
        } else if (dims.hasAttr(_Unicode(sizex)) && dims.hasAttr(_Unicode(sizey)) && dims.hasAttr(_Unicode(sizez))) {
            shape = "box";
        } else {
            std::cerr << func << " Error: Cannot determine shape of the calorimeter. "
                                 "Add shape (box, or disk) into dimensions\n";
            return det;
        }
    }

    Volume envVol(detName + "_envelope");
    std::string filler = dd4hep::getAttrOrDefault(detElem, _Unicode(filler), "Air");
    envVol.setMaterial(desc.material(filler));
    envVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

    // convert to lower case
    std::transform(shape.begin(), shape.end(), shape.begin(), [](unsigned char c){ return std::tolower(c); });

    if (shape == "box") {
        double sx = dims.attr<double>(_Unicode(sizex));
        double sy = dims.attr<double>(_Unicode(sizey));
        double sz = dims.attr<double>(_Unicode(sizez));
        envVol.setSolid(Box(sx/2., sy/2., sz/2.));
    } else if (shape == "disk") {
        double rmin = dims.rmin();
        double rmax = dims.rmax();
        double length = dims.length();
        double phimin = dd4hep::getAttrOrDefault<double>(dims, _Unicode(phimin), 0.);
        double phimax = dd4hep::getAttrOrDefault<double>(dims, _Unicode(phimax), 2*M_PI);
        envVol.setSolid(Tube(rmin, rmax, length/2., phimin, phimax));
    } else {
        std::cerr << func << " Error: Unsupported shape " << shape << ", use box or disk\n";
        return det;
    }

    // module placement
    xml::Component plm = detElem.child(_Unicode(placements));
    int sector = 1;
    for (xml::Collection_t arr(plm, _Unicode(array)); arr; ++arr) {
        add_array(desc, envVol, arr, sens, sector++);
    }
    for (xml::Collection_t mod(plm, _Unicode(blocks)); mod; ++mod) {
        add_blocks(desc, envVol, mod, sens, sector++);
    }
    for (xml::Collection_t disk(plm, _Unicode(disk)); disk; ++disk) {
        add_disk(desc, envVol, disk, sens, sector++);
    }

    // detector position and rotation
    auto pos = get_xml_xyz(detElem, _Unicode(position));
    auto rot = get_xml_xyz(detElem, _Unicode(rotation));
    Volume motherVol = desc.pickMotherVolume(det);
    Transform3D tr = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
    PlacedVolume envPV = motherVol.placeVolume(envVol, tr);
    envPV.addPhysVolID("system", detID);
    det.setPlacement(envPV);
    return det;
}

// helper function to build module with or w/o wrapper
Volume build_module(Detector &desc, xml::Collection_t &plm, SensitiveDetector &sens, Position &dim)
{
    auto mod = plm.child(_Unicode(module));
    auto sx = mod.attr<double>(_Unicode(sizex));
    auto sy = mod.attr<double>(_Unicode(sizey));
    auto sz = mod.attr<double>(_Unicode(sizez));
    dim = Position{sx, sy, sz};
    Box modShape(sx/2., sy/2., sz/2.);
    auto modMat = desc.material(mod.attr<std::string>(_Unicode(material)));
    Volume modVol("module_vol", modShape, modMat);
    modVol.setSensitiveDetector(sens);
    modVol.setVisAttributes(desc.visAttributes(mod.attr<std::string>(_Unicode(vis))));

    // no wrapper
    if (!plm.hasChild(_Unicode(wrapper))) {
        return modVol;
    // build wrapper
    } else {
        auto wrp = plm.child(_Unicode(wrapper));
        auto thickness = wrp.attr<double>(_Unicode(thickness));
        auto wrpMat = desc.material(wrp.attr<std::string>(_Unicode(material)));
        Box wrpShape((sx + thickness)/2., (sy + thickness)/2., sz/2.);
        Volume wrpVol("wrapper_vol", wrpShape, wrpMat);
        wrpVol.placeVolume(modVol, Position(0., 0., 0.));
        wrpVol.setVisAttributes(desc.visAttributes(wrp.attr<std::string>(_Unicode(vis))));
        dim = Position{sx + thickness, sy + thickness, sz};
        return wrpVol;
    }
}

// place array of modules
static void add_array(Detector& desc, Volume &env, xml::Collection_t &plm, SensitiveDetector &sens, int sid)
{
    Position modSize;
    auto modVol = build_module(desc, plm, sens, modSize);
    int sector_id = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
    int nrow = plm.attr<int>(_Unicode(nrow));
    int ncol = plm.attr<int>(_Unicode(ncol));

    // compute array position
    double begx = -modSize.x()*ncol/2. + modSize.x()/2.;
    double begy = modSize.y()*nrow/2. - modSize.y()/2.;

    std::vector<std::pair<int, int>> removals;
    // get the removal list
    for (xml::Collection_t rm(plm, _Unicode(removal)); rm; ++rm) {
        removals.push_back({rm.attr<int>(_Unicode(row)), rm.attr<int>(_Unicode(col))});
    }

    // placement to mother
    auto pos = get_xml_xyz(plm, _Unicode(position));
    auto rot = get_xml_xyz(plm, _Unicode(rotation));
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            if (std::find(removals.begin(), removals.end(), std::pair<int, int>(i, j)) != removals.end()) {
                continue;
            }
            double px = begx + modSize.x()*j;
            double py = begy - modSize.y()*i;
            Transform3D tr = Translation3D(pos.x() + px, pos.y() + py, pos.z())
                           * RotationZYX(rot.z(), rot.y(), rot.x());
            auto modPV = env.placeVolume(modVol, tr);
            modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", i*ncol + j);
        }
    }
}

// place modules, id must be provided
static void add_blocks(Detector& desc, Volume &env, xml::Collection_t &plm, SensitiveDetector &sens, int sid)
{
    Position modSize;
    auto modVol = build_module(desc, plm, sens, modSize);
    int sector_id = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);

    for (xml::Collection_t pl(plm, _Unicode(placement)); pl; ++pl) {
        Position pos(dd4hep::getAttrOrDefault<double>(pl, _Unicode(x), 0.),
                     dd4hep::getAttrOrDefault<double>(pl, _Unicode(y), 0.),
                     dd4hep::getAttrOrDefault<double>(pl, _Unicode(z), 0.));
        Position rot(dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotx), 0.),
                     dd4hep::getAttrOrDefault<double>(pl, _Unicode(roty), 0.),
                     dd4hep::getAttrOrDefault<double>(pl, _Unicode(rotz), 0.));
        auto mid = pl.attr<int>(_Unicode(id));
        Transform3D tr = Translation3D(pos.x(), pos.y(), pos.z())
                       * RotationZYX(rot.z(), rot.y(), rot.x());
        auto modPV = env.placeVolume(modVol, tr);
        modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", mid);
    }
}

// place disk of modules
static void add_disk(Detector& desc, Volume &env, xml::Collection_t &plm, SensitiveDetector &sens, int sid)
{
    Position modSize;
    auto modVol = build_module(desc, plm, sens, modSize);
    int sector_id = dd4hep::getAttrOrDefault<int>(plm, _Unicode(sector), sid);
    double rmin = plm.attr<double>(_Unicode(rmin));
    double rmax = plm.attr<double>(_Unicode(rmax));
    double phimin = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimin), 0.);
    double phimax = dd4hep::getAttrOrDefault<double>(plm, _Unicode(phimax), 2.*M_PI);

    auto points = ref::utils::fillSquares({0., 0.}, modSize.x(), rmin, rmax, phimin, phimax);
    // placement to mother
    auto pos = get_xml_xyz(plm, _Unicode(position));
    auto rot = get_xml_xyz(plm, _Unicode(rotation));
    int mid = 1;
    for (auto &p : points) {
        Transform3D tr = Translation3D(pos.x() + p.x(), pos.y() + p.y(), pos.z())
                       * RotationZYX(rot.z(), rot.y(), rot.x());
        auto modPV = env.placeVolume(modVol, tr);
        modPV.addPhysVolID("sector", sector_id).addPhysVolID("module", mid++);
    }
}
//@}
DECLARE_DETELEMENT(HomogeneousCalorimeter, create_detector)

