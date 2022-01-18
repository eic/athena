//==========================================================================
//  Scintillating fiber calorimeter with tower shape blocks
//  reference: https://github.com/adamjaro/lmon/blob/master/calo/src/WScFiZXv3.cxx
//  Support disk placement
//--------------------------------------------------------------------------
//  Author: Chao Peng (ANL)
//  Date: 07/19/2021
//==========================================================================

#include "GeometryHelpers.h"
#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <math.h>

using namespace dd4hep;
using Point =  ROOT::Math::XYPoint;

std::tuple<Volume, Position> build_module(const Detector &desc, const xml::Component &mod_x, SensitiveDetector &sens);

// helper function to get x, y, z if defined in a xml component
template<class XmlComp>
Position get_xml_xyz(const XmlComp &comp, dd4hep::xml::Strng_t name)
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
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    DetElement det(detName, detID);
    sens.setType("calorimeter");
    auto dim = detElem.dimensions();
    auto rmin = dim.rmin();
    auto rmax = dim.rmax();
    auto length = dim.length();
    auto phimin = dd4hep::getAttrOrDefault<double>(dim, _Unicode(phimin), 0.);
    auto phimax = dd4hep::getAttrOrDefault<double>(dim, _Unicode(phimax), 2.*M_PI);
    // envelope
    Tube envShape(rmin, rmax, length/2., phimin, phimax);
    Volume env(detName + "_envelope", envShape, desc.material("Air"));
    env.setVisAttributes(desc.visAttributes(detElem.visStr()));

    // build module
    auto [modVol, modSize] = build_module(desc, detElem.child(_Unicode(module)), sens);
    double modSizeR = std::sqrt(modSize.x() * modSize.x() + modSize.y() * modSize.y());
    double assembly_rwidth = modSizeR*2.;
    int nas = int((rmax - rmin) / assembly_rwidth) + 1;
    std::vector<Assembly> assemblies;
    // calorimeter block z-offsets (as blocks are shorter than the volume length)
    const double block_offset = -0.5*(length - modSize.z());
    for (int i = 0; i < nas; ++i) {
        Assembly assembly(detName + Form("_ring%d", i + 1));
        auto assemblyPV = env.placeVolume(assembly, Position{0., 0., block_offset});
        assemblyPV.addPhysVolID("ring", i + 1);
        assemblies.emplace_back(std::move(assembly));
    }

    int modid = 1;
    for (int ix = 0; ix < int(2.*rmax / modSize.x()) + 1; ++ix) {
        double mx = modSize.x() * ix - rmax;
        for (int iy = 0; iy < int(2.*rmax / modSize.y()) + 1; ++iy) {
            double my = modSize.y() * iy - rmax;
            double mr = std::sqrt(mx*mx + my*my);
            if (mr - modSizeR >= rmin && mr + modSizeR <= rmax) {
                int ias = int((mr - rmin) / assembly_rwidth);
                auto &assembly = assemblies[ias];
                auto modPV = assembly.placeVolume(modVol, Position(mx, my, 0.));
                modPV.addPhysVolID("module", modid++);
            }
        }
    }

    desc.add(Constant(detName + "_NModules", std::to_string(modid - 1)));

    for (auto &assembly : assemblies) {
        assembly.ptr()->Voxelize("");
    }

    // detector position and rotation
    auto pos = get_xml_xyz(detElem, _Unicode(position));
    auto rot = get_xml_xyz(detElem, _Unicode(rotation));
    Volume motherVol = desc.pickMotherVolume(det);
    Transform3D tr = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
    PlacedVolume envPV = motherVol.placeVolume(env, tr);
    envPV.addPhysVolID("system", detID);
    det.setPlacement(envPV);
    return det;
}

// helper function to build module with scintillating fibers
std::tuple<Volume, Position> build_module(const Detector &desc, const xml::Component &mod_x, SensitiveDetector &sens)
{
    auto sx = mod_x.attr<double>(_Unicode(sizex));
    auto sy = mod_x.attr<double>(_Unicode(sizey));
    auto sz = mod_x.attr<double>(_Unicode(sizez));
    auto modMat = desc.material(mod_x.attr<std::string>(_Unicode(material)));

    auto fiber_x  = mod_x.child(_Unicode(fiber));
    auto fr       = fiber_x.attr<double>(_Unicode(radius));
    auto fsx      = fiber_x.attr<double>(_Unicode(spacex));
    auto fsy      = fiber_x.attr<double>(_Unicode(spacey));
    auto foff     = dd4hep::getAttrOrDefault<double>(fiber_x, _Unicode(offset), 0.5*mm);
    auto fiberMat = desc.material(fiber_x.attr<std::string>(_Unicode(material)));

    // Fibers are placed in a honeycomb with the radius = sqrt(3)/2. * hexagon side length
    // So each fiber is fully contained in a regular hexagon, which are placed as
    //           ______________________________________
    //           |          ____        ____          |
    // even:     |         /    \      /    \         |
    //           |    ____/      \____/      \____    |
    //           |   /    \      /    \      /    \   |
    // odd:      |  /      \____/      \____/      \  |
    //           |  \      /    \      /    \      /  |
    //           |   \____/      \____/      \____/   |
    // even:     |        \      /    \      /        |
    //           |         \____/      \____/      ___|___
    //           |____________________________________|___offset
    //                                              | |
    //                                              |offset
    // the parameters space x and space y are used to add additional spaces between the hexagons
    double fside  = 2. / std::sqrt(3.) * fr;
    double fdistx = 2. * fside + fsx;
    double fdisty = 2. * fr + fsy;

    // maximum numbers of the fibers, help narrow the loop range
    int nx = int(sx / (2.*fr)) + 1;
    int ny = int(sy / (2.*fr)) + 1;

    // Instead of placing fibers, we make the whole module as a scintillator
    // and use daughters to cover the insensitive area by radiators
    Box modShape(sx/2., sy/2., sz/2.);
    Volume modVol("module_vol", modShape, fiberMat);
    modVol.setSensitiveDetector(sens);
    if (mod_x.hasAttr(_Unicode(vis)))
        modVol.setVisAttributes(desc.visAttributes(mod_x.attr<std::string>(_Unicode(vis))));

    Box fiberOuterShape(fdistx/2., fdisty/2., sz/2.);
    Tube fiberInnerShape(0., fr, sz/2.);
    SubtractionSolid fiberShape(fiberOuterShape, fiberInnerShape);
    Volume fiberVol("fiber_vol", fiberShape, modMat);

    double y0 = foff + fside;
    double yb = y0 - fdisty/2.;
    double yt = sy - yb - fdisty * std::floor((sy - yb*2.) / fdisty);
    Box modBottomShape(sx/2., yb/2., sz/2.);
    Box modTopShape(sx/2., yt/2., sz/2.);
    Volume modBottomVol("modBottom_vol", modBottomShape, modMat);
    Volume modTopVol("modTop_vol", modTopShape, modMat);

    double x0[2] = {foff + fside, foff + fside + fdistx/2.};
    double xl[2], xr[2];
    Box modLeftShape[2], modRightShape[2];
    Volume modLeftVol[2], modRightVol[2];
    for (int ieo = 0; ieo < 2; ++ieo) {
        xl[ieo] = x0[ieo] - fdistx/2.;
        xr[ieo] = sx - xl[ieo] - fdistx * std::floor((sx - xl[ieo]*2.) / fdistx);
        modLeftShape[ieo] = Box(xl[ieo]/2., fdisty/2., sz/2.);
        modRightShape[ieo] = Box(xr[ieo]/2., fdisty/2., sz/2.);
        modLeftVol[ieo] = Volume(Form("modLeft%d_vol",ieo), modLeftShape[ieo], modMat);
        modRightVol[ieo] = Volume(Form("modRight%d_vol",ieo), modRightShape[ieo], modMat);
    }

    // place the fibers
    int nfibers = 0, nleft = 0, nright = 0;
    modVol.placeVolume(modBottomVol, 0, Position{0, -sy/2.+yb/2., 0});
    for (int iy = 0; iy < ny; ++iy) {
        double y = y0 + fdisty * iy;
        // about to touch the boundary
        if (sy - y < y0) {
            modVol.placeVolume(modTopVol, 0, Position{0, sy/2.-yt/2., 0});
            break;
        }
        int ieo = iy % 2;
        modVol.placeVolume(modLeftVol[ieo], nleft++, Position{-sx/2.+xl[ieo]/2., -sy/2.+y, 0});
        for (int ix = 0; ix < nx; ++ix) {
            double x = x0[ieo] + fdistx * ix;
            // about to touch the boundary
            if (sx - x < x0[ieo]) {
                modVol.placeVolume(modRightVol[ieo], nright++, Position{sx/2.-xr[ieo]/2., -sy/2.+y, 0});
                break;
            }
            modVol.placeVolume(fiberVol, nfibers++, Position{-sx/2.+x, -sy/2.+y, 0});
        }
    }

    return std::make_tuple(modVol, Position{sx, sy, sz});
}

DECLARE_DETELEMENT(ScFiCalorimeter, create_detector)
