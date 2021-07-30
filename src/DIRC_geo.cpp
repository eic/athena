#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
#include "XML/Layering.h"
//////////////////////////////////
// Central Barrel DIRC
//////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  xml_dim_t dim   = x_det.dimensions();
  xml_dim_t pos   = x_det.position();
  double    RIn   = dim.rmin();
  double    ROut  = dim.rmax();
  double    SizeZ = dim.length();

  Material Vacuum = desc.material("Vacuum");
  Material mat = desc.material("Quartz");
  Tube     cb_DIRC_Barrel_GVol_Solid(RIn, ROut, SizeZ / 2.0, 0., 360.0 * deg);
  Volume   detVol("cb_DIRC_GVol_Solid_Logic", cb_DIRC_Barrel_GVol_Solid, Vacuum);
  detVol.setVisAttributes(desc.invisible());

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, pos.z()));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  detPV.addPhysVolID("barrel", 1);
  det.setPlacement(detPV);

  //////////////////
  // DIRC Bars
  //////////////////

  double dR                = dim.radius();
  double cb_DIRC_bars_DZ   = SizeZ; 
  double cb_DIRC_bars_DY   = dim.dy();
  double cb_DIRC_bars_DX   = dim.dx();
  int    NUM               = dim.number();
  double cb_DIRC_bars_DPhi = dim.deltaphi();

  Material cb_DIRC_bars_Material = desc.material("Quartz");

  Box    cb_DIRC_bars_Solid(cb_DIRC_bars_DX / 2., cb_DIRC_bars_DY / 2., cb_DIRC_bars_DZ / 2.);
  Volume cb_DIRC_bars_Logic("cb_DIRC_bars_Logix", cb_DIRC_bars_Solid, cb_DIRC_bars_Material);
  cb_DIRC_bars_Logic.setVisAttributes(desc.visAttributes(x_det.visStr()));
  sens.setType("photoncounter");
  cb_DIRC_bars_Logic.setSensitiveDetector(sens);

  int count = 0;
  for (xml_coll_t mod(x_det, _U(module)); mod; ++mod) {
    xml_comp_t x_mod = mod;
    Transform3D  tr(RotationZ(x_mod.phi()), Position(-x_mod.R() * cos(x_mod.phi()), -x_mod.R() * sin(x_mod.phi()), 0));
    PlacedVolume barPV = detVol.placeVolume(cb_DIRC_bars_Logic, tr);
    barPV.addPhysVolID("module", count++);
  }
  
  return det;
}

DECLARE_DETELEMENT(DIRC, createDetector)
