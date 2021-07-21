#include <XML/Helper.h>
#include "TMath.h"
#include "TString.h"
#include "GeometryHelpers.h"
#include "Math/Point2D.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;


// builders
namespace sphere_test_builder {

  void build_sphere(Detector &desc, DetElement &sdet, Volume &env, xml::Component plm) {

    // attributes
    auto mat = desc.material(plm.attr<std::string>(_Unicode(material)));
    auto vis = desc.visAttributes(plm.attr<std::string>(_Unicode(vis)));
    double centerx = plm.attr<double>(_Unicode(centerx));
    double centery = plm.attr<double>(_Unicode(centery));
    double centerz = plm.attr<double>(_Unicode(centerz));
    double radius = plm.attr<double>(_Unicode(radius));
    double thickness = plm.attr<double>(_Unicode(thickness));
    double thetamin = plm.attr<double>(_Unicode(thetamin));
    double thetamax = plm.attr<double>(_Unicode(thetamax));
    double phimin = plm.attr<double>(_Unicode(phimin));
    double phimax = plm.attr<double>(_Unicode(phimax));
    double yrot = plm.attr<double>(_Unicode(yrot));

    // optical surface
    //OpticalSurfaceManager surfMgr = desc.surfaceManager();
    //auto surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(plm, _Unicode(surface), "MirrorOpticalSurface"));

    // volume
    Volume vol("mirror_v");
    vol.setMaterial(mat);
    vol.setVisAttributes(vis);
    vol.setSolid(Sphere( radius-thickness, radius, thetamin, thetamax, phimin, phimax ));

    // placement
    auto pv = env.placeVolume(vol,
      Translation3D(0,0,0) * RotationY(yrot)
    );

    // apply surface properties
    DetElement de(sdet, "mirror_de", 0);
    de.setPlacement(pv);
    //SkinSurface skin(desc, de, "mirror_optical_surface", surf, vol);
    //skin.isValid();
  };


  void build_box(Detector &desc, DetElement &sdet, Volume &env, xml::Component plm, SensitiveDetector &sens) {

    // attributes
    auto mat = desc.material(plm.attr<std::string>(_Unicode(material)));
    auto vis = desc.visAttributes(plm.attr<std::string>(_Unicode(vis)));
    double centerx = plm.attr<double>(_Unicode(centerx));
    double centery = plm.attr<double>(_Unicode(centery));
    double centerz = plm.attr<double>(_Unicode(centerz));
    double lx = plm.attr<double>(_Unicode(halflengthx));
    double ly = plm.attr<double>(_Unicode(halflengthy));
    double lz = plm.attr<double>(_Unicode(halflengthz));

    // volume
    Volume vol("sensor_v");
    vol.setMaterial(mat);
    vol.setVisAttributes(vis);
    vol.setSensitiveDetector(sens);
    vol.setSolid(Box(lx,ly,lz));

    // placement
    auto pv = env.placeVolume(vol,
      Translation3D(0,0,0) * RotationX(0)
    );

    pv.addPhysVolID("sensor", 0);
    DetElement de(sdet, "sensor_de", 0);
    de.setPlacement(pv);
  };

};

// ====================================================


// create the detector
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {

  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();

  DetElement det(detName, detID);
  xml::Component dims = detElem.dimensions();

  // envelope
  double lx = dims.attr<double>(_Unicode(halflengthx));
  double ly = dims.attr<double>(_Unicode(halflengthy));
  double lz = dims.attr<double>(_Unicode(halflengthz));
  Box envShape(lx,ly,lz);
  auto envMat = desc.material(dd4hep::getAttrOrDefault(detElem, _Unicode(material), "AirOptical"));
  Volume envVol(detName + "_envelope", envShape, envMat);
  envVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

  // children
  sphere_test_builder::build_sphere(desc, det, envVol, detElem.child(_Unicode(sphere)));
  sphere_test_builder::build_box(desc, det, envVol, detElem.child(_Unicode(sensor)),sens);

  // place mother volume
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume envPV = motherVol.placeVolume(envVol, Position(0, 0, 0));
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);

  return det;
};


//@}

// clang-format off
DECLARE_DETELEMENT(athena_sphere_test, create_detector)

