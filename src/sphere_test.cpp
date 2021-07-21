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

// create the detector
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {

  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();

  DetElement det(detName, detID);
  xml::Component dims = detElem.dimensions();

  // envelope
  double envx = dims.attr<double>(_Unicode(halflengthx));
  double envy = dims.attr<double>(_Unicode(halflengthy));
  double envz = dims.attr<double>(_Unicode(halflengthz));
  Box envShape(envx,envy,envz);
  auto envMat = desc.material(dd4hep::getAttrOrDefault(detElem, _Unicode(material), "AirOptical"));
  Volume envVol(detName + "_envelope", envShape, envMat);
  envVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

  // build sphere --------------------------
  // attributes
  auto sphereMat = desc.material(detElem.child(_Unicode(sphere)).attr<std::string>(_Unicode(material)));
  auto sphereVis = desc.visAttributes(detElem.child(_Unicode(sphere)).attr<std::string>(_Unicode(vis)));
  double sphx = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(centerx));
  double sphy = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(centery));
  double sphz = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(centerz));
  double radius = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(radius));
  double thickness = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(thickness));
  double thetamin = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(thetamin));
  double thetamax = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(thetamax));
  double phimin = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(phimin));
  double phimax = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(phimax));
  double yrot = detElem.child(_Unicode(sphere)).attr<double>(_Unicode(yrot));

  // build box --------------------------
  // attributes
  auto boxMat = desc.material(detElem.child(_Unicode(sensor)).attr<std::string>(_Unicode(material)));
  auto boxVis = desc.visAttributes(detElem.child(_Unicode(sensor)).attr<std::string>(_Unicode(vis)));
  double boxx = detElem.child(_Unicode(sensor)).attr<double>(_Unicode(centerx));
  double boxy = detElem.child(_Unicode(sensor)).attr<double>(_Unicode(centery));
  double boxz = detElem.child(_Unicode(sensor)).attr<double>(_Unicode(centerz));
  double lx = detElem.child(_Unicode(sensor)).attr<double>(_Unicode(halflengthx));
  double ly = detElem.child(_Unicode(sensor)).attr<double>(_Unicode(halflengthy));
  double lz = detElem.child(_Unicode(sensor)).attr<double>(_Unicode(halflengthz));

  Sphere mirror_solid( radius-thickness, radius, thetamin, thetamax, phimin, phimax );
  IntersectionSolid mirror_solid2(mirror_solid, Box(lx,ly,lz), Position(0,0,radius-thickness));

  // volume
  Volume sphereVol("mirror_v",mirror_solid2,sphereMat);

  // placement
  auto spherePV = envVol.placeVolume(sphereVol, Transform3D(RotationY(yrot)));
  DetElement sphereDE(det, "mirror_de", 0);
  sphereDE.setPlacement(spherePV);



  Box sensor_solid(lx,ly,lz);
  // volume
  Volume boxVol("sensor_v",sensor_solid,boxMat);
  boxVol.setVisAttributes(boxVis);
  boxVol.setSensitiveDetector(sens);

  // placement
  auto boxPV = envVol.placeVolume(boxVol);
  boxPV.addPhysVolID("sensor", 0);
  DetElement boxDE(det, "sensor_de", 0);
  boxDE.setPlacement(boxPV);


  // place mother volume ---------
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume envPV = motherVol.placeVolume(envVol, Position(0, 0, 0));
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);

  return det;
};


//@}

// clang-format off
DECLARE_DETELEMENT(athena_sphere_test, create_detector)

