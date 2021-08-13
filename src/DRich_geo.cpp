//==========================================================================
//  dRICh: Dual Ring Imaging Cherenkov Detector
//--------------------------------------------------------------------------
//
// Author: Christopher Dilks (Duke University)
//
// - Design Adapted from Standalone Fun4all and GEMC implementations 
//   [ Evaristo Cisbani, Cristiano Fanelli, Alessio Del Dotto, et al. ]
//
//==========================================================================

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
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();

  DetElement det(detName, detID);
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();

  // attributes -----------------------------------------------------------
  // TODO [low priority]: some attributes have default values, some do not;
  //                      make this more consistent
  // - vessel
  double vesselZ0 = dims.attr<double>(_Unicode(z0));
  double vesselLength = dims.attr<double>(_Unicode(length));
  double vesselRmin0 = dims.attr<double>(_Unicode(rmin0));
  double vesselRmin1 = dims.attr<double>(_Unicode(rmin1));
  double vesselRmax0 = dims.attr<double>(_Unicode(rmax0));
  double vesselRmax1 = dims.attr<double>(_Unicode(rmax1));
  double vesselRmax2 = dims.attr<double>(_Unicode(rmax2));
  double snoutLength = dims.attr<double>(_Unicode(snout_length));
  int nSectors = getAttrOrDefault<int>(dims, _Unicode(nsectors), 6);
  double wallThickness = getAttrOrDefault<double>(dims, _Unicode(wall_thickness), 0.5*cm);
  double windowThickness = getAttrOrDefault<double>(dims, _Unicode(window_thickness), 0.1*cm);
  auto vesselMat = desc.material(getAttrOrDefault(detElem, _Unicode(material), "Aluminum"));
  auto gasvolMat = desc.material(getAttrOrDefault(detElem, _Unicode(gas), "AirOptical"));
  auto vesselVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_vessel)));
  auto gasvolVis = desc.visAttributes(detElem.attr<std::string>(_Unicode(vis_gas)));
  // - radiator (applies to aerogel and filter)
  auto radiatorElem = detElem.child(_Unicode(radiator));
  double radiatorRmin = getAttrOrDefault<double>(radiatorElem, _Unicode(rmin), 10.*cm);
  double radiatorRmax = getAttrOrDefault<double>(radiatorElem, _Unicode(rmax), 80.*cm);
  double radiatorPhiw = getAttrOrDefault<double>(radiatorElem, _Unicode(phiw), 2*M_PI/nSectors);
  double radiatorPitch = getAttrOrDefault<double>(radiatorElem, _Unicode(pitch), 0.*degree);
  double radiatorFrontplane = getAttrOrDefault<double>(radiatorElem, _Unicode(frontplane), 2.5*cm);
  // - aerogel
  auto aerogelElem = radiatorElem.child(_Unicode(aerogel));
  auto aerogelMat = desc.material(aerogelElem.attr<std::string>(_Unicode(material)));
  auto aerogelVis = desc.visAttributes(aerogelElem.attr<std::string>(_Unicode(vis)));
  double aerogelThickness = getAttrOrDefault<double>(aerogelElem, _Unicode(thickness), 2.5*cm);
  // - filter
  auto filterElem = radiatorElem.child(_Unicode(filter));
  auto filterMat = desc.material(filterElem.attr<std::string>(_Unicode(material)));
  auto filterVis = desc.visAttributes(filterElem.attr<std::string>(_Unicode(vis)));
  double filterThickness = getAttrOrDefault<double>(filterElem, _Unicode(thickness), 2.5*cm);
  // - mirror
  auto mirrorElem = detElem.child(_Unicode(mirror));
  auto mirrorMat = desc.material(mirrorElem.attr<std::string>(_Unicode(material)));
  auto mirrorVis = desc.visAttributes(mirrorElem.attr<std::string>(_Unicode(vis)));
  auto mirrorSurf = surfMgr.opticalSurface(getAttrOrDefault(mirrorElem, _Unicode(surface), "MirrorOpticalSurface"));
  double mirrorBackplane = getAttrOrDefault<double>(mirrorElem, _Unicode(backplane), 240.*cm);
  double mirrorThickness = getAttrOrDefault<double>(mirrorElem, _Unicode(thickness), 2.*mm);
  double mirrorRadius = getAttrOrDefault<double>(mirrorElem, _Unicode(radius), 190*cm);
  double mirrorCenterX = getAttrOrDefault<double>(mirrorElem, _Unicode(centerx), 95*cm);
  double mirrorRmin = getAttrOrDefault<double>(mirrorElem, _Unicode(rmin), 10.*cm);
  double mirrorRmax = getAttrOrDefault<double>(mirrorElem, _Unicode(rmax), 150.*cm);
  double mirrorPhiw = getAttrOrDefault<double>(mirrorElem, _Unicode(phiw), 2*M_PI/nSectors);
  int mirrorDebug = getAttrOrDefault<int>(mirrorElem, _Unicode(debug), 0);
  // - sensor module
  auto sensorElem = detElem.child(_Unicode(sensors)).child(_Unicode(module));
  auto sensorMat = desc.material(sensorElem.attr<std::string>(_Unicode(material)));
  auto sensorVis = desc.visAttributes(sensorElem.attr<std::string>(_Unicode(vis)));
  auto sensorSurf = surfMgr.opticalSurface(getAttrOrDefault(sensorElem, _Unicode(surface), "MirrorOpticalSurface"));
  double sensorSide = sensorElem.attr<double>(_Unicode(side));
  double sensorGap = sensorElem.attr<double>(_Unicode(gap));
  double sensorThickness = sensorElem.attr<double>(_Unicode(thickness));
  // - sensor sphere
  auto sensorSphElem = detElem.child(_Unicode(sensors)).child(_Unicode(sphere));
  double sensorSphRadius = sensorSphElem.attr<double>(_Unicode(radius));
  double sensorSphCenterX = sensorSphElem.attr<double>(_Unicode(centerx));
  double sensorSphCenterY = sensorSphElem.attr<double>(_Unicode(centery));
  double sensorSphCenterZ = sensorSphElem.attr<double>(_Unicode(centerz));
  int sensorSphDebug = getAttrOrDefault<int>(sensorSphElem, _Unicode(debug), 0);
  // - sensor sphere patch cuts
  auto sensorSphPatchElem = detElem.child(_Unicode(sensors)).child(_Unicode(sphericalpatch));
  double sensorSphPatchThetaMin = sensorSphPatchElem.attr<double>(_Unicode(thetamin));
  double sensorSphPatchThetaMax = sensorSphPatchElem.attr<double>(_Unicode(thetamax));
  double sensorSphPatchWidthFactor = sensorSphPatchElem.attr<double>(_Unicode(widthfactor));
  double sensorSphPatchTaper = sensorSphPatchElem.attr<double>(_Unicode(taper));

  
  // BUILD VESSEL ====================================================================
  /* - `vessel`: aluminum enclosure, the mother volume of the dRICh
   * - `gasvol`: gas volume, which fills `vessel`; all other volumes defined below
   *   are children of `gasvol`
   * - the dRICh vessel geometry has two regions: the snout refers to the conic region 
   *   in the front, housing the aerogel, while the tank refers to the cylindrical
   *   region, housing the rest of the detector components
   */

  // snout solids
  double boreDelta = vesselRmin1 - vesselRmin0;
  /*
  Cone vesselSnout(
      snoutLength/2.0,
      vesselRmin0,
      vesselRmax0,
      vesselRmin0 + boreDelta * snoutLength / vesselLength,
      vesselRmax1
      );
  Cone gasvolSnout(
      snoutLength/2.0,
      vesselRmin0 + wallThickness,
      vesselRmax0 - wallThickness,
      vesselRmin0 + boreDelta * (snoutLength-windowThickness) / vesselLength + wallThickness,
      vesselRmax1 - wallThickness + windowThickness * (vesselRmax1 - vesselRmax0) / snoutLength
      );

  // tank solids
  Cone vesselTank(
      (vesselLength - snoutLength)/2.0,
      vesselSnout.rMin2(),
      vesselRmax2,
      vesselRmin1,
      vesselRmax2
      );
  Cone gasvolTank(
      (vesselLength - snoutLength)/2.0 - windowThickness,
      gasvolSnout.rMin2(),
      vesselRmax2 - wallThickness,
      vesselRmin1 + wallThickness,
      vesselRmax2 - wallThickness
      );

  // snout + tank solids
  UnionSolid vesselSolid(
      vesselTank,
      vesselSnout,
      Position(0., 0., -vesselLength/2.)
      );
  UnionSolid gasvolSolid(
      gasvolTank,
      gasvolSnout,
      Position(0., 0., -vesselLength/2. + windowThickness)
      );
  */

  // instead of a vessel, make a giant box of C2F6, which contains
  // the interaction point; this is so we can fire opticalphotons 
  // from the interaction point
  Box vesselSolid(500, 500, 500);
  Box gasvolSolid(499, 499, 499);

  // volumes
  Volume vesselVol(detName, vesselSolid, vesselMat);
  Volume gasvolVol(detName + "_gas", gasvolSolid, gasvolMat);
  vesselVol.setVisAttributes(vesselVis);
  gasvolVol.setVisAttributes(gasvolVis);

  // reference position
  auto snoutFront = Position(0., 0., -(vesselLength + snoutLength)/2.);


  // sensitive detector type
  sens.setType("photoncounter");


  // SECTOR LOOP //////////////////////////////////
  for(int isec=0; isec<nSectors; isec++) {

    if(mirrorDebug*isec>0 || sensorSphDebug*isec>0) continue; // if debugging, draw only 1 sector
    double sectorRotation = isec * 360/nSectors * degree; // sector rotation about z axis
    if(isec!=0&&isec!=3) continue; // draw left and right sectors only, for clean OGL top view

    // BUILD RADIATOR ====================================================================

    // derived attributes
    auto radiatorPos = Position(0., 0., radiatorFrontplane) + snoutFront;

    // solid and volume: create aerogel and filter sectors
    Tube aerogelSolid(radiatorRmin, radiatorRmax, aerogelThickness/2, -radiatorPhiw/2.0, radiatorPhiw/2.0);
    Tube filterSolid( radiatorRmin, radiatorRmax, filterThickness/2,  -radiatorPhiw/2.0, radiatorPhiw/2.0);
    Volume aerogelVol("aerogel_v", aerogelSolid, aerogelMat);
    Volume filterVol( "filter_v",  filterSolid,  filterMat);
    aerogelVol.setVisAttributes(aerogelVis);
    filterVol.setVisAttributes(filterVis);

    // placement
    auto aerogelPV = vesselVol.placeVolume(aerogelVol,
          RotationZ(sectorRotation) // rotate about beam axis to sector
        * Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()) // re-center to snoutFront
        * RotationY(radiatorPitch) // change polar angle to specified pitch
        );
    /* // disable filter and surface properties
    auto filterPV = vesselVol.placeVolume(filterVol,
            RotationZ(sectorRotation) // rotate about beam axis to sector
        * Translation3D(radiatorPos.x(), radiatorPos.y(), radiatorPos.z()) // re-center to snoutFront
        * RotationY(radiatorPitch) // change polar angle
        * Translation3D(0., 0., (aerogelThickness+filterThickness)/2.) // move to aerogel backplane
        );

    // properties
    // TODO [critical]: define skin properties for aerogel and filter
    DetElement aerogelDE(det, Form("aerogel_de%d", isec), isec);
    aerogelDE.setPlacement(aerogelPV);
    //SkinSurface aerogelSkin(desc, aerogelDE, Form("mirror_optical_surface%d", isec), aerogelSurf, aerogelVol);
    //aerogelSkin.isValid();
    DetElement filterDE(det, Form("filter_de%d", isec), isec);
    filterDE.setPlacement(filterPV);
    //SkinSurface filterSkin(desc, filterDE, Form("mirror_optical_surface%d", isec), filterSurf, filterVol);
    //filterSkin.isValid();
    */


    // BUILD MIRROR ====================================================================

    // derived attributes
    auto mirrorPos = Position(mirrorCenterX, 0., mirrorBackplane) + snoutFront;
    double mirrorThetaRot = std::asin(mirrorCenterX/mirrorRadius);
    double mirrorTheta1 = mirrorThetaRot - std::asin((mirrorCenterX-mirrorRmin)/mirrorRadius);
    double mirrorTheta2 = mirrorThetaRot + std::asin((mirrorRmax-mirrorCenterX)/mirrorRadius);

    // if debugging, draw full sphere
    if(mirrorDebug>0) { mirrorTheta1=0; mirrorTheta2=M_PI; mirrorPhiw=2*M_PI; };

    // solid and volume: create sphere at origin, with specified angular limits
    Sphere mirrorSolid(
        mirrorRadius-mirrorThickness,
        mirrorRadius,
        mirrorTheta1,
        mirrorTheta2,
        -mirrorPhiw/2.0,
        mirrorPhiw/2.0
        );
    Volume mirrorVol("mirror_v", mirrorSolid, mirrorMat);
    mirrorVol.setVisAttributes(mirrorVis);

    // placement (note: transformations are in reverse order)
    auto mirrorPV = vesselVol.placeVolume(mirrorVol,
          RotationZ(sectorRotation) // rotate about beam axis to sector
        * Translation3D(0,0,-mirrorRadius) // move longitudinally so it intersects snoutFront
        * Translation3D(mirrorPos.x(), mirrorPos.y(), mirrorPos.z()) // re-center to snoutFront
        * RotationY(-mirrorThetaRot) // rotate about origin
        );

    // properties
    DetElement mirrorDE(det, Form("mirror_de%d", isec), isec);
    mirrorDE.setPlacement(mirrorPV);
    SkinSurface mirrorSkin(desc, mirrorDE, Form("mirror_optical_surface%d", isec), mirrorSurf, mirrorVol);
    mirrorSkin.isValid();


    // BUILD SENSORS ====================================================================

    // if debugging sphere properties, restrict number of sensors drawn
    if(sensorSphDebug>0) { sensorSide = 2*M_PI*sensorSphRadius / 64; };

    // solid and volume: single sensor module
    Box sensorSolid(sensorSide/2., sensorSide/2., sensorThickness/2.);
    Volume sensorVol("sensor_v", sensorSolid, sensorMat);
    sensorVol.setVisAttributes(sensorVis);

    // sensitivity
    sensorVol.setSensitiveDetector(sens);

    // SENSOR MODULE LOOP ------------------------
    /* ALGORITHM: generate sphere of positions
     * - NOTE: there are two coordinate systems here:
     *   - "global" the main ATHENA coordinate system
     *   - "generator" (vars end in `Gen`) is a local coordinate system for
     *     generating points on a sphere; it is related to the global system by
     *     a rotation; we do this so the "patch" (subset of generated
     *     positions) of sensors we choose to build is near the equator, where
     *     point distribution is more uniform
     * - PROCEDURE: loop over `thetaGen`, with subloop over `phiGen`, each divided evenly
     *   - the number of points to generate depends how many sensors (+`sensorGap`)
     *     can fit within each ring of constant `thetaGen` or `phiGen`
     *   - we divide the relevant circumference by the sensor
     *     size(+`sensorGap`), and this number is allowed to be a fraction,
     *     because likely we don't care about generating a full sphere and
     *     don't mind a "seam" at the overlap point
     *   - if we pick a patch of the sphere near the equator, and not near
     *     the poles or seam, the sensor distribution will appear uniform
     */

    // initialize module number for this sector
    int imod=0;

    // thetaGen loop: iterate less than "0.5 circumference / sensor size" times
    double nTheta = M_PI*sensorSphRadius / (sensorSide+sensorGap); 
    for(int t=0; t<(int)(nTheta+0.5); t++) {
      double thetaGen = t/((double)nTheta) * M_PI;

      // phiGen loop: iterate less than "circumference at this latitude / sensor size" times
      double nPhi = 2*M_PI * sensorSphRadius * std::sin(thetaGen) / (sensorSide+sensorGap); 
      for(int p=0; p<(int)(nPhi+0.5); p++) { 
        double phiGen = p/((double)nPhi) * 2*M_PI - M_PI; // shift to [-pi,pi]

        // determine global phi and theta
        // - convert {radius,thetaGen,phiGen} -> {xGen,yGen,zGen}
        double xGen = sensorSphRadius * std::sin(thetaGen) * std::cos(phiGen);
        double yGen = sensorSphRadius * std::sin(thetaGen) * std::sin(phiGen);
        double zGen = sensorSphRadius * std::cos(thetaGen);
        // - convert {xGen,yGen,zGen} -> global {x,y,z} via rotation
        double x = zGen;
        double y = xGen;
        double z = yGen;
        // - convert global {x,y,z} -> global {phi,theta}
        double phi = std::atan2(y,x);
        double theta = std::acos(z/sensorSphRadius);

        // cut spherical patch
        // TODO [low priority]: instead of cutting a patch with complicated parameters,
        //                      can we use something simpler such as Union or
        //                      Intersection of solids?
        // - applied on global coordinates
        // - theta cuts are signed, since we will offset the patch along +x;
        //   from the offset position, toward barrel is positive, toward beam is negative
        double thetaSigned = (x<0?-1:1) * theta;
        // - position of yz planes, associated to theta cuts
        double xmin = sensorSphRadius * std::sin(sensorSphPatchThetaMin/2);
        double xmax = sensorSphRadius * std::sin(sensorSphPatchThetaMax/2);
        // - we want a phi wedge, but offset from the origin to allow more width, so
        //   define phiCheck to account for the offset; the amount of the offset,
        //   and hence the width, is controlled by `sensorSphPatchWidthFactor`
        double phiCheck = std::atan2(y,(x+sensorSphCenterX)/sensorSphPatchWidthFactor);
        // - apply cuts (only if not debugging)
        // - NOTE: use `x<xmax` for straight cuts, or `theta<sensorSphPatchThetaMax` for
        //   rounded cuts (which allows for more sensors)
        if( ( std::fabs(phiCheck)<sensorSphPatchTaper && x>xmin && theta<sensorSphPatchThetaMax && z>0 )
            || sensorSphDebug>0
        ) {

          // placement (note: transformations are in reverse order)
          // - transformations operate on global coordinates; the corresponding
          //   generator coordinates are provided in the comments
          auto sensorPV = vesselVol.placeVolume(sensorVol,
                RotationZ(sectorRotation) // rotate about beam axis to sector
              * Translation3D(sensorSphCenterX, sensorSphCenterY, sensorSphCenterZ) // move sphere to specified center
              * Translation3D(snoutFront.x(), snoutFront.y(), snoutFront.z()) // move sphere to reference position
              * RotationX(phiGen) // rotate about `zGen`
              * RotationZ(thetaGen) // rotate about `yGen`
              * Translation3D(sensorSphRadius, 0., 0.) // push radially to spherical surface
              * RotationY(M_PI/2) // rotate sensor to be compatible with generator coords
              * RotationZ(-M_PI/2) // correction for readout segmentation mapping
              );

          // generate LUT for module number -> sensor position, for readout mapping tests
          //if(isec==0) printf("%d %f %f\n",imod,sensorPV.position().x(),sensorPV.position().y());

          // properties
          sensorPV.addPhysVolID("sector", isec).addPhysVolID("module", imod);
          DetElement sensorDE(det, Form("sensor_de%d_%d", isec, imod), 10000*isec+imod);
          sensorDE.setPlacement(sensorPV);
          //SkinSurface sensorSkin(desc, sensorDE, Form("sensor_optical_surface%d", isec), sensorSurf, sensorVol);
          //sensorSkin.isValid();

          // increment sensor module number
          imod++;

        }; // end patch cuts
      }; // end phiGen loop
    }; // end thetaGen loop
    // END SENSOR MODULE LOOP ------------------------


  }; // END SECTOR LOOP //////////////////////////


  // place gas volume
  /*
  PlacedVolume gasvolPV = vesselVol.placeVolume(gasvolVol,Position(0, 0, 0));
  DetElement gasvolDE(det, "gasvol_de", 0);
  gasvolDE.setPlacement(gasvolPV);
  */

  // place mother volume (vessel)
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol,
      Position(0, 0, vesselZ0) - snoutFront
      );
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
};

// clang-format off
DECLARE_DETELEMENT(athena_DRICH, createDetector)
