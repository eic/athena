R__LOAD_LIBRARY(libDDCore.so)
// R__LOAD_LIBRARY(libActsPluginDD4hep.so)
R__LOAD_LIBRARY(libDDG4.so)
R__LOAD_LIBRARY(libDDG4IO.so)
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Detector.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TGeoMedium.h"
#include "TGeoManager.h"
#include "DDRec/MaterialScan.h"
#include "DDRec/MaterialManager.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Printout.h"
#include "fmt/core.h"

#include <iostream>
#include <fstream>

// #include "Acts/Geometry/TrackingGeometry.hpp"
// #include "Acts/Geometry/TrackingVolume.hpp"
// #include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
using namespace dd4hep;
using namespace dd4hep::rec;

void test_matscan(const char* compact = "athena.xml", TString face="z"){

  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(compact);
  MaterialScan matscan(detector); 

  fmt::print("\n");
  fmt::print("All detector subsystem names:\n");
  for(const auto&  d : detector.detectors() ) {
    fmt::print("  {}\n", d.first);
  }
  return;

  TString det_list[14]={"TrackerBarrel_Inner","TrackerBarrel_Outer","TrackerEndcapN_Inner","TrackerEndcapN_Outer","TrackerEndcapP_Inner","TrackerEndcapP_Outer","TrackerSubAssembly_Inner","TrackerSubAssembly_Outer","VertexBarrel","VertexBarrelSubAssembly","VertexEndcapN","VertexEndcapP","VertexEndcapSubAssembly","cb_DIRC"};
  for (int dd=0;dd<14;dd++){
  TString detname = det_list[dd];
  matscan.setDetector(detname.Data());      

  double x0=0,y0=0,z0=0,x1,y1,z1;  // cm
  double epsilon=1e-4; // (mm) default 1e-4: Materials with a thickness smaller than epsilon (default 1e-4=1mu
  const char* fmt1 = "%7.2f %7.2f %7.2f %5d %-20s %3.0f %8.3f %8.4f %11.4f %11.4f %10.3f %8.2f %11.6f %11.6f %7.2f %7.2f %7.2f\n";
  const char* fmt2 = "%7.2f %7.2f %7.2f %5d %-20s %3.0f %8.3f %8.4f %11.6g %11.6g %10.3f %8.2f %11.6f %11.6f %7.2f %7.2f %7.2f\n";


  // x1 = 100; y1 = 100; z1 = 100; 
  // y1 = 100;
  double a1,a2,a3;
  for(a3=-100;a3<101;a3=a3+200){
    TString fname = Form("/global/u2/s/shujie/eic/output/matscan/%s_%s%g.dat",detname.Data(),face.Data(),a3);
    FILE * pFile;
    pFile = fopen (fname,"w");

  for(a1=-100;a1<100;a1=a1+1){
  for(a2=-100;a2<100;a2=a2+1){
    if (face=="x"){
      y1=a1; z1=a2; x1=a3;
    }
    else if (face=="y"){
      x1=a1; z1=a2; y1=a3;
    }
    else if (face=="z"){
      x1=a1; y1=a2; z1=a3;
    }

    Vector3D p0(x0, y0, z0), p1(x1, y1, z1);
    Vector3D end, direction;
    direction = (p1-p0).unit();
  	const auto& placements = matscan.scan(x0,y0,z0,x1,y1,z1,epsilon); 
  	// matscan.print(x0,y0,z0,x1,y1,z1,epsilon);
    // return;
    // Vector3D end, direction;
    double sum_x0 = 0;
    double sum_lambda = 0;
    double path_length = 0, total_length = 0;

    // TString fname = Form("/global/u2/s/shujie/eic_dir/matscan_%g_%g_%g.dat",x1,y1,z1);

  	for (unsigned i=0;i<placements.size();i++){

      TGeoMaterial* mat=placements[i].first->GetMaterial();
      double length = placements[i].second;
      double nx0     = length / mat->GetRadLen();
      double nLambda = length / mat->GetIntLen();
      sum_x0        += nx0;
      sum_lambda    += nLambda;
      path_length   += length;
      total_length  += length;
      end = p0 + total_length * direction;
      const char* fmt = mat->GetRadLen() >= 1e5 ? fmt2 : fmt1;

      // fprintf(pFile, "%d\n",i+1);
      fprintf(pFile, fmt,x1, y1, z1, i+1, mat->GetName(), mat->GetZ(), mat->GetA(),
                mat->GetDensity(), mat->GetRadLen(), mat->GetIntLen(), 
                length, path_length, sum_x0, sum_lambda, end[0], end[1], end[2]);

      // ::printf(fmt, i+1, mat->GetName(), mat->GetZ(), mat->GetA(),
      //           mat->GetDensity(), mat->GetRadLen(), mat->GetIntLen(), 
      //           length, path_length, sum_x0, sum_lambda, end[0], end[1], end[2]);
      //  // mat->Print();

  	} 
  }
    cout<<detname<<"  "<<x1<<","<<y1<<","<<z1<<endl;
    // cout<<x1<<","<<y1<<","<<z1<<": "<<placements.size()<<"  "<<sum_x0<<"  "<<total_length<<endl;
  }	
  fclose (pFile);

  }
}
}

/*
  <!-- <include ref="ecal_barrel_hybrid.xml"/> -->
<TGeoManager::CountLevels>: max level = 5, max placements = 2796
Error in <TGeoVoxelFinder::SortAll>: Volume B0Tracker: Cannot make slices on any axis
Error in <TGeoVoxelFinder::SortAll>: Volume ForwardRomanPot_Station_1: Cannot make slices on any axis


{  B0APF_BeamlineMagnet
B0PF_BeamlineMagnet
B0Tracker
B1APF_BeamlineMagnet
B1PF_BeamlineMagnet
B2PF_BeamlineMagnet
BPFR1_BeamlineMagnet
BackwardTOF
BarrelTOF
BeamPipe
EcalEndcapN
EcalEndcapP
ForwardOffMTracker
ForwardRomanPot_Station_1
ForwardTOF
ForwardTRD
GaseousRICH
HcalBarrel
HcalEndcapN
HcalEndcapP
Q1APF_BeamlineMagnet
Q1BPF_BeamlineMagnet
Q2PF_BeamlineMagnet
QPFC1_BeamlineMagnet
QPFC2_BeamlineMagnet
QPFC3_BeamlineMagnet
QPFC4_BeamlineMagnet
QPFR1_BeamlineMagnet
QPFR2_BeamlineMagnet
SolenoidCoilBarrel
SolenoidCoilEndcapN
SolenoidCoilEndcapP
TOFSubAssembly
TrackerBarrel_Inner
TrackerBarrel_Outer
TrackerEndcapN_Inner
TrackerEndcapN_Outer
TrackerEndcapP_Inner
TrackerEndcapP_Outer
TrackerSubAssembly_Inner
TrackerSubAssembly_Outer
VertexBarrel
VertexBarrelSubAssembly
VertexEndcapN
VertexEndcapP
VertexEndcapSubAssembly
cb_DIRC
ce_MRICH
ffi_ZDC_ECAL
ffi_ZDC_HCAL}


*/
