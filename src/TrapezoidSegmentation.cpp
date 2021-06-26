// Trapezoid Segmentation for GlueX type barrel calorimeter
// Author: M. Zurek, C. Peng (ANL)
// Date: 06/25/2021

#include "DD4hep/Factories.h"
#include "DD4hep/detail/SegmentationsInterna.h"
#include "DDSegmentation/Segmentation.h"


using namespace dd4hep::DDSegmentation;

namespace athena::seg {
  class TrapezoidGrid: public Segmentation {
    public:
      // Default constructor used by derived classes passing the encoding string
      TrapezoidGrid(const std::string& cellEncoding = "")
        : Segmentation(cellEncoding)
      {
        // define type and description
        _type = "TrapezoidGrid";
        _description = "Trapezoid segmentation in the local YZ-plane";
        registerParameters();
      }

      // Default constructor used by derived classes passing an existing decoder
      TrapezoidGrid(const BitFieldCoder* decoder)
        : Segmentation(decoder)
      {
        // define type and description
        _type = "TrapezoidGrid";
        _description = "Trapezoid segmentation in the local YZ-plane";
        registerParameters();
      }

      // Destructor
      virtual ~TrapezoidGrid() {}

      void registerParameters()
      {
        // @TODO: register all necessary parameters
        // mandatory
        registerParameter("grid_size_y", "Cell size in Y", _gridSizeY, 1., SegmentationParameter::LengthUnit);
        // optional
        registerIdentifier("identifier_y", "Cell ID identifier for Y", _yId, "y");
      }

      // override
      virtual Vector3D position(const CellID& cellID) const
      {
        // @TODO
        /* PolarGridRPhi
        Vector3D cellPosition;
        double R =   binToPosition(_decoder->get(cID,_rId),   _gridSizeR,   _offsetR);
        double phi = binToPosition(_decoder->get(cID,_phiId), _gridSizePhi, _offsetPhi);

        cellPosition.X = R * cos(phi);
        cellPosition.Y = R * sin(phi);

        return cellPosition;
        */
        return Vector3D();
      }

      // override
      virtual CellID cellID(const Vector3D& localPosition, const Vector3D& globalPosition, const VolumeID& volumeID) const
      {
        // @TODO
        /* PolarGridRPhi
        double phi = atan2(localPosition.Y,localPosition.X);
        double R = sqrt( localPosition.X * localPosition.X + localPosition.Y * localPosition.Y );
        CellID cID = vID ;
        _decoder->set(cID,_rId  , positionToBin(R, _gridSizeR, _offsetR));
        _decoder->set(cID,_phiId, positionToBin(phi, _gridSizePhi, _offsetPhi));
        */
        return volumeID;
      }

      // override
      virtual std::vector<double> cellDimensions(const CellID& cellID) const
      {
        // @TODO
        /* PolarGridRPhi
        const double rPhiSize = binToPosition(_decoder->get(cID,_rId), _gridSizeR, _offsetR)*_gridSizePhi;
        return {_gridSizeR, rPhiSize};
        */
        return std::vector<double>{0., 0., 0.};
      }


      // public methods to get member value
      double gridSizeY() const { return _gridSizeY; }
      const std::string& fieldNameY() const { return _yId; }
      // public methods to set member value
      void setGridSizeY(double cellSize) { _gridSizeY = cellSize; }
      void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

    protected:
      // @TODO: add parameters accordingly
      double _gridSizeY;
      std::string _yId;


  }; // class TrapezoidGrid
} // namespace athena::seg

DECLARE_SEGMENTATION(TrapezoidGrid, new dd4hep::SegmentationWrapper<athena::seg::TrapezoidGrid>)

