// STL
#include <vector>

// GEAR
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/BField.h>
#include <gearimpl/TPCModuleImpl.h>
#include <gearxml/GearXML.h>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

// #include <streamlog/streamlog.h>

#include "kaldet/LCTPCKalDetector.h"
#include "kaldet/ILDCylinderMeasLayer.h"

namespace kaldet
{

LCTPCKalDetector::LCTPCKalDetector(const gear::GearMgr& gearMgr)
                  : TVKalDetector(210)
{
  gear::TPCParameters const & tpcParameters = gearMgr.getTPCParameters();

  // Try to read the material properties from the GEAR file.
  // If it cannot be found use an Argon based mixture instead.
  Double_t A, Z, density, radlen;

  try{ 
    // One try/catch block for all.
    // It does not make sense to only replace part of the parameters,
    // they all have to be there.
    A = tpcParameters.getDoubleVal("TPCGas_A");
    Z = tpcParameters.getDoubleVal("TPCGas_Z");
    density = tpcParameters.getDoubleVal("TPCGas_density");
    radlen = tpcParameters.getDoubleVal("TPCGas_radlen");
  }
  catch( gear::UnknownParameterException & )
  {
    // streamlog_out(MESSAGE) << "LCTPCKalDetector: No TPCGas parameters found in the gear file."
    //            << " Using Ar/CH4 90/10." << std::endl;

    // Is this supposed to be Ar/CH4 90/10?
    // I think it's wrong. A CH4 module does not consist of 0.2 carbon
    // and 0.8 hydrogen atoms, but 1 C and 4 H, so it should be
    // A =  39.948 * 0.9 + (12.011 + 1.00794 * 4) * 0.1;
    A       = 39.948 * 0.9 + (12.011 * 0.2 + 1.00794 * 0.8) * 0.1;
    Z       = 16.4; // now how does this calculate?
    density = 0.749e-3; // in which units?
    radlen  = 1.196e4 * 2; // in which units?
  }
  
  TMaterial &gas = *new TMaterial("TPCGas", "", A, Z, density, radlen, 0.);

  // FIXME: what about the cathode tickness. Sensitivity gap in the middle?
  // And what about the LP? There it's not the half length...
  Double_t lhalf
    = tpcParameters.getMaxDriftLength();     // half length
  
  // FIXME:: needed more careful B calculation ??
  const Double_t bz = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;

//  gear::BField const & bField = gearMgr.getBField();
//  // get the BField at 0,0,0. Check that there are no transverse components
//  // FIXME: Event though there are no transverse components at 0, 0, 0 does not mean
//  // there are no transverse components somewhere else.
//  gear::Vector3D bFieldVector =  bField.at(gear::Vector3D(0., 0., 0.));
//  if (bFieldVector[0]!=0 || bFieldVector[1]!=0)
//  {
//    streamlog_out(ERROR) << "B-field has transverse components."
//           << " LCTPCKalDetector only works with homogeneous B-field"
//           << " in z direction" << std::endl;
//    throw gear::Exception("Magnetic field not homogeneous");
//  }

  // set the protected member variable of EXVKalDetector
  //  fBfield = bFieldVector[2];

  // FIXME: Just don't put dummy modules into the gear file?
  // No, we might need the radiation length
  // Damn it, the description in GEAR is not complete, grrrr.
  // Bool_t active = GearTPCMeasLayer::kActive;
  // Bool_t dummy  = GearTPCMeasLayer::kDummy;

  // The index where in the TObjectArray the next element will be stored. Unfortunately we have to
  // do the bookkeeping manually :-(
  Int_t arrayIndex = 0;

  // ILD-type Encoder
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;

  // A map to store if a layer, which is a full cylinder, already exists. If it has the same 
  // offset and the same radius, this is the case. Then to not add a Kalman layer but add the module
  // to the layer. 
  std::map< std::pair<double, double> , Int_t > uniqueLayerMap;// <offset, r>, LayerIndex

  for (std::vector<gear::TPCModule *>::const_iterator moduleIter = tpcParameters.getModules().begin();
       moduleIter < tpcParameters.getModules().end(); ++moduleIter)
  {
    gear::TPCModule *module = *moduleIter;

    // FIXME: Implementation for RectangularPadRowLayout missing
    switch (module->getLocalPadLayout().getCoordinateType())
    {
      case gear::PadRowLayout2D::POLAR :
        //Unfortunately cylinder segments are not implemented in KalDet yet.
        //Perhaps in a future implementation...
        // Get the phi min and phi max of the module, plus the right rotation
        // in the global coordinate system. This is not the global angle because it is
        // relative to the origin of the module, which has a shift in global coordinates.
        //Double_t phimin = module->getLocalModuleExtent()[2] + module->getAngle();
        //Double_t phimax = module->getLocalModuleExtent()[3] + module->getAngle();

        // the offset is either in r/phi or in x/y, depending on the global coordinate system
        double xOffset, yOffset;
        switch (module->getCoordinateType())
        {
          case gear::PadRowLayout2D::POLAR :
            xOffset=module->getOffset()[0]*cos(module->getOffset()[1]);
            yOffset=module->getOffset()[0]*sin(module->getOffset()[1]);
            break;
          case  gear::PadRowLayout2D::CARTESIAN :
            xOffset=module->getOffset()[0];
            yOffset=module->getOffset()[1];
            break;
          default:
            throw gear::UnknownParameterException("Unknown global coordinate system in LCTPCKalDetector");

        }

        if (yOffset!=0)
        {
          // streamlog_out(ERROR) << "yOffset not 0. The current implementation of LCTPCKalDetector"
          //              << "only allows an offset in x direction." << std::endl;
          throw gear::UnknownParameterException("Offset in y is not 0.");
        }

        // Loop the pad rows and place one cylinder segment each
        for (int row = 0; row < module->getNRows(); ++row)
        {
          // To get the radius we have to ask for the center of one one the pads,
          // just take pad 0 in this row.
          int padIndex  = module->getPadIndex(row, 0);

          // The radius only makes sense in local module coordinates.
          // So we have to ask the local pad layout. The module always answers in global coordinates.
          Double_t r  = module->getLocalPadLayout().getPadCenter(padIndex)[0];

          std::map< std::pair<double, double> , Int_t >::iterator uniqueLayerIterator=
            uniqueLayerMap.find( std::pair<double, double>(xOffset, r) );

          if ( uniqueLayerIterator==uniqueLayerMap.end() )
          {
            // streamlog_out(DEBUG2) << "adding new layer "<<arrayIndex<< " (xOffset="<<xOffset
            //           <<", r="<<r
            //           <<") with module "<< module->getModuleID()
            //           <<" and row "<<row << std::endl;
            // add the measurement layer to this object, and remember where we put it
            encoder.reset() ;  // reset to 0

            encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::TPC ;

            int mod = module -> getModuleID() ;
            int row_global =
        	// Modules 0 and 1 get the same row index
        	( mod==0 || mod==1 ) ? row :
		/* For modules 2, 3 and 4 the row gets increased
		 * by the number of rows in modules 0/1
		 */
              ( ( mod==2 || mod==3 || mod==4) ? tpcParameters.getModule(0).getNRows() + row :
  		/* For modules 5 and 6 the row gets increased
  		 * by the number of rows in modules 0/1 and 2/3/4
  		 */
        	tpcParameters.getModule(0).getNRows() + tpcParameters.getModule(2).getNRows() + row );

            encoder[lcio::ILDCellID0::layer] = row_global ;

	    //	    encoder[lcio::ILDCellID0::module] = mod ;

            int CellID = encoder.lowWord() ;

            Add(new ILDCylinderMeasLayer(gas, gas,
                    r, lhalf,
                    xOffset, yOffset, 0, // FIXME: cathode thickness?
                    bz,
                    true, // active
                    CellID,
                    "ILDCylinderMeasLayer"));

            // add the layer to the uniqueLayerMap
            uniqueLayerMap[std::pair<double, double>(xOffset, r)] = arrayIndex;

            // count up the array index
            ++arrayIndex;
          }
          else // layer already exists
          {
            // streamlog_out(DEBUG2)
            //         << "A new layer will not be added as the layer with the same parameters already exist"<< std::endl;
          }

        }// loop rows in module
        break;

      case  gear::PadRowLayout2D::CARTESIAN :
        throw gear::NotImplementedException("Cartesian local coordinates not yet supported by LCTPCKalDetector");

      default:
        throw gear::UnknownParameterException("Unknown local coordinate system in LCTPCKalDetector");
    }//switch coordinate type
    
  }

  SetOwner();// make the KalDetector (TObjectArry) owner of the measurement layers so they are
             // deleted when the detector is deleted

}

LCTPCKalDetector::~LCTPCKalDetector()
{
}

}//namespace kaldet
