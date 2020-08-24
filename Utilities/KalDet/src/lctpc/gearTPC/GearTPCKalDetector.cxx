// STL
#include <vector>

// GEAR
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/BField.h>
#include <gearimpl/TPCModuleImpl.h>
#include <gearxml/GearXML.h>

// #include <streamlog/streamlog.h>

#include "GearTPCKalDetector.h"
#include "GearTPCCylinderMeasLayer.h"
//#include "GearTPCHit.h"
#include <kaltest/TKalDetCradle.h>
#include <TRandom.h>
#include <TMath.h>

//#include "TTUBE.h"
//#include "TNode.h"
//#include "TVirtualPad.h"

namespace kaldet
{

GearTPCKalDetector::GearTPCKalDetector(const gear::GearMgr& gearMgr)
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
    // streamlog_out(MESSAGE) << "GearTPCKalDetector: No TPCGas parameters found in the gear file."
    //     		   << " Using Ar/CH4 90/10." << std::endl;

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
  
//  gear::BField const & bField = gearMgr.getBField();
//  // get the BField at 0,0,0. Check that there are no transverse components
//  // FIXME: Event though there are no transverse components at 0, 0, 0 does not mean
//  // there are no transverse components somewhere else.
//  gear::Vector3D bFieldVector =  bField.at(gear::Vector3D(0., 0., 0.));
//  if (bFieldVector[0]!=0 || bFieldVector[1]!=0)
//  {
//    streamlog_out(ERROR) << "B-field has transverse components."
//			 << " GearTPCKalDetector only works with homogeneous B-field"
//			 << " in z direction" << std::endl;
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

  // A map to store if a layer, which is a full cylinder, already exists. If it has the same 
  // offset and the same radius, this is the case. Then to not add a Kalman layer but add the module
  // to the layer. 
  std::map< std::pair<double, double> , Int_t > uniqueLayerMap;// <offset, r>, LayerIndex

  for (std::vector<gear::TPCModule *>::const_iterator moduleIter = tpcParameters.getModules().begin();
       moduleIter < tpcParameters.getModules().end(); ++moduleIter)
  {
    gear::TPCModule *module = *moduleIter;

    // The resolution parameters can vary from module to module.
    // In case they are not there use the jgem values
    Double_t sigmax0, sigmax1, sigmaz0, sigmaz1;
    try{
      sigmax0 = module->getDoubleVal("sigmax0");
      sigmax1 = module->getDoubleVal("sigmax1");
      sigmaz0 = module->getDoubleVal("sigmaz0");
      sigmaz1 = module->getDoubleVal("sigmaz1");
    }
    catch( gear::UnknownParameterException & )
    {
      // streamlog_out(MESSAGE) << "No resolution parameters found for module "
      //   		     << module->getModuleID()<<"."
      //   		     << " Using jgem settings." << std::endl;

      // FIXME: n_eff of argon per distance, adapt to pad length. Or read from GEAR file?
      Double_t neff    = 22.7;

      // FIXME : I think the sqrt(10) is a leftover from the old cm units....
      sigmax0 = 38.3e-3;
      sigmax1 = 101.5e-3 / TMath::Sqrt(10.) / TMath::Sqrt(neff);
      sigmaz0 = 500.e-3 ;
      sigmaz1 = 154.e-3  / TMath::Sqrt(10.) / TMath::Sqrt(neff);
    }

    // FIXME: Implementation for RectangularPadRowLayout missing
    switch (module->getLocalPadLayout().getCoordinateType())
    {
      case gear::PadRowLayout2D::POLAR :
      {
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
	    throw gear::UnknownParameterException("Unknown global coordinate system in GearTPCKalDet");
	    
	}

	if (yOffset!=0)
	{
	  // streamlog_out(ERROR) << "yOffset not 0. The current implementation of GearTPCKalDetector"
	  //       	       << "only allows an offset in x direction." << std::endl;
	  throw gear::UnknownParameterException("Offset in y is not 0.");
	}
	
	TVector3 cylinderCentre(xOffset,
				yOffset,
				0.); // FIXME: cathode thickness

	// Loop the pad rows and place one cylinder segement each
	for (int row = 0; row < module->getNRows(); ++row)
	{
	  // To get the radius we have to ask for the centre of one one the pads,
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
	    //     		  <<", r="<<r
	    //     		  <<") with module "<< module->getModuleID()
	    //     		  <<" and row "<<row << std::endl;
	    // add the measurement layer to this object, and remember where we put it
	    Add(new GearTPCCylinderMeasLayer(gas, gas, 
					     module->getModuleID(), row,
					     r, lhalf, 
					     cylinderCentre,
					     true, true, // perfect and active
					     sigmax0, sigmax1, sigmaz0, sigmaz1));
	    // add the moduleRow to the map with the array indices
	    moduleRowToMeasurementLayerMap[ std::pair<int, int>(module->getModuleID(), row) ]
	      = arrayIndex;

	    // add the layer to the uniqueLayerMap
	    uniqueLayerMap[std::pair<double, double>(xOffset, r)] = arrayIndex;

	    // count up the array index
	    ++arrayIndex;
	  }
	  else // layer already exists
	  {
	    Int_t existingIndex = uniqueLayerIterator->second;
	    // streamlog_out(DEBUG2) << "adding module "<< module->getModuleID()
	    //     		  <<", row "<<row <<" to existing layer " << existingIndex
	    //     		  << " (xOffset="<<xOffset<<", row="<<row<<")"
	    //     		  << std::endl;
	    
	    GearTPCCylinderMeasLayer * measLayer = 
	      dynamic_cast<GearTPCCylinderMeasLayer *>( At(existingIndex));
	    
	    // If the cast fails something went terribly wrong. This layer should be a 
	    // GearTPCCylinderMeasLayers, otherwise we cannot add a cylindrical module
	    if (measLayer==0)
	    {
	      // streamlog_out(ERROR) << "Something went terribly wrong. Could not cast a member of "
	      //   		   << " GerTPCKalDetector to GearTPCCylinderMeasLayer" << std::endl;
	      throw std::bad_cast();
	    }

	    measLayer->AddModuleRow( module->getModuleID(), row );

	    moduleRowToMeasurementLayerMap[ std::pair<int, int>(module->getModuleID(), row) ]
	      = existingIndex;
	    
	  }
	    
	}// loop rows in module
      }
      break;
      case  gear::PadRowLayout2D::CARTESIAN :
	throw gear::NotImplementedException("Cartesian local coordinates not yet supported by GearTPCKalDetector");
      default:
	throw gear::UnknownParameterException("Unknown local coordinate system in GearTPCKalDet");
      }//switch coordinate type
    
  }

  SetOwner();// make the KalDetector (TObjectArry) owner of the measurement layers so they are
             // deleted when the detector is deleted

}

GearTPCKalDetector::~GearTPCKalDetector()
{
}

 
GearTPCMeasLayer const * GearTPCKalDetector::GetMeasLayer(int moduleID, int row) const
{
  std::map< std::pair<int, int >, Int_t >::const_iterator indexIter = 
    moduleRowToMeasurementLayerMap.find( std::pair<int, int>(moduleID, row) );

  if (indexIter==moduleRowToMeasurementLayerMap.end())
  {
    // streamlog_out(ERROR) << "GearTPCKalDetector::getMeasLayer: "
    //     		 << "layer with moduleID=" << moduleID 
    //     		 << " and row=" << row << " not defined." << std::endl;
    throw gear::Exception("Unknown row on moduleID");
  }
  
  // second is the key of the iterator (the arrow dereferencing the iterator to a key-value pair)
  GearTPCMeasLayer const * measurementLayer = 
    dynamic_cast<GearTPCMeasLayer const * >( At(indexIter->second) );
  if (measurementLayer==0)
  {
    // streamlog_out(ERROR) << "GearTPCKalDetector::getMeasLayer: "
    //     		 << " cannot cast object this->At(" << indexIter->second
    //     		 << ") to GearTPCMeasLayer const * " << std::endl;
  }
  return measurementLayer;
}

//_________________________________________________________________________
//  --------------
//  Utility Method
//  --------------
//_________________________________________________________________________
//  --------------------------------------
//  Draw: Drawing method for event display
//  --------------------------------------
//
//void GearTPCKalDetector::Draw(Int_t color, const Char_t *opt)
//{
//  if (! gPad) return;
//  TNode *nodep = GetNodePtr();
//  nodep->cd();
//
//  if (! fNodePtr) {
//    GearTPCMeasLayer *inp  = static_cast<GearTPCMeasLayer *>(First());
//    GearTPCMeasLayer *outp = static_cast<GearTPCMeasLayer *>(Last());
//    Double_t rin  = inp->GetR();
//    Double_t rout = outp->GetR();
//    Double_t hlen = outp->GetZmax();
//    const Char_t *name  = "TPC";
//    const Char_t *nname = "TPCNode";
//    TTUBE *tubep = new TTUBE(name, name, "void", rin, rout, hlen);
//    tubep->SetBit(kCanDelete);
//    fNodePtr = new TNode(nname, nname, name);
//    fNodePtr->SetLineColor(color);
//    fNodePtr->SetLineWidth(1); // line width given in number of pixels
//  }
//  EXVKalDetector::Draw(color, opt);
//}

}//namespace kaldet
