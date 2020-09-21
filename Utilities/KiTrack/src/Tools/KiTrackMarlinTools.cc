#include "Tools/KiTrackMarlinTools.h"

#include <sstream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#include "UTIL/LCTrackerConf.h"

using namespace KiTrackMarlin;


std::string KiTrackMarlin::getCellID0Info( int cellID0 ){
   
   std::stringstream s;
   
   //find out layer, module, sensor
   UTIL::BitField64  cellID( UTIL::LCTrackerCellID::encoding_string() );
   cellID.setValue( cellID0 );
   
   int subdet = cellID[ UTIL::LCTrackerCellID::subdet() ] ;
   int side   = cellID[ UTIL::LCTrackerCellID::side() ];
   int module = cellID[ UTIL::LCTrackerCellID::module() ];
   int sensor = cellID[ UTIL::LCTrackerCellID::sensor() ];
   int layer  = cellID[ UTIL::LCTrackerCellID::layer() ];
   
   s << "(su" << subdet << ",si" << side << ",la" << layer << ",mo" << module << ",se" << sensor << ")";
   
   return s.str();
   
}

int KiTrackMarlin::getCellID0Layer( int cellID0 ){
   
   
   UTIL::BitField64  cellID( UTIL::LCTrackerCellID::encoding_string() );
   cellID.setValue( cellID0 );
   
   return cellID[ UTIL::LCTrackerCellID::layer() ];
   
   
}




void KiTrackMarlin::setUpRootFile( std::string rootNamePath, std::string treeName ,  std::set<std::string> branchNames , bool createNew ){
   
   
   std::string fileNamePath = rootNamePath.substr( 0 , rootNamePath.find_last_of(".")  );
   
   
   
   std::ifstream rf ((fileNamePath + ".root").c_str());       //rf for RootFile
   if ((rf) && (createNew)) { // The file already exists and we don't want to append
   
    int i=0;
    std::stringstream newname;
    while (rf){         //Try adding a number starting from 1 to the filename until no file with this name exists and use this.
      
      rf.close();
      i++;
      newname.str("");
      newname << fileNamePath << i << ".root";
      
      rf.open( newname.str().c_str() );
      
    }
    rename ( (fileNamePath + ".root").c_str() , newname.str().c_str());      //renames the file in the way,so that our new file can have it's name
    //and not ovrewrite it.
    
   }
   
   float x = 0;
   
   std::string modus = "RECREATE";
   if ( !createNew ) modus = "UPDATE";
   
   TFile* myRootFile = new TFile((fileNamePath + ".root").c_str(), modus.c_str() );        //Make new file or update it
   TTree* myTree;
   
   myTree = new TTree(treeName.c_str(),"My tree"); //make a new tree
   
   //create the branches:
   
   std::set < std::string >::iterator it;
   
   
   for ( it = branchNames.begin() ; it != branchNames.end() ; it++ ){
      
      
      myTree->Branch( (*it).c_str() , &x );
      
   }
   
   
   myTree->Write("",TObject::kOverwrite);
   myRootFile->Close();
   
}





void KiTrackMarlin::saveToRoot( std::string rootFileName, std::string treeName , std::map < std::string , float > map_name_data ){
  std::map < std::string , float >::iterator it;
     
  TFile*   myRootFile = new TFile( rootFileName.c_str(), "UPDATE"); //add values to the root file
  TTree*   myTree = dynamic_cast <TTree*>( myRootFile->Get( treeName.c_str()) );
   
  if( myTree == NULL ){
    std::cout << "ERROR: KiTrackMarlin::saveToRoot could not get tree " << treeName << std::endl;
    return;
  }
      
  for( it = map_name_data.begin() ; it != map_name_data.end() ; it++){
    myTree->SetBranchAddress( it->first.c_str(), & it->second );   
  }
     
  myTree->Fill();
  myTree->Write("",TObject::kOverwrite);
  
  myRootFile->Close();
}

void KiTrackMarlin::saveToRoot( std::string rootFileName, std::string treeName , std::vector < std::map < std::string , float > > rootDataVec ){
  std::map < std::string , float >::iterator it;
     
  TFile*   myRootFile = new TFile( rootFileName.c_str(), "UPDATE"); //add values to the root file
  TTree*   myTree = dynamic_cast <TTree*>( myRootFile->Get( treeName.c_str()) );
  
  if( myTree == NULL ){
    std::cout << "ERROR: KiTrackMarlin::saveToRoot could not get tree " << treeName << std::endl;
    return;
  }
     
  for( unsigned i=0; i<rootDataVec.size(); i++ ){ //for all entries
    std::map < std::string , float > map_name_data = rootDataVec[i];
    for( it = map_name_data.begin() ; it != map_name_data.end() ; it++){ // for all data in the entrie
      myTree->SetBranchAddress( it->first.c_str(), & it->second );   
    }
    myTree->Fill();
  }
  myTree->Write("",TObject::kOverwrite);   
  myRootFile->Close();
}


//bool KiTrackMarlin::compare_TrackerHit_z( edm4hep::ConstTrackerHit a, edm4hep::ConstTrackerHit b ){
//  return ( fabs(a.getPosition()[2]) < fabs( b.getPosition()[2]) ); //compare their z values
//}

bool KiTrackMarlin::compare_TrackerHit_z( edm4hep::ConstTrackerHit& a, edm4hep::ConstTrackerHit& b ){
  return ( fabs(a.getPosition()[2]) < fabs( b.getPosition()[2]) );
}

bool KiTrackMarlin::compare_TrackerHit_R( edm4hep::ConstTrackerHit& a, edm4hep::ConstTrackerHit& b ){
  double Rad_a2 = (a.getPosition()[0]*a.getPosition()[0]) + (a.getPosition()[1]*a.getPosition()[1]) ;
  double Rad_b2 = (b.getPosition()[0]*b.getPosition()[0]) + (b.getPosition()[1]*b.getPosition()[1]) ;
  
  return ( Rad_a2 < Rad_b2 ); //compare their radii
}



FTDHitSimple* KiTrackMarlin::createVirtualIPHit( int side , const SectorSystemFTD* sectorSystemFTD ){
   
   unsigned layer = 0;
   unsigned module = 0;
   unsigned sensor = 0;
   FTDHitSimple* virtualIPHit = new FTDHitSimple( 0.,0.,0., side , layer , module , sensor , sectorSystemFTD );
   
   
   virtualIPHit->setIsVirtual ( true );
   
   return virtualIPHit;
   
}


VXDHitSimple* KiTrackMarlin::createVirtualIPHit( const SectorSystemVXD* sectorSystemVXD ){
   
   int layer = 0 ;
   int phi = 0 ;
   int theta = 0 ;

   VXDHitSimple* virtualIPHit = new VXDHitSimple( 0.,0.,0., layer, phi, theta, sectorSystemVXD );

   virtualIPHit->setIsVirtual ( true );
   
   return virtualIPHit;
   
}


std::string KiTrackMarlin::getPositionInfo( edm4hep::ConstTrackerHit hit ){
   
   std::stringstream info;
   
   double x = hit.getPosition()[0];
   double y = hit.getPosition()[1];
   double z = hit.getPosition()[2];
   
   info << "(" << x << "," << y << "," << z << ")";
   
   return info.str();
   
}

std::string KiTrackMarlin::getPositionInfo( IHit* hit ){
   
   std::stringstream info;
   
   double x = hit->getX();
   double y = hit->getY();
   double z = hit->getZ();
   
   info << "(" << x << "," << y << "," << z << ")";
   
   return info.str();
   
}


std::string KiTrackMarlin::getTrackHitInfo( ITrack* track){
   
   std::stringstream info;
   
   std::vector< IHit* > hits = track->getHits();
   
   for( unsigned i=0; i < hits.size(); i++ ){
      
      info << getPositionInfo( hits[i] );
      
   };
   
   return info.str();
   
   
}

std::string KiTrackMarlin::getTrackHitInfo( edm4hep::Track* track){
  std::stringstream info;
  //std::vector< edm4hep::ConstTrackerHit > hits;
  unsigned int nHits = track->trackerHits_size();
  for(unsigned i=0; i<nHits; i++){
    edm4hep::ConstTrackerHit hit = track->getTrackerHits(i);
    info << getPositionInfo(hit);
  }
   
  //for( unsigned i=0; i < hits.size(); i++ ){
  //  info << getPositionInfo( hits[i] );
  //};
   
  return info.str();
}


