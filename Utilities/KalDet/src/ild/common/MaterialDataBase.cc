
#include "MaterialDataBase.h"

#include <stdexcept>
#include <vector>
#include <algorithm>

#include <gear/GEAR.h>
#include "gearimpl/Util.h"
#include <gear/SimpleMaterial.h>
#include <gearimpl/GearParametersImpl.h>
#include "DetInterface/IGeomSvc.h"

#include "TMaterial.h"

MaterialDataBase::~MaterialDataBase(){
  
  std::map<std::string,TMaterial* >::iterator it = _material_map.begin();
  std::vector<TMaterial*> deleted_objects;
  
  for( /**/; it!=_material_map.end(); ++it) 
    
    if( std::find( deleted_objects.begin(), deleted_objects.end(), (*it).second ) != deleted_objects.end() ) {
      delete (*it).second ;
      deleted_objects.push_back((*it).second) ;
    }
}

TMaterial* MaterialDataBase::getMaterial(std::string mat_name){
  
  if (_isInitialised == false) {    
    MaterialDataBaseException exp ( mat_name ) ;
    
    
    throw exp;
  }
  
  std::map<std::string,TMaterial* >::iterator it = _material_map.find(mat_name) ;        
  
  if ( it == _material_map.end() ) { 
    MaterialDataBaseException exp( mat_name ) ;
    throw exp ; 
  } 
  else { 
    return (*it).second ; 
  }
  
}

void MaterialDataBase::initialise( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc ){
  
  if( !_isInitialised ){
    this->createMaterials(gearMgr, geoSvc); 
    _isInitialised = true ;
    _gearMgr = &gearMgr;
  }
  
}

void MaterialDataBase::registerForService(const gear::GearMgr& gearMgr, IGeomSvc* geoSvc ) {
  
  if( !_isInitialised ){
    //std::cout << "debug fucd: " << "--------------------" << std::endl;
    this->initialise(gearMgr, geoSvc);
  }
  
  else {
    if ( _gearMgr != &gearMgr ) {
      MaterialDataBaseException exp( " MaterialDataBase::registerForService : _gearMgr != &gearMgr !  " ) ;
      throw exp;
    }
  }
  
}


void MaterialDataBase::addMaterial(TMaterial* mat, std::string name) {
  std::map<std::string, TMaterial*>::iterator it = _material_map.find(name) ; 
  //std::cout << name << " " << mat << std::endl;
  std::string what( name ) ;
  what += " - already exists [MaterialDataBase::addMaterial() ] " ;
  MaterialDataBaseException exp( what ) ;

  if ( it != _material_map.end() ) { 
    //std::cout << what << std::endl;
    throw exp; 
  } 
  else { 
    _material_map[name] = mat  ; 
  }
}


void MaterialDataBase::createMaterials(const gear::GearMgr& gearMgr, IGeomSvc* geoSvc ){
  
  Double_t A, Z, density, radlen ;
  std::string name;
  
  // Beam
  A       = 14.00674 * 0.7 + 15.9994 * 0.3 ;
  Z       = 7.3 ;
  density = 1.0e-25 ; // density set to very low value
  radlen  = 1.0e25 ;  // give huge radiation length 
  name    = "beam" ;
  
  TMaterial& beam = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&beam, name);
  
  // Air
  A       = 14.00674 * 0.7 + 15.9994 * 0.3 ;
  Z       = 7.3 ;
  density = 1.205e-3 ; // g/cm^3
  radlen  = 3.42e4 ;   // cm // SJA:FIXME using the standard formular for the radiation length and the above values gives 3.06e4
  name    = "air" ;
  
  TMaterial &air = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&air, name) ;
  
  // Si
  A       = 28.09 ;
  Z       = 14.0 ;
  density = 2.33 ; // g/cm^3
  radlen  = 9.36607 ;   // cm 
  name    = "silicon";
  
  TMaterial &silicon = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&silicon, name);
  
  // C
  A       = 12.01 ;
  Z       = 6.0 ;
  density = 2.00 ; // g/cm^3
  radlen  = 21.3485 ;   // cm 
  name    = "carbon" ;
  
  TMaterial &carbon = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&carbon, name);
  
  // Aluminium
  A       = 26.9815 ;
  Z       = 13.0 ;
  density = 2.699 ; // g/cm^3
  radlen  = 8.9 ;   // cm 
  name    = "aluminium" ;
  
  TMaterial &aluminium = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&aluminium, name);

  // Beryllium
  A       = 9.012 ;
  Z       = 4.0 ;
  density = 1.85 ; // g/cm^3
  radlen  = 35.28 ;   // cm 
  name    = "beryllium" ;
  
  TMaterial &beryllium = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&beryllium, name);

  
  // TPC Gas
  A       = 39.948*0.9+(12.011*0.2+1.00794*0.8)*0.1;
  Z       = 16.4;
  density = 1.749e-3 ;
//  radlen  =  1.196e4*2;
  radlen  =  0.5*1.196e4*2; // SJA:FIXME: reduce by a factor of 10%
  name    = "tpcgas" ;
  
  TMaterial &tpcgas = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.);
  this->addMaterial(&tpcgas, name);
  
  // TPC Field Cage
  A       = air.GetA()*0.97 + aluminium.GetA()*0.03 ; // SJA:FIXME just use this simple approximation for now
  Z       = air.GetZ()*0.97 + aluminium.GetZ()*0.03 ; // SJA:FIXME just use this simple approximation for now 
  density = 0.0738148 ;
  radlen  =  489.736 * 0.5 ; // SJA:FIXME just use factor of two for now as the amount differs by ~ factor of 2 from observation in GEANT4 
  name    = "tpcinnerfieldcage" ;
  
  TMaterial &tpcinnerfieldcage = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.);
  this->addMaterial(&tpcinnerfieldcage, name);
  
  // TPC Field Cage
  A       = air.GetA()*0.97 + aluminium.GetA()*0.03 ; // SJA:FIXME just use this simple approximation for now
  Z       = air.GetZ()*0.97 + aluminium.GetZ()*0.03 ; // SJA:FIXME just use this simple approximation for now
  density = 0.0738148 ;
  radlen  =  489.736 ; 
  name    = "tpcouterfieldcage" ;
  
  TMaterial &tpcouterfieldcage = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.);
  this->addMaterial(&tpcouterfieldcage, name);
  
  
  // FTD Support Material
  // Needed because of the lack of space frame in the Tracking Geometry description
  // Carbon with density and rad-length changed by a factor of 2
  
  A       = 12.01 ;
  Z       = 6.0 ;
  density = 0.5 * 2.00 ; // g/cm^3
  radlen  = 2.0 * 21.3485 ;   // cm
  name    = "FTDSupportMaterial" ;
  
  TMaterial &ftdsupport = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&ftdsupport, name);

  
  // VXD Support Material
  if(geoSvc){
    TMaterial* vxdsupport = geoSvc->getMaterial("VXDSupportMaterial");
    if(vxdsupport) this->addMaterial(vxdsupport, "VXDSupportMaterial");
    else std::cout << "Material VXDSupportMaterial not found" << std::endl;
  }
  else{
    try{
      const gear::SimpleMaterial& vxd_sup_mat = gearMgr.getSimpleMaterial("VXDSupportMaterial");
      A       = vxd_sup_mat.getA();
      Z       = vxd_sup_mat.getZ();
      density = vxd_sup_mat.getDensity() * (1000.0/ 1000000.0); // kg/m^3 -> g/cm^3
      radlen  = vxd_sup_mat.getRadLength() / 10.0 ; // mm -> cm
      name    = vxd_sup_mat.getName() ;
      TMaterial &vxdsupport = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.);
      this->addMaterial(&vxdsupport, name);
    }
    catch( gear::UnknownParameterException& e){   
      std::cout << "Error while read material from GeomSvc!" << std::endl;
    }
  }
}

