#ifndef MaterialDataBase_h
#define MaterialDataBase_h

/** MaterialDataBase: Class to hold and manage collection of materials  
 *
 * @author S.Aplin DESY
 */

#include <string>
#include <map>
#include <exception>

#include "lcio.h"
#include "Exceptions.h"

class TMaterial;

namespace gear{
  class GearMgr ;
}
class IGeomSvc;
// fg: define the MaterialDataBaseException as an lcio Exception to allow for 
//     messages to be printed in what() 
typedef lcio::Exception MaterialDataBaseException ;

class MaterialDataBase {
  
public:
  
  /** Accessor Method */
  static MaterialDataBase& Instance() {
    
    static MaterialDataBase singleton;
        
    return singleton;
    
  }
  
  // Other non-static member functions
  
public:
  
  /** Destructor */
  ~MaterialDataBase();   
  
  /** Get Material via name */
  TMaterial* getMaterial(std::string mat_name) ;  
  
  void registerForService(const gear::GearMgr& gearMgr, IGeomSvc* geoSvc=0) ;
  
  
private:
 
  void initialise(const gear::GearMgr& gearMgr, IGeomSvc* geoSvc) ;
  
  MaterialDataBase() { _material_map.clear(); _isInitialised = false ; _gearMgr = 0; }                               // Private constructor
  
  
  MaterialDataBase(const MaterialDataBase&) ;                 // Prevent copy-construction
  MaterialDataBase& operator=(const MaterialDataBase&) ;      // Prevent assignment
  
  void addMaterial(TMaterial* mat, std::string name); 
  void createMaterials(const gear::GearMgr& gearMgr, IGeomSvc* geoSvc);
  
  // private member variables
  std::map<std::string,TMaterial* > _material_map;
  
  bool _isInitialised;
  
  const gear::GearMgr* _gearMgr;
  
};



#endif
