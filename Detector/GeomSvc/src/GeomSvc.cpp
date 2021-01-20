#include "GeomSvc.h"
#include "gearimpl/GearParametersImpl.h"
#include "TMath.h"
#include "TMaterial.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "DD4hep/Detector.h"
#include "DD4hep/Plugins.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"
#include "DDRec/DetectorData.h"

#include <iomanip>
#include <iostream>

DECLARE_COMPONENT(GeomSvc)

GeomSvc::GeomSvc(const std::string& name, ISvcLocator* svc)
: base_class(name, svc), m_dd4hep_geo(nullptr), m_vxdData(nullptr), m_beamPipeData(nullptr){

}

GeomSvc::~GeomSvc() {

}

StatusCode
GeomSvc::initialize() {
  StatusCode sc = Service::initialize();

  m_dd4hep_geo = &(dd4hep::Detector::getInstance());
  // if failed to load the compact, a runtime error will be thrown.
  m_dd4hep_geo->fromCompact(m_dd4hep_xmls.value());
  
  dd4hep::DetElement world = m_dd4hep_geo->world();
  //info() << world.type() << " " << world.path() << " " << world.placementPath() << endmsg;
  const std::map<std::string, dd4hep::DetElement>& subs = world.children();
  for(std::map<std::string, dd4hep::DetElement>::const_iterator it=subs.begin();it!=subs.end();it++){
    dd4hep::DetElement sub = it->second;
    info() << it->first << " " << sub.path() << " " << sub.placementPath() << endmsg;
    if(it->first=="Tube"){
      dd4hep::Volume vol = sub.volume();
      dd4hep::Solid  solid = vol.solid();
      //info() << " " << solid.type() << " " << solid.name() << endmsg;
      const std::map<std::string, dd4hep::DetElement>& pipes = sub.children();
      for(std::map<std::string, dd4hep::DetElement>::const_iterator it=pipes.begin();it!=pipes.end();it++){
	dd4hep::DetElement pipe = it->second;
	//info() << " " << it->first << " " << pipe.id() << " " << pipe.path() << " " << pipe.placementPath() << endmsg;
      }
      try{
	m_beamPipeData = sub.extension<dd4hep::rec::ConicalSupportData>();
      }
      catch(std::runtime_error& e){
	info() << e.what() << " " << m_beamPipeData << endmsg;
      }
    }
    if(it->first=="VXD"){
      sc = convertVXD(sub);
    }
  }
  return sc;
}

StatusCode
GeomSvc::finalize() {
  StatusCode sc;
  if(m_vxdParameters) delete m_vxdParameters;
  return sc;
}

dd4hep::DetElement
GeomSvc::getDD4HepGeo() {
    if (lcdd()) {
        return lcdd()->world();
    }
    return dd4hep::DetElement();
}

dd4hep::Detector*
GeomSvc::lcdd() {
    return m_dd4hep_geo;
}


IGeomSvc::Decoder*
GeomSvc::getDecoder(const std::string& readout_name) {

    IGeomSvc::Decoder* decoder = nullptr;

    if (!lcdd()) {
        error() << "Failed to get lcdd()" << endmsg;
        return decoder;
    }

    auto readouts = m_dd4hep_geo->readouts();
    if (readouts.find(readout_name) == readouts.end()) {
        error() << "Failed to find readout name '" << readout_name << "'"
                << " in DD4hep::readouts. "
                << endmsg;
        return decoder;
    }
    
    dd4hep::Readout readout = lcdd()->readout(readout_name);
    auto m_idspec = readout.idSpec(); 

    decoder = m_idspec.decoder();

    if (!decoder) {
        error() << "Failed to get the decoder with readout '"
                << readout_name << "'" << endmsg;
    }

    return decoder;

}


const std::map<std::string,double>& GeomSvc::getDetParameters(std::string name){
  if(m_detParameters.find(name)!=m_detParameters.end()) return m_detParameters[name];
  else{
    char message[200];
    sprintf(message,"GeomSvc has not the parameter set named %s", name); 
    throw std::runtime_error(message);
  }
}

double GeomSvc::getDetParameter(std::string set_name, std::string par_name){
  std::map<std::string, std::map<std::string,double> >::iterator it=m_detParameters.find(set_name);
  if(it!=m_detParameters.end()){
    if(it->second.find(par_name)!=it->second.end()) return it->second[par_name];  
  }
  char message[200];
  sprintf(message,"GeomSvc has not the parameter named %s in set %s", par_name, set_name);
  throw std::runtime_error(message);
}

StatusCode GeomSvc::convertVXD(dd4hep::DetElement& vxd){
  StatusCode sc;
  //fucd: another method to obtain parameters, but not fully for KalDet  
  bool extensionDataValid = true;
  try{
    m_vxdData = vxd.extension<dd4hep::rec::ZPlanarData>();
  }
  catch(std::runtime_error& e){
    extensionDataValid = false;
    info() << e.what() << " " << m_vxdData << endmsg;
  }
  
  std::vector<helpLayer> helpSensitives;
  std::vector<helpLayer> helpLadders;
  std::vector<int>       helpNumberLadders;
  std::vector<double>    helpPhi0; 
  int                    helpCount=0;
  int                    type=0;
  double shellInnerRadius, shellOuterRadius, shellHalfLength, gap, shellRadLength;
  int    nLadders=0;
  double phi0=0;
  helpLayer thisLadder;
  double beryllium_ladder_block_length=0,end_electronics_half_z=0,side_band_electronics_width=0;
  double rAlu, drAlu, rSty, drSty, dzSty, rInner, aluEndcapZ, aluHalfZ, alu_RadLen, Cryostat_dEdx;
  double VXDSupportDensity, VXDSupportZeff, VXDSupportAeff, VXDSupportRadLen;
  dd4hep::Volume vxd_vol = vxd.volume();
  for(int i=0;i<vxd_vol->GetNdaughters();i++){
    TGeoNode* daughter = vxd_vol->GetNode(i);
    std::string nodeName = daughter->GetName();
    //info() << daughter->GetName() << endmsg; 
    if(nodeName=="VXD_support_assembly_0"){
      TGeoNode* shell = FindNode(daughter, "SupportShell");
      if(shell){
	const TGeoShape* shape_shell = shell->GetVolume()->GetShape();
	//fucd: IsA() method does not always work for TGeoTube, sometimes, strange?
	//if(shape_shell->IsA()==TGeoTube::Class()){
	if(shape_shell->TestShapeBit(TGeoTube::kGeoTube)){
	  const TGeoTube* tube = (const TGeoTube*) shape_shell;
	  shellInnerRadius = tube->GetRmin()*CLHEP::cm;
	  shellOuterRadius = tube->GetRmax()*CLHEP::cm;
	  shellHalfLength  = tube->GetDz()*CLHEP::cm;
	}
	else{
	  error() << shell->GetName() << " is not a TGeoTube!  Shape bits = " << shape_shell->TestShapeBits(0xFFFFFFFF) << endmsg; 
	}
	TGeoMaterial* mat = shell->GetMedium()->GetMaterial();
	shellRadLength = mat->GetRadLen()*CLHEP::cm;
      }
      TGeoNode* block = FindNode(daughter, "BerylliumAnnulusBlock");
      if(block){
	const TGeoShape* shape_block = block->GetVolume()->GetShape();
	if(shape_block->IsA()==TGeoBBox::Class()){
	  const TGeoBBox* box = (const TGeoBBox*) shape_block;
	  beryllium_ladder_block_length = box->GetDY()*CLHEP::cm;
	}
	else{
	  error() << block->GetName() << " is not a TGeoTube!  Shape bits = " << shape_block->TestShapeBits(0xFFFFFFFF) << endmsg;
	}
      }
      TGeoNode* skin = FindNode(daughter, "CryostatAluSkinBarrel");
      if(skin){
	const TGeoShape* shape_skin = skin->GetVolume()->GetShape();
	if(shape_skin->TestShapeBit(TGeoTube::kGeoTube)){
          const TGeoTube* tube = (const TGeoTube*) shape_skin;
          rAlu  = tube->GetRmin()*CLHEP::cm;
          drAlu = tube->GetRmax()*CLHEP::cm - rAlu;
          aluHalfZ = tube->GetDz()*CLHEP::cm;
	  //info() << rmin << "," << rmax << "," << zhalf << endmsg;
        }
        else{
          error() << skin->GetName() << " is not a TGeoTube! Shape bits = " <<  shape_skin->TestShapeBits(0xFFFFFFFF) << endmsg;
        }
      }
      TGeoNode* foam = FindNode(daughter, "CryostatFoamBarrel");
      if(foam){
        const TGeoShape* shape_foam = foam->GetVolume()->GetShape();
        if(shape_foam->TestShapeBit(TGeoTube::kGeoTube)){
          const TGeoTube* tube = (const TGeoTube*) shape_foam;
          rSty = tube->GetRmin()*CLHEP::cm;
          drSty = tube->GetRmax()*CLHEP::cm - rSty;
          dzSty = tube->GetDz()*CLHEP::cm;
	  //info() << rmin << "," << rmax << "," << zhalf << endmsg;
        }
        else{
          error() << foam->GetName() << " is not a TGeoTube! Shape bits = " << shape_foam->TestShapeBits(0xFFFFFFFF) << endmsg;
        }
      }
      TGeoNode* skinEnd = FindNode(daughter, "CryostatAluSkinEndPlateInner");
      if(skinEnd){
        const TGeoShape* shape_skinEnd = skinEnd->GetVolume()->GetShape();
        if(shape_skinEnd->TestShapeBit(TGeoTube::kGeoTube)){
          const TGeoTube* tube = (const TGeoTube*) shape_skinEnd;
          rInner = tube->GetRmin()*CLHEP::cm;
          double rmax = tube->GetRmax()*CLHEP::cm;
          drAlu = tube->GetDz()*CLHEP::cm*2;
	  //info() << rmin << "," << rmax << "," << zhalf << endmsg;
        }
        else{
          error() << skinEnd->GetName() << " is not a TGeoTube! Shape bits = " << shape_skinEnd->TestShapeBits(0xFFFFFFFF) << endmsg;
        }
      }
      TGeoNode* shellEnd = FindNode(daughter, "EndPlateShell_outer");
      if(shellEnd){
        const TGeoShape* shape_shellEnd = shellEnd->GetVolume()->GetShape();
	if(shape_shellEnd->TestShapeBit(TGeoTube::kGeoTube)){
	  const TGeoTube* tube = (const TGeoTube*) shape_shellEnd;
          double rmin = tube->GetRmin()*CLHEP::cm;
          double rmax = tube->GetRmax()*CLHEP::cm;
          double zhalf = tube->GetDz()*CLHEP::cm;
          //info() << rmin << "," << rmax << "," << zhalf << endmsg;
	}
	else{
          error() << shellEnd->GetName() << " is not a TGeoTube! Shape bits = " << shape_shellEnd->TestShapeBits(0xFFFFFFFF) << endmsg;
        }
      }
      
    }
    else if(nodeName=="layer_assembly_0_1"){
      if(TGeoNode* side_band = FindNode(daughter, "ElectronicsBand")){
	const TGeoShape* shape_band = side_band->GetVolume()->GetShape();
	if(shape_band->IsA()==TGeoBBox::Class()){
	  const TGeoBBox* box = (const TGeoBBox*) shape_band;
	  side_band_electronics_width = box->GetDX()*CLHEP::cm*2;
	  //info() << "fucd: "<< box->GetDX() << " " << box->GetDY() << " " << box->GetDZ() <<endmsg;
	}
	else{
	  error() << "ElectronicsBand is not a TGeoBBox!!!"<< endmsg;
	}
      }
      if(TGeoNode* end = FindNode(daughter, "ElectronicsEnd")){
	const TGeoShape* shape_end = end->GetVolume()->GetShape();
	if(shape_end->IsA()==TGeoBBox::Class()){
	  const TGeoBBox* box = (const TGeoBBox*) shape_end;
	  end_electronics_half_z = box->GetDY()*CLHEP::cm;
	  //info() << "fucd: " << box->GetDX() << " " << box->GetDY() << " " << box->GetDZ() << endmsg;                                                                                                                              
	}
	else{
	  error() << "ElectronicsEnd is not a TGeoBBox!!!"<< endmsg;
	}
      }
    }
    
    /*
    for(int j=0;j<daughter->GetNdaughters();j++){
      TGeoNode* next = daughter->GetDaughter(j);
      info() << "fucd:  " << next->GetName() << endmsg;
    }
    */
  }
  
  const std::map<std::string, dd4hep::DetElement>& components = vxd.children();
  for(std::map<std::string, dd4hep::DetElement>::const_iterator it=components.begin();it!=components.end();it++){
    dd4hep::DetElement component = it->second;
    dd4hep::Volume vol = component.volume();
    dd4hep::PlacedVolume phys = component.placement();
    TGeoMaterial* mat = vol->GetMaterial();
    const TGeoShape* shape = vol->GetShape();
    const dd4hep::PlacedVolumeExtension::VolIDs& ids = phys.volIDs();
    //info() << " " << it->first << " " << vol->GetName() << " " << component.id() << " " << component.path() << " " << component.placementPath() << endmsg;
    //info() << "   " << shape->GetName() << " " << vol.solid().name() << endmsg; 
    //info() << "   " << mat->GetName() <<  " " << mat->GetRadLen() << endmsg;
    //info() << "   " << ids.str() << endmsg;
    //info() << "   " << vol->GetNdaughters() << endmsg; 
    //double radL = mat->GetRadLen();
    //dd4hep::Solid  solid = vol.solid();
    //info() << "  " <<  sh->TestShapeBit(TGeoShape::kGeoBox) <<  " " << sh->GetName() << " " << phys.material().radLength() << endmsg;
    if(vol.isSensitive()&&shape->IsA()==TGeoBBox::Class()){
      int iLayer  = ids.find("layer")->second;
      int iModule = ids.find("module")->second;
      int iSide   = ids.find("side")->second;
      //info() << "  layer=" << iLayer << " module=" << iModule << mat->GetName() << endmsg;
      if(iModule==0&&iLayer==helpCount+1){
	helpCount++;
	helpSensitives.push_back(thisLadder);
	helpLadders.push_back(thisLadder);
	helpNumberLadders.push_back(nLadders);
	helpPhi0.push_back(phi0);
	nLadders = 0;
	thisLadder.length = 0;
      }
      if(iLayer == helpCount){
	if(iModule == 0){
	  const TGeoBBox* box = (const TGeoBBox*) shape;
	  double width     = box->GetDX()*CLHEP::cm;
	  double length    = box->GetDY()*CLHEP::cm;
	  double thickness = box->GetDZ()*CLHEP::cm;
	  TGeoMatrix* matrix = phys->GetMatrix();
	  const double* pos = matrix->GetTranslation();
	  const double* rot_data = matrix->GetRotationMatrix();
	  TGeoRotation rot;
	  rot.SetMatrix(rot_data);
	  double theta,phi,psi;
	  rot.GetAngles(phi,theta,psi);
	  phi *= TMath::DegToRad();
	  theta *= TMath::DegToRad();
	  psi *= TMath::DegToRad();
	  phi0 = -TMath::PiOver2()+phi;
	  double distance = fabs(cos(phi0)*sin(theta)*pos[0]+sin(phi0)*sin(theta)*pos[1]+cos(theta)*pos[2]);
	  double offset = sqrt(pos[0]*pos[0]+pos[1]*pos[1]-distance*distance)*pos[0]/fabs(pos[0])*CLHEP::cm;
	  distance *= CLHEP::cm;
	  distance -= thickness;
	  double radL = mat->GetRadLen()*CLHEP::cm;
	  //info() << " ->   " << helpCount << ": " << distance << " " << cos(atan2(pos[1],pos[0])-phi)*sqrt(pos[0]*pos[0]+pos[1]*pos[1]) << endmsg;
	  thisLadder.distance  = distance;
	  thisLadder.offset    = offset;
	  thisLadder.thickness = thickness;
	  thisLadder.length   += length;
	  thisLadder.width     = width;
	  thisLadder.radLength = radL;
	  thisLadder.z         = pos[2]*CLHEP::cm;
	}
	if(iModule==nLadders) nLadders++; 
      }
    }
    //info() << "  " << vol.solid().type() << " " << vol.solid().name() << " " << vol->GetNdaughters() << endmsg;
    else if(it->first=="VXD_support"){
      helpCount++;
      helpSensitives.push_back(thisLadder);
      helpLadders.push_back(thisLadder);
      helpNumberLadders.push_back(nLadders);
      helpPhi0.push_back(phi0);
      nLadders = 0;
      if(vol->GetNdaughters()==0) error() << "!!!!!!!!!" << endmsg;

      int nFlexCable = 0, nFoamSpacer=0, nMetalTraces=0;
      int currentLayer = -1;
      double tFlexCable, tFoamSpacer, tMetalTraces;
      double radLFlexCable, radLFoamSpacer, radLMetalTraces;
      double intLFlexCable, intLFoamSpacer, intLMetalTraces;
      double dFlexCable, dFoamSpacer, dMetalTraces;
      double metalZeff, metalZAeff, foamZeff, foamZAeff, flexZeff, flexZAeff;
      for(int i=0;i<vol->GetNdaughters();i++){
	TGeoNode* daughter = vol->GetNode(i);
	TGeoMaterial* matDaughter = daughter->GetMedium()->GetMaterial();
	const TGeoShape* shape_sup = daughter->GetVolume()->GetShape();
	TGeoMatrix* matrix = daughter->GetMatrix();
	const double* pos = matrix->GetTranslation();
	const double* rot_data = matrix->GetRotationMatrix();
	TGeoRotation rot;
	rot.SetMatrix(rot_data);
	double theta,phi,psi;
	rot.GetAngles(phi,theta,psi);
	phi *= TMath::DegToRad();
	theta *= TMath::DegToRad();
	psi *= TMath::DegToRad();
	phi0 = -TMath::PiOver2()+phi;
	std::string phy_name = daughter->GetName();
	if(phy_name.find("FoamSpacer")==-1&&phy_name.find("FlexCable")==-1&&phy_name.find("MetalTraces")==-1){
	  //info() << phy_name <<endmsg;
	  continue;
	}
	int iLayer = atoi(phy_name.substr(phy_name.find("_")+1,2).c_str());
	if(iLayer!=currentLayer){
	  //info() << tFoamSpacer << "," << tFlexCable << "," << tMetalTraces << endmsg;
	  helpLadders[currentLayer].thickness = tFoamSpacer+tFlexCable+tMetalTraces;
	  helpLadders[currentLayer].radLength = helpLadders[currentLayer].thickness / (tFoamSpacer/radLFoamSpacer+tFlexCable/radLFlexCable+tMetalTraces/radLMetalTraces);
	  nFlexCable = 0;
	  nFoamSpacer=0;
	  nMetalTraces=0;
	  currentLayer=iLayer;
	}
	//info() << "ss   pos=" << pos[0] << "," << pos[1] << "," << pos[2] << " distance=" << sqrt(pos[0]*pos[0]+pos[1]*pos[1]) << endmsg;
	//info() << "ss   rot=" << phi << "," << theta << "," << psi << endmsg;
	//info() << "ss " << daughter->GetName() << " " << daughter->GetVolume()->GetName() << " " << endmsg;
	if(shape_sup->IsA()==TGeoBBox::Class()&&(nFoamSpacer==0||nFlexCable==0||nMetalTraces==0)){
	  const TGeoBBox* box = (const TGeoBBox*) shape_sup;
	  //info() << phy_name.substr(phy_name.find("_")+1,2) << " " << iLayer << " " << box->GetDX() << "," << box->GetDY() << "," << box->GetDZ() << endmsg;
	  //info() << "fucd: pos " << pos[0] << " " << pos[1] << " " << pos[2] << endmsg; 
	  double distance = fabs(cos(phi0)*sin(theta)*pos[0]+sin(phi0)*sin(theta)*pos[1]+cos(theta)*pos[2]);
          double offset = sqrt(pos[0]*pos[0]+pos[1]*pos[1]-distance*distance)*pos[0]/fabs(pos[0])*CLHEP::cm;
          distance -= box->GetDZ();
	  distance *= CLHEP::cm;
	  if(helpLadders[iLayer].distance == helpSensitives[iLayer].distance) helpLadders[iLayer].distance = distance;
	  else helpLadders[iLayer].distance = TMath::Min(helpLadders[iLayer].distance, distance);
	  //info() << phy_name << " " << distance << " " << offset << endmsg;
	  if(phy_name.find("FoamSpacer")!=-1&&nFoamSpacer==0){
	    helpLadders[iLayer].offset    = offset;
	    tFoamSpacer = box->GetDZ()*CLHEP::cm;
	    radLFoamSpacer = matDaughter->GetRadLen()*CLHEP::cm;
	    intLFoamSpacer = matDaughter->GetIntLen()*CLHEP::cm;
	    dFoamSpacer = matDaughter->GetDensity()*CLHEP::g/CLHEP::cm3;
	    //fucd: A calculated by TGeoMaterial class is not equal to Zeff/ZAeff, Zeff = sum(Zi*Ai/sumA), ZAeff = sum(Zi/Ai*A/totalA)
	    //      use Zeff and ZAeff to keep same with Mokka case  
	    //foamZ = matDaughter->GetZ();
	    //foamA = matDaughter->GetA();
	    double totalA = 0, Zeff = 0, ZAeff = 0; 
	    for(int iEle = 0; iEle<matDaughter->GetNelements(); iEle++){
	      totalA += matDaughter->GetElement(iEle)->A();
	    }
	    for(int iEle = 0; iEle<matDaughter->GetNelements(); iEle++){
	      double A, Z, w;
	      // by fucd: w!=A/totalA, strange! to fix
	      matDaughter->GetElementProp(A,Z,w,iEle);
              Zeff  += Z*w;
	      ZAeff += Z/A*w;
	      //info() << std::setprecision(16) << Z << " " << A << " " << A/totalA << " " << w << endmsg;
            }
	    foamZeff = Zeff;
	    foamZAeff = ZAeff;
	    nFoamSpacer++;
	  }
	  if(phy_name.find("FlexCable")!=-1&&nFlexCable==0){
	    helpLadders[iLayer].width     = box->GetDX()*CLHEP::cm;
            helpLadders[iLayer].length    = box->GetDY()*CLHEP::cm-beryllium_ladder_block_length*2-end_electronics_half_z*2;
	    tFlexCable = box->GetDZ()*CLHEP::cm;
	    radLFlexCable = matDaughter->GetRadLen()*CLHEP::cm;
	    intLFlexCable = matDaughter->GetIntLen()*CLHEP::cm;
	    dFlexCable = matDaughter->GetDensity()*CLHEP::g/CLHEP::cm3;
	    double Zeff = 0, ZAeff = 0;
            for(int iEle = 0; iEle<matDaughter->GetNelements(); iEle++){
              double A, Z, w;
	      matDaughter->GetElementProp(A,Z,w,iEle);
              Zeff  += Z*w;
              ZAeff += Z/A*w;
	      //std::cout << std::setprecision(16) << Z << " " << A << " " << w << std::endl;
            }
	    flexZeff  = Zeff;
	    flexZAeff = ZAeff;
	    nFlexCable++;
	  }
	  if(phy_name.find("MetalTraces")!=-1&&nMetalTraces==0){
	    tMetalTraces = box->GetDZ()*CLHEP::cm;
	    radLMetalTraces = matDaughter->GetRadLen()*CLHEP::cm;
	    intLMetalTraces = matDaughter->GetIntLen()*CLHEP::cm;
	    dMetalTraces = matDaughter->GetDensity()*CLHEP::g/CLHEP::cm3;
	    double totalA = 0, Zeff = 0, ZAeff = 0;
            for(int iEle = 0; iEle<matDaughter->GetNelements(); iEle++){
              totalA += matDaughter->GetElement(iEle)->A();
            }
            for(int iEle = 0; iEle<matDaughter->GetNelements(); iEle++){
              double A, Z, w;
	      matDaughter->GetElementProp(A,Z,w,iEle);
	      Zeff  += Z*w;
              ZAeff += Z/A*w;
	      //info() << Z << " " << A << " " << w << endmsg;
            }
            metalZeff  = Zeff;
            metalZAeff = ZAeff;
	    nMetalTraces++;
	  }
	}
      }
      {
	//info() << tFoamSpacer << "," << tFlexCable << "," << tMetalTraces << endmsg;
	double tSupport = tMetalTraces + tFoamSpacer + tFlexCable;
	helpLadders[currentLayer].thickness = tSupport;
	helpLadders[currentLayer].radLength = helpLadders[currentLayer].thickness / (tFoamSpacer/radLFoamSpacer+tFlexCable/radLFlexCable+tMetalTraces/radLMetalTraces);
	nFlexCable = 0;
	nFoamSpacer=0;
	nMetalTraces=0;

	//calculations of thickness fractions of each layer of the support
	double metalTF = tMetalTraces / tSupport;
	double foamTF = tFoamSpacer / tSupport;
	double flexTF = tFlexCable / tSupport;
	//info() << foamTF << "," << flexTF << "," << metalTF << endmsg;
	//info() << dFoamSpacer/(CLHEP::kg/CLHEP::cm3) << "," << dFlexCable/(CLHEP::kg/CLHEP::cm3) << "," << dMetalTraces/(CLHEP::kg/CLHEP::cm3) << endmsg;
	//info() << foamZeff << " " << flexZeff << " " << metalZeff << endmsg;
	//info() << foamZAeff << " " << flexZAeff << " " << metalZAeff << endmsg;
	G4double elemVol = 1*CLHEP::cm3;
	G4double VXDSupportMass = (foamTF*dFoamSpacer + flexTF*dFlexCable + metalTF*dMetalTraces)*elemVol;
	VXDSupportDensity = VXDSupportMass/elemVol;
	
	G4double foamFM = 100. * ((foamTF*(elemVol)*dFoamSpacer) / VXDSupportMass) ;
	G4double kaptonFM = 100. * ((flexTF*(elemVol)*dFlexCable) / VXDSupportMass) ;
	G4double metalFM = 100. * ((metalTF*(elemVol)*dMetalTraces) / VXDSupportMass) ;

	VXDSupportRadLen = helpLadders[currentLayer].radLength;

	VXDSupportZeff = (metalFM/100.)*metalZeff + (kaptonFM/100.)*flexZeff + (foamFM/100.)*foamZeff;
	G4double VXDSupportZAeff = (metalFM/100.)*metalZAeff + (kaptonFM/100.)*flexZAeff + (foamFM/100.)*foamZAeff;
	VXDSupportAeff = VXDSupportZeff / VXDSupportZAeff;
	G4double VXDSupportIntLength = 1. / ((metalTF/intLMetalTraces) + (flexTF/intLFlexCable) + (foamTF/intLFoamSpacer));
	VXDSupportDensity = 1000*VXDSupportDensity/(CLHEP::g/CLHEP::cm3);
	//info() << "fucd: " << VXDSupportZeff << " " << VXDSupportAeff << " " << VXDSupportRadLen << " " << VXDSupportIntLength << " " << VXDSupportDensity << endmsg; 
	//info() << intLMetalTraces << " " << intLFlexCable << " " << intLFoamSpacer <<endmsg;
      }
    }
    //info() << it->first << endmsg; 
  }
  if(end_electronics_half_z>0 && side_band_electronics_width==0) type = gear::ZPlanarParametersImpl::CCD  ;
  if(side_band_electronics_width>0 && end_electronics_half_z==0 ) type = gear::ZPlanarParametersImpl::CMOS ;
  if(side_band_electronics_width>0 && end_electronics_half_z>0) type = gear::ZPlanarParametersImpl::HYBRID ;
  
  m_vxdParameters = new gear::ZPlanarParametersImpl(type, shellInnerRadius, shellOuterRadius, shellHalfLength, gap, shellRadLength );
  // by fucd: debug info, if validated enough, merge them in future
  info() << "=====================from convertor==============================" << endmsg;
  info() << type << " " << shellInnerRadius << " " << shellOuterRadius << " " << shellHalfLength << " " << gap << " " << shellRadLength << endmsg;
  for(int i=0;i<helpCount;i++){
    m_vxdParameters->addLayer(helpNumberLadders[i] , helpPhi0[i] ,
			    helpLadders[i].distance , helpLadders[i].offset, helpLadders[i].thickness*2 ,
			    helpLadders[i].length , helpLadders[i].width*2 , helpLadders[i].radLength ,
			    helpSensitives[i].distance, helpSensitives[i].offset , helpSensitives[i].thickness*2 ,
			    helpSensitives[i].length , helpSensitives[i].width*2 , helpSensitives[i].radLength ) ;
    info() << "fucd " << i << ": " << helpNumberLadders[i] << ", " << helpPhi0[i] << ", "
	   << helpLadders[i].distance << ", " << helpLadders[i].offset << ", " << helpLadders[i].thickness*2 << ", " << helpLadders[i].length << ", "
	   << helpLadders[i].width*2 << ", " << helpLadders[i].radLength << ", " << helpSensitives[i].distance << ", " << helpSensitives[i].offset << ", "
	   << helpSensitives[i].thickness*2 << ", " << helpSensitives[i].length << ", " << helpSensitives[i].width*2 << ", " << helpSensitives[i].radLength << endmsg;
  }
  //m_vxdInfra = new gear::GearParametersImpl;
  //CryostatAlRadius, CryostatAlThickness, CryostatAlInnerR, CryostatAlZEndCap, CryostatAlHalfZ
  //m_vxdInfra->setDoubleVal("CryostatAlRadius",rAlu);
  //m_vxdInfra->setDoubleVal("CryostatAlThickness",drAlu);
  //m_vxdInfra->setDoubleVal("CryostatAlInnerR",rInner);
  //m_vxdInfra->setDoubleVal("CryostatAlZEndCap",aluEndcapZ=dzSty + drSty + drAlu / 2);
  //m_vxdInfra->setDoubleVal("CryostatAlHalfZ",aluHalfZ= dzSty + drSty);
  // change GearParametersImpl to map
  std::map<std::string,double> vxdInfra;
  vxdInfra["CryostatAlRadius"]    = rAlu;
  vxdInfra["CryostatAlThickness"] = drAlu;
  vxdInfra["CryostatAlInnerR"]    = rInner;
  vxdInfra["CryostatAlZEndCap"]   = dzSty+drSty+drAlu/2;
  vxdInfra["CryostatAlHalfZ"]     = dzSty+drSty;
  m_detParameters["VXDInfra"] = vxdInfra;
  //effective A different with what in Mokka, fix them as Mokka's 
  //m_materials["VXDSupportMaterial"] = new TMaterial("VXDSupportMaterial", "", VXDSupportAeff, VXDSupportZeff, VXDSupportDensity, VXDSupportRadLen, 0.);
  m_materials["VXDSupportMaterial"] = new TMaterial("VXDSupportMaterial", "", 2.075865162e+01, 1.039383117e+01, 2.765900000e+02/1000, 1.014262421e+02, 0.);

  info() << "=====================from ZPlanarData==============================" << endmsg;
  if(m_vxdData){
    info() << m_vxdData->rInnerShell << " " << m_vxdData->rOuterShell << " " << m_vxdData->zHalfShell << " " << m_vxdData->gapShell << endmsg;
    const std::vector<dd4hep::rec::ZPlanarData::LayerLayout>& layers = m_vxdData->layers;
    for(int i=0;i<layers.size();i++){
      const dd4hep::rec::ZPlanarData::LayerLayout& thisLayer = layers[i];
      info() << i << ": " << thisLayer.ladderNumber << "," << thisLayer.phi0 << "," << thisLayer.distanceSupport << "," << thisLayer.offsetSupport << ","
	     << thisLayer.thicknessSupport << "," << thisLayer.zHalfSupport << "," << thisLayer.widthSupport << "," << "NULL," 
	     << thisLayer.distanceSensitive << "," << thisLayer.offsetSensitive << "," << thisLayer.thicknessSensitive << "," << thisLayer.zHalfSensitive << ","
	     << thisLayer.widthSensitive << ",NULL" << endmsg;
    }
  }
  info() << rAlu << " " << drAlu << " " << rInner << " " << aluEndcapZ << " " << aluHalfZ << endmsg; 
  //info() << m_materials["VXDSupportMaterial"] << endmsg; 
  return sc;
}

TGeoNode* GeomSvc::FindNode(TGeoNode* mother, char* name){
  TGeoNode* next = 0;
  if(mother->GetNdaughters()!=0){ 
    for(int i=0;i<mother->GetNdaughters();i++){
      TGeoNode* daughter = mother->GetDaughter(i);
      std::string s = daughter->GetName();
      //info() << "current: " << s << " search for" << name << endmsg;
      if(s.find(name)!=-1){
	next = daughter;
	break;
      }
      else{
	next = FindNode(daughter, name);
      }
    }
  }
  return next;
}

TMaterial* GeomSvc::getMaterial(std::string name){
  std::map<std::string, TMaterial*>::const_iterator it = m_materials.find(name);
  if(it!=m_materials.end()) return it->second;
  else return 0;     

}
