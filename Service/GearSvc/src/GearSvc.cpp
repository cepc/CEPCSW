#include "GearSvc.h"
#include "DetInterface/IGeomSvc.h"
#include "DetSegmentation/GridDriftChamber.h"

#include "gearxml/GearXML.h"
#include "gearimpl/GearMgrImpl.h"
#include "gearimpl/ConstantBField.h"
#include "gearimpl/ZPlanarParametersImpl.h"
#include "gearimpl/ZPlanarLayerLayoutImpl.h"
#include "gearimpl/FTDParametersImpl.h"
#include "gearimpl/TPCParametersImpl.h"
#include "gearimpl/FixedPadSizeDiskLayout.h"
#include "gearimpl/CalorimeterParametersImpl.h"
#include "gearimpl/SimpleMaterialImpl.h"
#include "gearxml/tinyxml.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h"
#include "CLHEP/Units/SystemOfUnits.h"

struct helpLayer {
  double distance =0;
  double offset =0;
  double thickness =0;
  double length =0;
  double width =0;
  double radLength =0;
  double z =0;
  double foam_spacer_radLength =0;
};

static const double deg_to_rad = dd4hep::degree/CLHEP::rad;
static const double rad_to_deg = dd4hep::rad/CLHEP::degree;

DECLARE_COMPONENT(GearSvc)

GearSvc::GearSvc(const std::string& name, ISvcLocator* svc)
    : base_class(name, svc),
      m_gearMgr(nullptr)
{
}

GearSvc::~GearSvc()
{
}

gear::GearMgr* GearSvc::getGearMgr()
{
    return m_gearMgr;
}

StatusCode GearSvc::initialize()
{
  StatusCode sc;
  
  if ( m_gearFile.size() > 0 ) {
    info() << "instantiated GEAR from file " << m_gearFile << endmsg;
    m_gearMgr = gear::GearXML(m_gearFile).createGearMgr();
  }
  else {
    warning() << "no GEAR XML file given ..." << endmsg;
    m_gearMgr = new gear::GearMgrImpl;
    
    auto geomSvc = service<IGeomSvc>("GeomSvc");
    if ( !geomSvc ) {
      info() << "Failed to find GeomSvc ..." << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "Fill GEAR data from GeomSvc" << endmsg;
    m_gearMgr->setDetectorName("CRD_o1_v01");

    const dd4hep::Direction& field = geomSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
    gear::ConstantBField* b = new gear::ConstantBField(gear::Vector3D(field.x()/dd4hep::tesla, field.y()/dd4hep::tesla, field.z()/dd4hep::tesla));
    m_gearMgr->setBField(b);

    dd4hep::DetElement world = geomSvc->getDD4HepGeo();
    const std::map<std::string, dd4hep::DetElement>& subs = world.children();
    for(std::map<std::string, dd4hep::DetElement>::const_iterator it=subs.begin();it!=subs.end();it++){
      dd4hep::DetElement sub = it->second;
      info() << it->first << " " << sub.path() << " " << sub.placementPath() << endmsg;
      if(it->first=="Tube"||it->first=="BeamPipe"){
	sc = convertBeamPipe(sub);
      }
      else if(it->first=="VXD"){
	sc = convertVXD(sub);
      }
      else if(it->first=="FTD"){
        sc = convertFTD(sub);
      }
      else if(it->first=="SIT"){
	sc = convertSIT(sub);
      }
      else if(it->first=="TPC"){
	sc = convertTPC(sub);
      }
      else if(it->first=="DriftChamber"){
        sc = convertDC(sub);
      }
      else if(it->first=="SET"){
	sc = convertSET(sub);
      }
      else{
	info() << it->first << " will convert in future! now fake parameters" << endmsg;
      }
      if(sc==StatusCode::FAILURE) return sc;
    }

    gear::CalorimeterParametersImpl* barrelParam = new gear::CalorimeterParametersImpl(1847.415655, 2350., 8, 0.);
    gear::CalorimeterParametersImpl* endcapParam = new gear::CalorimeterParametersImpl(400., 2088.8, 2450., 2, 0.);
    for(int i=0;i<29;i++){
      if(i<19){
	barrelParam->layerLayout().positionLayer(0, 5.25, 1.016666667e+01, 1.016666667e+01, 2.1);
	endcapParam->layerLayout().positionLayer(0, 5.25, 1.016666667e+01, 1.016666667e+01, 2.1);
      }
      else if(i<20){
	barrelParam->layerLayout().positionLayer(0, 6.3, 1.016666667e+01, 1.016666667e+01, 2.1);
	endcapParam->layerLayout().positionLayer(0, 6.3, 1.016666667e+01, 1.016666667e+01, 2.1);
      }
      else{
	barrelParam->layerLayout().positionLayer(0, 4.2, 1.016666667e+01, 1.016666667e+01, 4.2);
	endcapParam->layerLayout().positionLayer(0, 4.2, 1.016666667e+01, 1.016666667e+01, 4.2);
      }
    }
    m_gearMgr->setEcalBarrelParameters(barrelParam);
    m_gearMgr->setEcalEndcapParameters(endcapParam);

    gear::CalorimeterParametersImpl* barrelYokeParam = new gear::CalorimeterParametersImpl(4173.929932, 4072., 12, 0.0);
    gear::CalorimeterParametersImpl* endcapYokeParam = new gear::CalorimeterParametersImpl(320., 7414.929932, 4072., 2, 0.0);
    gear::CalorimeterParametersImpl* plugYokeParam   = new gear::CalorimeterParametersImpl(320., 2849.254326, 3781.43, 2, 0.0);
    plugYokeParam->setDoubleVal("YokePlugThickness", 290.57) ;
    m_gearMgr->setYokeBarrelParameters(barrelYokeParam) ;
    m_gearMgr->setYokeEndcapParameters(endcapYokeParam) ;
    m_gearMgr->setYokePlugParameters(plugYokeParam) ;
    gear::TiXmlDocument* doc = new gear::TiXmlDocument ;
    gear::GearXML::createXMLFile(m_gearMgr, "test.xml");
  }
  
  return StatusCode::SUCCESS;
}

StatusCode GearSvc::finalize()
{
    if ( m_gearMgr ) {
        delete m_gearMgr;
        m_gearMgr = nullptr;
    }

    return StatusCode::SUCCESS;
}

StatusCode GearSvc::convertBeamPipe(dd4hep::DetElement& pipe){
  StatusCode sc;

  dd4hep::rec::ConicalSupportData* beamPipeData = nullptr;
  try{
    beamPipeData = pipe.extension<dd4hep::rec::ConicalSupportData>();
  }
  catch(std::runtime_error& e){
    warning() << e.what() << " " << beamPipeData << endmsg;
    return StatusCode::FAILURE;
  }

  std::vector<double> gearValRInner;
  std::vector<double> gearValROuter;
  std::vector<double> gearValZ;
  const std::vector<dd4hep::rec::ConicalSupportData::Section>& sections = beamPipeData->sections;
  for(int i=0;i<sections.size();i++){
    gearValZ.push_back(sections[i].zPos*CLHEP::cm );
    gearValRInner.push_back(sections[i].rInner*CLHEP::cm );
    gearValROuter.push_back(sections[i].rOuter*CLHEP::cm );
  }

  gear::GearParametersImpl* gearParameters = new gear::GearParametersImpl;
  gearParameters -> setDoubleVals( "Z" , gearValZ ) ;
  gearParameters -> setDoubleVals( "RInner" , gearValRInner ) ;
  gearParameters -> setDoubleVals( "ROuter" , gearValROuter ) ;

  m_gearMgr->setGearParameters("BeamPipe", gearParameters ) ;

  return StatusCode::SUCCESS;
}

StatusCode GearSvc::convertVXD(dd4hep::DetElement& vxd){
  StatusCode sc;
  //fucd: another method to obtain parameters, but not fully for KalDet
  dd4hep::rec::ZPlanarData* vxdData = nullptr;
  bool extensionDataValid = true;
  try{
    vxdData = vxd.extension<dd4hep::rec::ZPlanarData>();
  }
  catch(std::runtime_error& e){
    extensionDataValid = false;
    info() << e.what() << " " << vxdData << endmsg;
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
  double rAlu=0, drAlu, rSty, drSty, dzSty, rInner, aluEndcapZ, aluHalfZ, alu_RadLen, Cryostat_dEdx;
  double VXDSupportDensity, VXDSupportZeff, VXDSupportAeff, VXDSupportRadLen, VXDSupportIntLen=0;
  double styDensity, styZeff, styAeff, styRadLen, styIntLen; 
  dd4hep::Volume vxd_vol = vxd.volume();
  for(int i=0;i<vxd_vol->GetNdaughters();i++){
    TGeoNode* daughter = vxd_vol->GetNode(i);
    std::string nodeName = daughter->GetName();
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
        }
        else{
          error() << foam->GetName() << " is not a TGeoTube! Shape bits = " << shape_foam->TestShapeBits(0xFFFFFFFF) << endmsg;
        }
	TGeoMaterial* mat = foam->GetMedium()->GetMaterial();
	double Zeff = 0, ZAeff = 0;
	for(int iEle = 0; iEle<mat->GetNelements(); iEle++){
	  double A, Z, w;
	  mat->GetElementProp(A,Z,w,iEle);
	  Zeff  += Z*w;
	  ZAeff += Z/A*w;
	  //std::cout << std::setprecision(16) << Z << " " << A << " " << w << std::endl;
	}
	styZeff    = Zeff;
	styAeff    = Zeff/ZAeff;
	styRadLen  = mat->GetRadLen()*CLHEP::cm;
	styIntLen  = mat->GetIntLen()*CLHEP::cm;
	styDensity = mat->GetDensity();
      }
      TGeoNode* skinEnd = FindNode(daughter, "CryostatAluSkinEndPlateInner");
      if(skinEnd){
        const TGeoShape* shape_skinEnd = skinEnd->GetVolume()->GetShape();
        if(shape_skinEnd->TestShapeBit(TGeoTube::kGeoTube)){
          const TGeoTube* tube = (const TGeoTube*) shape_skinEnd;
          rInner = tube->GetRmin()*CLHEP::cm;
          double rmax = tube->GetRmax()*CLHEP::cm;
          drAlu = tube->GetDz()*CLHEP::cm*2;
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
  }

  const std::map<std::string, dd4hep::DetElement>& components = vxd.children();
  for(std::map<std::string, dd4hep::DetElement>::const_iterator it=components.begin();it!=components.end();it++){
    dd4hep::DetElement component = it->second;
    dd4hep::Volume vol = component.volume();
    dd4hep::PlacedVolume phys = component.placement();
    TGeoMaterial* mat = vol->GetMaterial();
    const TGeoShape* shape = vol->GetShape();
    const dd4hep::PlacedVolumeExtension::VolIDs& ids = phys.volIDs();
    if(vol.isSensitive()&&shape->IsA()==TGeoBBox::Class()){
      int iLayer  = ids.find("layer")->second;
      int iModule = ids.find("module")->second;
      int iSide   = ids.find("side")->second;
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
          phi *= deg_to_rad;
          theta *= deg_to_rad;
          psi *= deg_to_rad;
          phi0 = -dd4hep::halfpi+phi;
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
        phi *= deg_to_rad;
        theta *= deg_to_rad;
        psi *= deg_to_rad;
        phi0 = -CLHEP::halfpi+phi;
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
	if(shape_sup->IsA()==TGeoBBox::Class()&&(nFoamSpacer==0||nFlexCable==0||nMetalTraces==0)){
          const TGeoBBox* box = (const TGeoBBox*) shape_sup;
	  double distance = fabs(cos(phi0)*sin(theta)*pos[0]+sin(phi0)*sin(theta)*pos[1]+cos(theta)*pos[2]);
          double offset = sqrt(pos[0]*pos[0]+pos[1]*pos[1]-distance*distance)*pos[0]/fabs(pos[0])*CLHEP::cm;
          distance -= box->GetDZ();
          distance *= CLHEP::cm;
          if(helpLadders[iLayer].distance == helpSensitives[iLayer].distance) helpLadders[iLayer].distance = distance;
          else helpLadders[iLayer].distance = TMath::Min(helpLadders[iLayer].distance, distance);
          if(phy_name.find("FoamSpacer")!=-1&&nFoamSpacer==0){
            helpLadders[iLayer].offset    = offset;
            tFoamSpacer = box->GetDZ()*CLHEP::cm;
            radLFoamSpacer = matDaughter->GetRadLen()*CLHEP::cm;
            intLFoamSpacer = matDaughter->GetIntLen()*CLHEP::cm;
            dFoamSpacer = matDaughter->GetDensity();
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
            dFlexCable = matDaughter->GetDensity();
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
            dMetalTraces = matDaughter->GetDensity();
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
        double elemVol = 1*CLHEP::cm3;
        double VXDSupportMass = (foamTF*dFoamSpacer + flexTF*dFlexCable + metalTF*dMetalTraces)*elemVol;
        VXDSupportDensity = VXDSupportMass/elemVol;
	double foamFM = 100. * ((foamTF*(elemVol)*dFoamSpacer) / VXDSupportMass) ;
        double kaptonFM = 100. * ((flexTF*(elemVol)*dFlexCable) / VXDSupportMass) ;
        double metalFM = 100. * ((metalTF*(elemVol)*dMetalTraces) / VXDSupportMass) ;

        VXDSupportRadLen = helpLadders[currentLayer].radLength;

        VXDSupportZeff = (metalFM/100.)*metalZeff + (kaptonFM/100.)*flexZeff + (foamFM/100.)*foamZeff;
        double VXDSupportZAeff = (metalFM/100.)*metalZAeff + (kaptonFM/100.)*flexZAeff + (foamFM/100.)*foamZAeff;
        VXDSupportAeff = VXDSupportZeff / VXDSupportZAeff;
        double VXDSupportIntLength = 1. / ((metalTF/intLMetalTraces) + (flexTF/intLFlexCable) + (foamTF/intLFoamSpacer));
        VXDSupportDensity = VXDSupportDensity*(CLHEP::g/CLHEP::cm3)/(CLHEP::kg/CLHEP::m3);
        //info() << "fucd: " << VXDSupportZeff << " " << VXDSupportAeff << " " << VXDSupportRadLen << " " << VXDSupportIntLength << " " << VXDSupportDensity << endmsg;
        //info() << intLMetalTraces << " " << intLFlexCable << " " << intLFoamSpacer <<endmsg;
      }
    }
    //info() << it->first << endmsg;
  }
  if(end_electronics_half_z>0 && side_band_electronics_width==0) type = gear::ZPlanarParametersImpl::CCD  ;
  if(side_band_electronics_width>0 && end_electronics_half_z==0 ) type = gear::ZPlanarParametersImpl::CMOS ;
  if(side_band_electronics_width>0 && end_electronics_half_z>0) type = gear::ZPlanarParametersImpl::HYBRID ;
  gear::ZPlanarParametersImpl* vxdParameters = new gear::ZPlanarParametersImpl(type, shellInnerRadius, shellOuterRadius, shellHalfLength, gap, shellRadLength);
  // by fucd: debug info, if validated enough, merge them in future
  info() << "=====================from convertor==============================" << endmsg;
  info() << type << " " << shellInnerRadius << " " << shellOuterRadius << " " << shellHalfLength << " " << gap << " " << shellRadLength << endmsg;
  for(int i=0;i<helpCount;i++){
    vxdParameters->addLayer(helpNumberLadders[i] , helpPhi0[i] ,
                            helpLadders[i].distance , helpLadders[i].offset, helpLadders[i].thickness*2 ,
                            helpLadders[i].length , helpLadders[i].width*2 , helpLadders[i].radLength ,
                            helpSensitives[i].distance, helpSensitives[i].offset , helpSensitives[i].thickness*2 ,
                            helpSensitives[i].length , helpSensitives[i].width*2 , helpSensitives[i].radLength ) ;
    info() << "fucd " << i << ": " << helpNumberLadders[i] << ", " << helpPhi0[i] << ", "
           << helpLadders[i].distance << ", " << helpLadders[i].offset << ", " << helpLadders[i].thickness*2 << ", " << helpLadders[i].length << ", "
           << helpLadders[i].width*2 << ", " << helpLadders[i].radLength << ", " << helpSensitives[i].distance << ", " << helpSensitives[i].offset << ", "
           << helpSensitives[i].thickness*2 << ", " << helpSensitives[i].length << ", " << helpSensitives[i].width*2 << ", " << helpSensitives[i].radLength << endmsg;
  }
  m_gearMgr->setVXDParameters(vxdParameters) ;
  if(rAlu!=0){
    // rAlu=0, denote no cryostat 
    gear::GearParametersImpl* gearParameters = new gear::GearParametersImpl;
    //CryostatAlRadius, CryostatAlThickness, CryostatAlInnerR, CryostatAlZEndCap, CryostatAlHalfZ
    gearParameters->setDoubleVal("CryostatAlRadius",    rAlu);
    gearParameters->setDoubleVal("CryostatAlThickness", drAlu);
    gearParameters->setDoubleVal("CryostatAlInnerR",    rInner);
    gearParameters->setDoubleVal("CryostatAlZEndCap",   aluEndcapZ = dzSty+drSty+drAlu/2);
    gearParameters->setDoubleVal("CryostatAlHalfZ",     dzSty+drSty);
    m_gearMgr->setGearParameters("VXDInfra", gearParameters);
  }
  //effective A different with what in Mokka, fix them as Mokka's
  gear::SimpleMaterialImpl* VXDFoamShellMaterial_old = new gear::SimpleMaterialImpl("VXDFoamShellMaterial_old", 1.043890843e+01, 5.612886646e+00, 2.500000000e+01, 1.751650267e+04, 0);
  m_gearMgr->registerSimpleMaterial(VXDFoamShellMaterial_old);
  gear::SimpleMaterialImpl* VXDFoamShellMaterial = new gear::SimpleMaterialImpl("VXDFoamShellMaterial", styAeff, styZeff, styDensity*(CLHEP::g/CLHEP::cm3)/(CLHEP::kg/CLHEP::m3),
										styRadLen, styIntLen);
  m_gearMgr->registerSimpleMaterial(VXDFoamShellMaterial);
  gear::SimpleMaterialImpl* VXDSupportMaterial_old = new gear::SimpleMaterialImpl("VXDSupportMaterial_old", 2.075865162e+01, 1.039383117e+01, 2.765900000e+02, 1.014262421e+03, 3.341388059e+03);
  m_gearMgr->registerSimpleMaterial(VXDSupportMaterial_old);
  gear::SimpleMaterialImpl* VXDSupportMaterial = new gear::SimpleMaterialImpl("VXDSupportMaterial", VXDSupportAeff, VXDSupportZeff, VXDSupportDensity, VXDSupportRadLen, VXDSupportIntLen);
  m_gearMgr->registerSimpleMaterial(VXDSupportMaterial);
  info() << "=====================from ZPlanarData==============================" << endmsg;
  if(vxdData){
    info() << vxdData->rInnerShell << " " << vxdData->rOuterShell << " " << vxdData->zHalfShell << " " << vxdData->gapShell << endmsg;
    const std::vector<dd4hep::rec::ZPlanarData::LayerLayout>& layers = vxdData->layers;
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

StatusCode GearSvc::convertFTD(dd4hep::DetElement& ftd){
  dd4hep::rec::ZDiskPetalsData* ftdData = nullptr;
  try{
    ftdData = ftd.extension<dd4hep::rec::ZDiskPetalsData>();
  }
  catch(std::runtime_error& e){
    warning() << e.what() << " " << ftdData << endmsg;
    return StatusCode::FAILURE;
  }

  std::vector<dd4hep::rec::ZDiskPetalsData::LayerLayout>& ftdlayers = ftdData->layers;
  int nLayers = ftdlayers.size();

  gear::FTDParametersImpl* ftdParam = new gear::FTDParametersImpl();
  ftdParam->setDoubleVal("strip_width_mm", ftdData->widthStrip*CLHEP::cm);
  ftdParam->setDoubleVal("strip_length_mm", ftdData->lengthStrip*CLHEP::cm);
  ftdParam->setDoubleVal("strip_pitch_mm", ftdData->pitchStrip*CLHEP::cm);
  ftdParam->setDoubleVal("strip_angle_deg", ftdData->angleStrip*rad_to_deg);
  for(int layer = 0; layer < nLayers; layer++){
    dd4hep::rec::ZDiskPetalsData::LayerLayout& ftdlayer = ftdlayers[layer];
    int nPetals = ftdlayer.petalNumber;
    double dphi = CLHEP::twopi/nPetals;
    double phi0 = ftdlayer.phi0;
    double alpha = ftdlayer.alphaPetal;
    double zposition = ftdlayer.zPosition*CLHEP::cm;
    double zoffset = ftdlayer.zOffsetSupport*CLHEP::cm;
    int signoffset = ftdlayer.zOffsetSupport>0?1:-1;
    zoffset *= signoffset;

    double supRinner    = ftdlayer.distanceSupport*CLHEP::cm;
    double supThickness = ftdlayer.thicknessSupport*CLHEP::cm;
    double supLengthMin = ftdlayer.widthInnerSupport*CLHEP::cm;
    double supLengthMax = ftdlayer.widthOuterSupport*CLHEP::cm;
    double supWidth     = ftdlayer.lengthSupport*CLHEP::cm;
    double senRinner    = ftdlayer.distanceSensitive*CLHEP::cm;
    double senThickness = ftdlayer.thicknessSensitive*CLHEP::cm;
    double senLengthMin = ftdlayer.widthInnerSensitive*CLHEP::cm;
    double senLengthMax = ftdlayer.widthOuterSensitive*CLHEP::cm;
    double senWidth     = ftdlayer.lengthSensitive*CLHEP::cm;

    bool isDoubleSided  = ftdlayer.typeFlags[dd4hep::rec::ZDiskPetalsData::SensorType::DoubleSided];
    bool isPixelReadout = (bool)ftdlayer.typeFlags[dd4hep::rec::ZDiskPetalsData::SensorType::Pixel];
    int sensorType = (isPixelReadout)?gear::FTDParameters::PIXEL:gear::FTDParameters::STRIP;
    int nSensors = ftdlayer.sensorsPerPetal;
    double phalfangle = ftdlayer.petalHalfAngle;

    ftdParam->addLayer(nPetals, nSensors, isDoubleSided, sensorType, phalfangle, phi0, alpha, zposition, zoffset, signoffset,
                       supRinner, supThickness, supLengthMin, supLengthMax, supWidth, 0,
                       senRinner, senThickness, senLengthMin, senLengthMax, senWidth, 0);
  }
  m_gearMgr->setFTDParameters(ftdParam);

  return StatusCode::SUCCESS;
}

StatusCode GearSvc::convertSIT(dd4hep::DetElement& sit){
  dd4hep::rec::ZPlanarData* sitData = nullptr;
  try{
    sitData = sit.extension<dd4hep::rec::ZPlanarData>();
  }
  catch(std::runtime_error& e){
    warning() << e.what() << " " << sitData << endmsg;
    return StatusCode::FAILURE;
  }

  std::vector<dd4hep::rec::ZPlanarData::LayerLayout>& sitlayers = sitData->layers;
  int nLayers = sitlayers.size();
  double strip_angle_deg = sitData->angleStrip*rad_to_deg;
  
  gear::ZPlanarParametersImpl* sitParams = new gear::ZPlanarParametersImpl(1, 0.0, 0.0, 0.0, 0.0, 0.0);
  sitParams->setDoubleVal("strip_width_mm",  sitData->widthStrip*CLHEP::cm);
  sitParams->setDoubleVal("strip_length_mm", sitData->lengthStrip*CLHEP::cm);
  sitParams->setDoubleVal("strip_pitch_mm",  sitData->pitchStrip*CLHEP::cm);
  sitParams->setDoubleVal("strip_angle_deg", strip_angle_deg);
  std::vector<int> n_sensors_per_ladder;
  for( int layer=0; layer < nLayers; layer++){
    dd4hep::rec::ZPlanarData::LayerLayout& layout = sitlayers[layer];

    int nLadders = layout.ladderNumber;
    double phi0 = layout.phi0;
    double supRMin = layout.distanceSupport*CLHEP::cm;
    double supOffset = layout.offsetSupport*CLHEP::cm;
    double supThickness = layout.thicknessSupport*CLHEP::cm;
    double supHalfLength = layout.zHalfSupport*CLHEP::cm;
    double supWidth = layout.widthSupport*CLHEP::cm;
    double senRMin = layout.distanceSensitive*CLHEP::cm;
    double senOffset = layout.offsetSensitive*CLHEP::cm;
    double senThickness = layout.thicknessSensitive*CLHEP::cm;
    double senHalfLength = layout.zHalfSensitive*CLHEP::cm;
    double senWidth = layout.widthSensitive*CLHEP::cm;
    int nSensorsPerLadder = layout.sensorsPerLadder;
    double stripAngle = strip_angle_deg*CLHEP::degree;
    n_sensors_per_ladder.push_back(nSensorsPerLadder);
    sitParams->addLayer(nLadders, phi0, supRMin, supOffset, supThickness, supHalfLength, supWidth, 0, senRMin, senOffset, senThickness, senHalfLength, senWidth, 0);
  }
  sitParams->setIntVals("n_sensors_per_ladder",n_sensors_per_ladder);
  m_gearMgr->setSITParameters( sitParams ) ;

  return StatusCode::SUCCESS;
}

StatusCode GearSvc::convertTPC(dd4hep::DetElement& tpc){
  dd4hep::rec::FixedPadSizeTPCData* tpcData = nullptr;
  try{
    tpcData = tpc.extension<dd4hep::rec::FixedPadSizeTPCData>();
  }
  catch(std::runtime_error& e){
    warning() << e.what() << " " << tpcData << endmsg;
    return StatusCode::FAILURE;
  }

  gear::TPCParametersImpl *tpcParameters = new gear::TPCParametersImpl();
  gear::PadRowLayout2D *padLayout = new gear::FixedPadSizeDiskLayout(tpcData->rMinReadout*CLHEP::cm, tpcData->rMaxReadout*CLHEP::cm,
                                                                     tpcData->padHeight*CLHEP::cm, tpcData->padWidth*CLHEP::cm, tpcData->maxRow, 0.0);
  tpcParameters->setPadLayout(padLayout);
  tpcParameters->setMaxDriftLength(tpcData->driftLength*CLHEP::cm);
  tpcParameters->setDriftVelocity(    0.0); // SJA: not set in Mokka so set to 0.0                                                                                                          
  tpcParameters->setReadoutFrequency( 0.0);
  tpcParameters->setDoubleVal( "tpcOuterRadius" , tpcData->rMax*CLHEP::cm ) ;
  tpcParameters->setDoubleVal( "tpcInnerRadius",  tpcData->rMin*CLHEP::cm ) ;
  tpcParameters->setDoubleVal( "tpcInnerWallThickness",  tpcData->innerWallThickness*CLHEP::cm ) ;
  tpcParameters->setDoubleVal( "tpcOuterWallThickness",  tpcData->outerWallThickness*CLHEP::cm ) ;

  m_gearMgr->setTPCParameters(tpcParameters);

  return StatusCode::SUCCESS;
}

StatusCode GearSvc::convertDC(dd4hep::DetElement& dc){
  dd4hep::rec::FixedPadSizeTPCData* dcData = nullptr;
  try{
    dcData = dc.extension<dd4hep::rec::FixedPadSizeTPCData>();
  }
  catch(std::runtime_error& e){
    warning() << e.what() << " " << dcData << ", to search volume" << endmsg;
    // before extension ready, force to convert from volumes
    auto geomSvc = service<IGeomSvc>("GeomSvc");
    dd4hep::Readout readout = geomSvc->lcdd()->readout("DriftChamberHitsCollection");
    dd4hep::Segmentation seg = readout.segmentation();
    dd4hep::DDSegmentation::GridDriftChamber* grid = dynamic_cast< dd4hep::DDSegmentation::GridDriftChamber* > ( seg.segmentation() ) ;
    
    dcData = new dd4hep::rec::FixedPadSizeTPCData;
    dcData->rMinReadout = 99999;
    dcData->rMaxReadout = 0;
    std::vector<double> innerRadiusWalls,outerRadiusWalls;
    bool is_convert = true;
    dd4hep::Volume dc_vol = dc.volume();
    for(int i=0;i<dc_vol->GetNdaughters();i++){
      TGeoNode* daughter = dc_vol->GetNode(i);
      std::string nodeName = daughter->GetName();
      //info << nodeName << endmsg;
      if(nodeName.find("chamber_vol")!=-1||nodeName.find("assembly")!=-1){
	if(grid){
	  // if more than one chamber, just use the outer, TODO
	  dcData->rMinReadout = grid->DC_rbegin();
	  dcData->rMaxReadout = grid->DC_rend();
	  dcData->driftLength = grid->detectorLength();
	  dcData->padHeight   = grid->layer_width();
      dcData->padWidth    = dcData->padHeight;
	}
	else{
	  TGeoNode* next = daughter;
	  if(nodeName.find("assembly")!=-1){
	    // if more than one chamber, just use the outer, TODO 
	    next = daughter->GetDaughter(1);
	    std::string s = next->GetName();
	    if(s.find("chamber_vol")==-1){
	      error() << s << " not chamber_vol" << endmsg;
	      is_convert = false;
	    }
	  }
	  
	  //info() << next->GetName() << endmsg;
	  
	  const TGeoShape* chamber_shape = next->GetVolume()->GetShape();
	  if(chamber_shape->TestShapeBit(TGeoTube::kGeoTube)){
	    const TGeoTube* tube = (const TGeoTube*) chamber_shape;
	    double innerRadius = tube->GetRmin();
	    double outerRadius = tube->GetRmax();
	    double halfLength  = tube->GetDz();
	    dcData->driftLength = halfLength;
	  }
	  else{
	    error() << next->GetName() << " not TGeoTube::kGeoTube" << endmsg;
	    is_convert = false;
	  }
	  
	  dcData->maxRow = next->GetNdaughters();
	  if(dcData->maxRow>512){
	    error() << " layer number > 512, something wrong!" << endmsg;
	    is_convert = false;
	  }
	  for(int i=0;i<next->GetNdaughters();i++){
	    TGeoNode* layer = next->GetDaughter(i);
	    const TGeoShape* shape = layer->GetVolume()->GetShape();
	    if(shape->TestShapeBit(TGeoTube::kGeoTube)){
	      const TGeoTube* tube = (const TGeoTube*) shape;
	      double innerRadius = tube->GetRmin();
	      double outerRadius = tube->GetRmax();
	      double halfLength  = tube->GetDz();
	      if(innerRadius<dcData->rMinReadout) dcData->rMinReadout = innerRadius;
	      if(outerRadius>dcData->rMaxReadout) dcData->rMaxReadout = outerRadius;
	    }
	    else{
	      error() << layer->GetName() << " not TGeoTube::kGeoTube" << endmsg;
	      is_convert = false;
	    }
	  }
	  dcData->padHeight = (dcData->rMaxReadout-dcData->rMinReadout)/dcData->maxRow;
	  dcData->padWidth  = dcData->padHeight;
	}
      }
      else if(nodeName.find("wall_vol")!=-1){
	const TGeoShape* wall_shape = daughter->GetVolume()->GetShape();
        if(wall_shape->TestShapeBit(TGeoTube::kGeoTube)){
          const TGeoTube* tube = (const TGeoTube*) wall_shape;
          double innerRadius = tube->GetRmin();
          double outerRadius = tube->GetRmax();
          double halfLength  = tube->GetDz();
          innerRadiusWalls.push_back(innerRadius);
	  outerRadiusWalls.push_back(outerRadius);
        }
	else{
	  error() << nodeName << " not TGeoTube::kGeoTube" << endmsg;
	  is_convert = false;
	}
      }
    }
    if(innerRadiusWalls.size()<2||outerRadiusWalls.size()<2){
      error() << "wall number < 2" << endmsg;
      is_convert = false;
    }
    if(!is_convert){
      error() << "Cannot convert DC volume to extension data!" << endmsg;
      delete dcData;
      return StatusCode::FAILURE;
    }
    if(innerRadiusWalls[0]<innerRadiusWalls[innerRadiusWalls.size()-1]){
      dcData->rMin = innerRadiusWalls[0];
      dcData->rMax = outerRadiusWalls[outerRadiusWalls.size()-1];
      dcData->innerWallThickness = outerRadiusWalls[0]-innerRadiusWalls[0];
      dcData->outerWallThickness = outerRadiusWalls[outerRadiusWalls.size()-1]-innerRadiusWalls[innerRadiusWalls.size()-1];
    }
    else{
      dcData->rMin = innerRadiusWalls[innerRadiusWalls.size()-1];
      dcData->rMax = outerRadiusWalls[0];
      dcData->innerWallThickness = outerRadiusWalls[outerRadiusWalls.size()-1]-innerRadiusWalls[innerRadiusWalls.size()-1];
      dcData->outerWallThickness = outerRadiusWalls[0]-innerRadiusWalls[0];
    }
    info() << (*dcData) << endmsg;
  }

  // regard as TPCParameters, TODO: drift chamber parameters
  gear::TPCParametersImpl *tpcParameters = new gear::TPCParametersImpl();
  gear::PadRowLayout2D *padLayout = new gear::FixedPadSizeDiskLayout(dcData->rMinReadout*CLHEP::cm, dcData->rMaxReadout*CLHEP::cm,
                                                                     dcData->padHeight*CLHEP::cm, dcData->padWidth*CLHEP::cm, dcData->maxRow, 0.0);
  tpcParameters->setPadLayout(padLayout);
  tpcParameters->setMaxDriftLength(dcData->driftLength*CLHEP::cm);
  tpcParameters->setDriftVelocity(    0.0);
  tpcParameters->setReadoutFrequency( 0.0);
  tpcParameters->setDoubleVal( "tpcOuterRadius" , dcData->rMax*CLHEP::cm ) ;
  tpcParameters->setDoubleVal( "tpcInnerRadius",  dcData->rMin*CLHEP::cm ) ;
  tpcParameters->setDoubleVal( "tpcInnerWallThickness",  dcData->innerWallThickness*CLHEP::cm ) ;
  tpcParameters->setDoubleVal( "tpcOuterWallThickness",  dcData->outerWallThickness*CLHEP::cm ) ;

  m_gearMgr->setTPCParameters(tpcParameters);

  return StatusCode::SUCCESS;
}

StatusCode GearSvc::convertSET(dd4hep::DetElement& set){
  dd4hep::rec::ZPlanarData* setData = nullptr;
  try{
    setData = set.extension<dd4hep::rec::ZPlanarData>();
  }
  catch(std::runtime_error& e){
    warning() << e.what() << " " << setData << endmsg;
    return StatusCode::FAILURE;
  }

  std::vector<dd4hep::rec::ZPlanarData::LayerLayout>& setlayers = setData->layers;
  int nLayers = setlayers.size();
  double strip_angle_deg = setData->angleStrip*rad_to_deg;

  gear::ZPlanarParametersImpl* setParams = new gear::ZPlanarParametersImpl(1, 0.0 , 0.0 , 0.0 , 0.0 , 0.0);
  setParams->setDoubleVal("strip_width_mm",  setData->widthStrip*CLHEP::cm);
  setParams->setDoubleVal("strip_length_mm", setData->lengthStrip*CLHEP::cm);
  setParams->setDoubleVal("strip_pitch_mm",  setData->pitchStrip*CLHEP::cm);
  setParams->setDoubleVal("strip_angle_deg", strip_angle_deg);
  std::vector<int> n_sensors_per_ladder;
  for( int layer=0; layer < nLayers; layer++){
    dd4hep::rec::ZPlanarData::LayerLayout& layout = setlayers[layer];
    
    int nLadders = layout.ladderNumber;
    double phi0 = layout.phi0;
    double supRMin = layout.distanceSupport*CLHEP::cm;
    double supOffset = layout.offsetSupport*CLHEP::cm;
    double supThickness = layout.thicknessSupport*CLHEP::cm;
    double supHalfLength = layout.zHalfSupport*CLHEP::cm;
    double supWidth = layout.widthSupport*CLHEP::cm;
    double senRMin = layout.distanceSensitive*CLHEP::cm;
    double senOffset = layout.offsetSensitive*CLHEP::cm;
    double senThickness = layout.thicknessSensitive*CLHEP::cm;
    double senHalfLength = layout.zHalfSensitive*CLHEP::cm;
    double senWidth = layout.widthSensitive*CLHEP::cm;
    int nSensorsPerLadder = layout.sensorsPerLadder;
    double stripAngle = strip_angle_deg*CLHEP::degree;
    n_sensors_per_ladder.push_back(nSensorsPerLadder);
    setParams->addLayer(nLadders, phi0, supRMin, supOffset, supThickness, supHalfLength, supWidth, 0, senRMin, senOffset, senThickness, senHalfLength, senWidth, 0);
  }
  setParams->setIntVals("n_sensors_per_ladder",n_sensors_per_ladder);
  m_gearMgr->setSETParameters( setParams ) ;

  return StatusCode::SUCCESS;
}

TGeoNode* GearSvc::FindNode(TGeoNode* mother, char* name){
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
