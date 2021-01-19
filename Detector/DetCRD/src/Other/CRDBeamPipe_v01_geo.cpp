//====================================================================
// CepC BeamPipe models in DD4hep 
//--------------------------------------------------------------------
#include "OtherDetectorHelpers.h"
#include "TGeoScaledShape.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/detail/Handle.inl"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "XML/Utilities.h"
#include <cmath>
#include <map>
#include <string>

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::PlacedVolume;
using dd4hep::Ref_t;
using dd4hep::SensitiveDetector;
using dd4hep::Volume;

using dd4hep::rec::ConicalSupportData;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::VolCylinderImpl;
using dd4hep::rec::VolSurface;
using dd4hep::rec::volSurfaceList;

//BeamPipe
static Ref_t create_detector(Detector& theDetector,
			     xml_h element,
			     SensitiveDetector /*sens*/) {

  std::cout << "This is the Beampipe:"  << std::endl;

  //Access to the XML File
  xml_det_t x_beampipe = element;
  const std::string name = x_beampipe.nameStr();

  DetElement tube(  name, x_beampipe.id()  ) ;

  // --- create an envelope volume and position it into the world ---------------------
  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , tube ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, tube ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return tube ;
  
  //-----------------------------------------------------------------------------------
  ConicalSupportData* beampipeData = new ConicalSupportData ;

  const double phi0 = 0 ;
  const double dPhi = 360.0*dd4hep::degree;
  
  //Parameters we have to know about
  dd4hep::xml::Component xmlParameter = x_beampipe.child(_Unicode(parameter));
  const double crossingAngle  = xmlParameter.attr< double >(_Unicode(crossingangle));
  std::cout << "Crossing angle = " << crossingAngle << std::endl;

  for(xml_coll_t si( x_beampipe ,Unicode("section")); si; ++si) {
    xml_comp_t x_section(si);
    
    CEPC::ECrossType type = CEPC::getCrossType(x_section.attr< std::string >(_Unicode(type)));

    const double zstart       = x_section.attr< double > (_Unicode(zStart));
    const double zend         = x_section.attr< double > (_Unicode(zEnd));
    const double rInnerStart  = x_section.attr< double > (_Unicode(rStart));
    double rInnerEnd=0, size=0, shift=0;
    try{
      rInnerEnd = x_section.attr< double > (_Unicode(rEnd));
    }
    catch(std::runtime_error& e){
      rInnerEnd = rInnerStart;
    }
    if(type==CEPC::kWaist || type==CEPC::kFatWaist || type==CEPC::kCrotch || type==CEPC::kCrotchAsymUp || type==CEPC::kCrotchAsymDn){
      try{
	size = x_section.attr< double > (_Unicode(size));
      }
      catch(std::runtime_error& e){
	std::cout << "The maximum distance of runway is not set, will be calculated automatically by crossing angle" <<std::endl;
      }
      try{
        shift = x_section.attr< double > (_Unicode(shift));
      }
      catch(std::runtime_error& e){
	shift = 0;
      }
    }
    
    const std::string volName      = "BeamPipe_" + x_section.nameStr();

    std::cout << "section: "
	      << std::setw(8) << zstart      /dd4hep::mm
	      << std::setw(8) << zend	     /dd4hep::mm
	      << std::setw(8) << rInnerStart /dd4hep::mm
	      << std::setw(8) << rInnerEnd   /dd4hep::mm
	      << std::setw(8) << size        /dd4hep::mm
	      << std::setw(8) << type
	      << std::setw(35) << volName
	      << std::endl;    

    const double angle   = crossingAngle;
    const double zHalf   = fabs(zend - zstart) * 0.5;
    const double zCenter = fabs(zend + zstart) * 0.5;
    dd4hep::Material beamMaterial    = theDetector.material("beam");
    
    double clipSize = 100*dd4hep::mm;
    if(type==CEPC::kFlareLegUp || type==CEPC::kFlareLegDn){
      double total = rInnerStart, totalEnd = rInnerEnd;
      for(xml_coll_t li(x_section,_U(layer)); li; ++li)  {
	xml_comp_t  x_layer(li);
	const double thickness = x_layer.thickness(); 
	total += thickness;
	double thicknessEnd = 0;
	try{
	  thicknessEnd = x_layer.attr< double > (_Unicode(thicknessEnd));
	}
	catch(std::runtime_error& e){
	  thicknessEnd = thickness;
	}
	totalEnd += thicknessEnd;
      }
      clipSize = std::max(total, totalEnd)*tan(0.5*angle);
    }

    int ilayer = 0;
    double radius = rInnerStart;
    double radiusEnd = rInnerEnd;
    double pipeRadius = 0;
    double pipeThickness = 0;
    double pipeThicknessRel = 0;
    double pipeRadiusEnd = 0;
    double pipeThicknessEnd = 0;
    double pipeThicknessRelEnd = 0;
    for(xml_coll_t li(x_section,_U(layer)); li; ++li, ++ilayer)  {
      xml_comp_t  x_layer(li);
      double thickness = x_layer.thickness();
      dd4hep::Material material  = theDetector.material( x_layer.materialStr() );
      double thicknessEnd = 0;
      try{
	thicknessEnd = x_layer.attr< double > (_Unicode(thicknessEnd));
      }
      catch(std::runtime_error& e){
	thicknessEnd = thickness;
      }
      std::cout << "  layer: " << std::setw(8) << thickness/dd4hep::mm << std::setw(8) << thicknessEnd/dd4hep::mm << std::setw(15) << material.name() << std::endl;
      
      char suffix[20];
      sprintf(suffix,"_%d",ilayer);
      std::string layerName = volName + suffix;
      if(type==CEPC::kCenter || type==CEPC::kCenterSide){
	dd4hep::ConeSegment subLayer(zHalf, radius, radius+thickness, radiusEnd, radiusEnd+thicknessEnd, phi0, dPhi);
	dd4hep::Volume subLayerLog(volName, subLayer, material);
	dd4hep::Transform3D transformer(dd4hep::RotationY(0), dd4hep::Position(0, 0, zCenter));
	dd4hep::Transform3D transmirror(dd4hep::RotationY(180*dd4hep::degree), dd4hep::RotateY(dd4hep::Position(0, 0, zCenter), 180*dd4hep::degree));
	envelope.placeVolume(subLayerLog,  transformer);
	envelope.placeVolume(subLayerLog,  transmirror);
	std::cout << "fucd debug: radL = " << material.radLength()/dd4hep::mm << " intL = " << material.intLength()/dd4hep::mm << std::endl; 
	if(material.radLength()<10000*dd4hep::mm) subLayerLog.setVisAttributes(theDetector, "TubeVis");
	else                                      subLayerLog.setVisAttributes(theDetector, "VacVis");
	
        if(material.radLength()<10000*dd4hep::mm){
          double tEff    = thickness/material.radLength()*theDetector.material("G4_Be").radLength();
	  double tEffEnd = thicknessEnd/material.radLength()*theDetector.material("G4_Be").radLength();
          if(pipeRadius==0)    pipeRadius    = radius;
	  if(pipeRadiusEnd==0) pipeRadiusEnd = radiusEnd;
	  pipeThickness       += tEff;
	  pipeThicknessEnd    += tEffEnd;
	  pipeThicknessRel    += thickness;
	  pipeThicknessRelEnd += thicknessEnd;
	}
      }
      else if(type==CEPC::kLegs){
        double clipAngle = 0.5*angle;
        double lowNorml[3]  = { sin(clipAngle), 0, -cos(clipAngle)};
        double highNorml[3] = {-sin(clipAngle), 0,  cos(clipAngle)};
	dd4hep::CutTube subLayer(radius, radius+thickness, zHalf/cos(clipAngle), phi0, dPhi, lowNorml[0], lowNorml[1], lowNorml[2], highNorml[0], highNorml[1], highNorml[2]);
	dd4hep::Volume subLayerLog(volName, subLayer, material);
	dd4hep::Transform3D plusUpTransformer(dd4hep::RotationY(clipAngle), dd4hep::RotateY(dd4hep::Position(0,0,zCenter/cos(clipAngle)), clipAngle));
	dd4hep::Transform3D plusDownTransformer(dd4hep::RotationZYX(180*dd4hep::degree, -clipAngle, 0), dd4hep::RotateY(dd4hep::Position(0,0,zCenter/cos(clipAngle)), -clipAngle));
	dd4hep::Transform3D minusUpTransformer(dd4hep::RotationY(clipAngle), dd4hep::RotateY(dd4hep::Position(0,0,-zCenter/cos(clipAngle)), clipAngle));
	dd4hep::Transform3D minusDownTransformer(dd4hep::RotationZYX(180*dd4hep::degree, -clipAngle, 0), dd4hep::RotateY(dd4hep::Position(0,0,-zCenter/cos(clipAngle)), -clipAngle));
        envelope.placeVolume(subLayerLog, plusUpTransformer);
        envelope.placeVolume(subLayerLog, plusDownTransformer);
	envelope.placeVolume(subLayerLog, minusUpTransformer);
	envelope.placeVolume(subLayerLog, minusDownTransformer);

	if(material.radLength()<10000*dd4hep::mm) subLayerLog.setVisAttributes(theDetector, "TubeVis");
        else                                      subLayerLog.setVisAttributes(theDetector, "VacVis");
      }
      else if(type==CEPC::kFlareLegUp || type==CEPC::kFlareLegDn){
        double clipAngle = (type==CEPC::kFlareLegUp)?0.5*angle:-0.5*angle;
        double rOuter = radius+thickness;
        double rOuterEnd = radiusEnd+thicknessEnd;
	dd4hep::Tube clipSolid(0, zHalf*tan(0.5*angle)+clipSize/sin(0.5*angle), zHalf, phi0, dPhi);
	dd4hep::Transform3D clipTransformer(dd4hep::RotationY(-clipAngle), dd4hep::Position(0, 0, 0));
	dd4hep::Transform3D placementTransformer(dd4hep::RotationY(clipAngle), dd4hep::RotateY(dd4hep::Position(0, 0, zCenter/cos(clipAngle)), clipAngle));
	dd4hep::Transform3D placementTransmirror(dd4hep::RotationZYX(0, clipAngle, 180*dd4hep::degree), dd4hep::RotateY(dd4hep::Position(0, 0, -zCenter/cos(clipAngle)), -clipAngle));
	dd4hep::ConeSegment wholeSolid(zHalf + clipSize, radius, rOuter, radiusEnd, rOuterEnd, phi0, dPhi);
	dd4hep::IntersectionSolid layerSolid(wholeSolid, clipSolid, clipTransformer);
	dd4hep::Volume subLayerLog(volName, layerSolid, material);

        envelope.placeVolume(subLayerLog, placementTransformer);
        envelope.placeVolume(subLayerLog, placementTransmirror);
	
	if(material.radLength()<10000*dd4hep::mm) subLayerLog.setVisAttributes(theDetector, "TubeVis");
        else                                      subLayerLog.setVisAttributes(theDetector, "VacVis");
      }
      else if(type==CEPC::kCrotch){
        double beamAngle = 0.5*angle;
        if(size==0) size = (zstart*tan(beamAngle)+radius)*2;
	double x1 = 0.5*size - radius;
        double x2 = zend*tan(beamAngle);
        double y1 = radius + thickness;
        double y2 = y1;
        double axisAngle = atan((x2-x1)/zHalf/2);
        if(fabs(beamAngle-axisAngle)>1e-12){
	  std::cout << "Warning! axis angle not equal to beam angle. beam=" << beamAngle << " VS axis=" << axisAngle << ", user defined design and workable" << std::endl;
        }
        double zSide = 2*zHalf/cos(axisAngle)+y1*tan(axisAngle)+y2*tan(axisAngle);
        double xshift = 0.5*(x1+x2);
	dd4hep::Trd2 body1(x1, x2, y1, y2, zHalf);
	dd4hep::Trd2 cut1(x1+y1/cos(axisAngle), x2+y2/cos(axisAngle), y1, y2, zHalf);
	dd4hep::EllipticalTube side1(y1*cos(axisAngle), y1, 0.5*zSide);
	dd4hep::Transform3D unionTransformer1(dd4hep::RotationY(axisAngle), dd4hep::Position(xshift, 0, 0));
	dd4hep::Transform3D unionTransformer2(dd4hep::RotationY(-axisAngle), dd4hep::Position(-xshift, 0, 0));
	dd4hep::Transform3D sameTransformer(dd4hep::RotationY(0), dd4hep::Position(0, 0, 0));
	dd4hep::UnionSolid tmp1Solid(body1, side1, unionTransformer1);
	dd4hep::UnionSolid tmp2Solid(tmp1Solid, side1, unionTransformer2);
	dd4hep::IntersectionSolid shell(tmp2Solid, cut1, sameTransformer);
	dd4hep::Volume shellLog(volName+"Shell", shell, material);
	envelope.placeVolume(shellLog, dd4hep::Position(0, 0, zCenter));
	envelope.placeVolume(shellLog, dd4hep::Transform3D(dd4hep::RotationY(180*dd4hep::degree), dd4hep::Position(0, 0, -zCenter)));

	double yHole = y1-thickness;
	dd4hep::Trd2 body2(x1, x2, yHole, yHole, zHalf);
	dd4hep::Trd2 cut2(0, x2, yHole, yHole, zHalf);
	dd4hep::SubtractionSolid tmp3Solid(body2, cut2, sameTransformer);
	dd4hep::EllipticalTube side2(yHole*cos(axisAngle), yHole, zSide);
	dd4hep::UnionSolid tmp4Solid(tmp3Solid, side2, unionTransformer1);
	dd4hep::UnionSolid tmp5Solid(tmp4Solid, side2, unionTransformer2);
        double x1shift = radius-shift;
        double crotchAngle = atan(0.5*(x2-x1shift)/zHalf);
	dd4hep::EllipticalTube side3(yHole*cos(crotchAngle), yHole, zSide);
	dd4hep::Transform3D unionTransformer3(dd4hep::RotationY(crotchAngle), dd4hep::Position(0.5*(x2+x1shift), 0, 0));
	dd4hep::Transform3D unionTransformer4(dd4hep::RotationY(-crotchAngle), dd4hep::Position(-0.5*(x2+x1shift), 0, 0));
	dd4hep::UnionSolid tmp6Solid(tmp5Solid, side3, unionTransformer3);
	dd4hep::UnionSolid tmp7Solid(tmp6Solid, side3, unionTransformer4);
	dd4hep::IntersectionSolid vacuumPipe(tmp7Solid, cut1, sameTransformer);
	dd4hep::Volume pipeLog(volName+"Vacuum", vacuumPipe, beamMaterial);
        shellLog.placeVolume(pipeLog, dd4hep::Position(0, 0, 0));
	
	shellLog.setVisAttributes(theDetector, "TubeVis");
        pipeLog.setVisAttributes(theDetector, "VacVis");
      }
      else if(type==CEPC::kCrotchAsymUp || type==CEPC::kCrotchAsymDn){
	double beamAngle = 0.5*angle;
        double xC2 = (shift==0)?zend*tan(beamAngle):shift;
        if(radiusEnd==0) radiusEnd = radius;
	if(thicknessEnd==0) thicknessEnd = thickness;
        if(size==0) size = 2*radius;
        double rOuter = radius+thickness;
        double rOuterEnd = radiusEnd+thicknessEnd;
        double xMaxEnd = xC2+rOuterEnd;
        double yMax = 0.5*size+thickness;
	dd4hep::Trd2 body(0, xC2, yMax, rOuterEnd, zHalf);

        double expandAngle = atan(xC2/(2*zHalf));
        double edge1ToZAngle = atan((xMaxEnd-rOuter)/(2*zHalf));
        double edge2ToZAngle = atan((xC2-rOuterEnd+rOuter)/(2*zHalf));
        double edge2ToXAngle = 90*dd4hep::degree - edge2ToZAngle;
        double bottomAngle = 0.5*(180*dd4hep::degree-(edge2ToZAngle-edge1ToZAngle));
        double rotateAngle = 0.5*(edge1ToZAngle+edge2ToZAngle);
        double edge1ToCAngle = asin(sin(90*dd4hep::degree+edge1ToZAngle)/(xC2/sin(expandAngle))*(rOuter-rOuterEnd));
        double CToEConeAxisAngle = edge1ToCAngle-0.5*(edge2ToZAngle-edge1ToZAngle);
        if(fabs(rotateAngle-(expandAngle-CToEConeAxisAngle))>1e-12){
	  std::cout << "Warning! rotate angle was not calculated rightly. Please check input parameters whether satisfy the Waist case." << std::endl;
        }
	double a1 = rOuter/sin(bottomAngle)*sin(90*dd4hep::degree-edge1ToZAngle);
        double a2 = rOuterEnd/sin(180*dd4hep::degree-bottomAngle)*sin(90*dd4hep::degree-edge2ToZAngle);
        double zC1 = rOuter/sin(edge1ToCAngle)*sin(90*dd4hep::degree+edge1ToZAngle)*cos(CToEConeAxisAngle);
        double zC2 = rOuterEnd/rOuter*zC1;
        double zBottom = a1*tan(bottomAngle);
        double aC1 = a1/zBottom*zC1;
        double aC2 = a1/zBottom*zC2;
        double xC1InECone = zC1*tan(CToEConeAxisAngle);
        double xC2InECone = zC2*tan(CToEConeAxisAngle);
        double bC1 = sqrt(rOuter*rOuter/(1-xC1InECone*xC1InECone/aC1/aC1));
        double bC2 = sqrt(rOuterEnd*rOuterEnd/(1-xC2InECone*xC2InECone/aC2/aC2));
        double b1 = bC1/zC1*zBottom;
        if(fabs(bC1/zC1-bC2/zC2)>1e-12){
	  std::cout << "Warning! bC1/zC1 not equal to bC2/zC2. Please tell Chengdong(fucd@ihep.ac.cn)." << std::endl;
        }
        double pzTopCut = 0.5*(a1-a2)*tan(bottomAngle);
        double thetaCut1 = atan((0.5*(xC2+rOuterEnd)-0.5*rOuter)/(2*zHalf));
        double xcenterCut1 = 0.5*(0.5*(xC2+rOuterEnd)+0.5*rOuter);
	dd4hep::Trap cut1(zHalf, thetaCut1, 0, yMax, 0.5*rOuter, 0.5*rOuter, 0, rOuterEnd, 0.5*(xC2+rOuterEnd), 0.5*(xC2+rOuterEnd), 0);
	TGeoCone* pCone1 = new TGeoCone(pzTopCut, 0, a1, 0, a2); 
	//double factor = 
        TGeoScale* pScale1 = new TGeoScale(1, b1/a1, 1);
        TGeoScaledShape* pScaledShape1 = new TGeoScaledShape(pCone1, pScale1);
	dd4hep::Solid_type<TGeoScaledShape> side1(pScaledShape1);

        double xshift = 0.5*(xMaxEnd-a2*cos(rotateAngle)-rOuter+a1*cos(bottomAngle-edge2ToXAngle));
        double zshift = 0.5*(a2-a1)*sin(rotateAngle);
	dd4hep::Transform3D unionTransformer1(dd4hep::RotationY(rotateAngle), dd4hep::Position(xshift, 0, zshift));
	dd4hep::Transform3D cutTransformer1(dd4hep::RotationY(0), dd4hep::Position(xcenterCut1, 0, 0));
	dd4hep::UnionSolid tmp1Solid(body, side1, unionTransformer1);
	dd4hep::IntersectionSolid shell(tmp1Solid, cut1, cutTransformer1);
	dd4hep::Volume shellLog(volName, shell, material);
	if(type==CEPC::kCrotchAsymUp){
          envelope.placeVolume(shellLog, dd4hep::Position(0, 0, zCenter));
          envelope.placeVolume(shellLog, dd4hep::Transform3D(dd4hep::RotationX(180*dd4hep::degree), dd4hep::Position(0, 0, -zCenter)));
	}
        else{
          envelope.placeVolume(shellLog, dd4hep::Transform3D(dd4hep::RotationZ(180*dd4hep::degree), dd4hep::Position(0, 0, zCenter)));
	  envelope.placeVolume(shellLog, dd4hep::Transform3D(dd4hep::RotationY(180*dd4hep::degree), dd4hep::Position(0, 0,-zCenter)));
        }
	
        double edge1ToZ = atan((xMaxEnd-thicknessEnd-radius)/(2*zHalf));
        double edge2ToZ = atan((xC2-radiusEnd+radius)/(2*zHalf));
        double edge2ToX = 90*dd4hep::degree - edge2ToZ;
        double bottom = 0.5*(180*dd4hep::degree-(edge2ToZ-edge1ToZ));
        double rotate = 0.5*(edge1ToZ+edge2ToZ);
        double edge1ToC = asin(sin(90*dd4hep::degree+edge1ToZ)/(xC2/sin(expandAngle))*(radius-radiusEnd));
        double CToEConeAxis = edge1ToC-0.5*(edge2ToZ-edge1ToZ);
        if(fabs(rotate-(expandAngle-CToEConeAxis))>1e-12){
	  std::cout << "Warning! rotate angle was not calculated rightly. Please check input parameters whether satisfy the Waist case." << std::endl;
        }
	double a1Hole = radius/sin(bottom)*sin(90*dd4hep::degree-edge1ToZ);
        double a2Hole = radiusEnd/sin(180*dd4hep::degree-bottom)*sin(90*dd4hep::degree-edge2ToZ);
        double zC1Hole = radius/sin(edge1ToC)*sin(90*dd4hep::degree+edge1ToZ)*cos(CToEConeAxis);
        double zC2Hole = radiusEnd/radius*zC1Hole;
        double zBottomHole = a1Hole*tan(bottom);
        double aC1Hole = a1Hole/zBottomHole*zC1Hole;
        double aC2Hole = a1Hole/zBottomHole*zC2Hole;
        double xC1InEConeHole = zC1Hole*tan(CToEConeAxis);
        double xC2InEConeHole = zC2Hole*tan(CToEConeAxis);
        double bC1Hole = sqrt(radius*radius/(1-xC1InEConeHole*xC1InEConeHole/aC1Hole/aC1Hole));
        double bC2Hole = sqrt(radiusEnd*radiusEnd/(1-xC2InEConeHole*xC2InEConeHole/aC2Hole/aC2Hole));
        double b1Hole = bC1Hole/zC1Hole*zBottomHole;
        if(fabs(bC1Hole/zC1Hole-bC2Hole/zC2Hole)>1e-12){
	  std::cout << "Warning! bC1/zC1 not equal to bC2/zC2 for Hole. Please tell Chengdong(fucd@ihep.ac.cn)." << std::endl;
        }
	double pzTopCutHole = 0.5*(a1Hole-a2Hole)*tan(bottom);
	dd4hep::Trd2 body2(0, xC2, yMax-thickness, radiusEnd, zHalf);
        double thetaCut2 = atan((xC2-0.5*radius)/(2*zHalf));
        double xcenterCut2 = 0.5*radius+0.5*(xC2-0.5*radius);
	dd4hep::Trap cut2(zHalf, thetaCut2, 0, yMax-thickness, 0.5*radius, 0.5*radius, 0, radiusEnd, radiusEnd, radiusEnd, 0);
	TGeoCone* pCone2 = new TGeoCone(pzTopCutHole, 0, a1Hole, 0, a2Hole);
        TGeoScale* pScale2 = new TGeoScale(1, b1Hole/a1Hole, 1);
        TGeoScaledShape* pScaledShape2 = new TGeoScaledShape(pCone2, pScale2);
	dd4hep::Solid_type<TGeoScaledShape> side2(pScaledShape2);
        double xshiftHole = 0.5*(xMaxEnd-thicknessEnd-a2Hole*cos(rotate)-radius+a1Hole*cos(bottom-edge2ToX));
        double zshiftHole = 0.5*(a2Hole-a1Hole)*sin(rotate);
	dd4hep::Transform3D cutTransformer2(dd4hep::RotationY(rotate), dd4hep::Position(xshiftHole-xcenterCut2, 0, zshiftHole));
	dd4hep::IntersectionSolid vacuumPipe(cut2, side2, cutTransformer2);
	dd4hep::Volume pipeLog(volName, vacuumPipe, beamMaterial);
        shellLog.placeVolume(pipeLog, dd4hep::Position(xcenterCut2, 0, 0));
	
	shellLog.setVisAttributes(theDetector, "TubeVis");
	pipeLog.setVisAttributes(theDetector, "VacVis");
      }
      else if(type==CEPC::kWaist){
        double beamAngle = 0.5*angle;
        if(radiusEnd==0) radiusEnd = radius;
        if(size==0) size = (zend*tan(beamAngle)+radiusEnd)*2;
        if(thicknessEnd==0) thicknessEnd = thickness;
        double xC2 = 0.5*size - radiusEnd;
        double rOuter = radius+thickness;
        double rOuterEnd = radiusEnd+thicknessEnd;
        double rMaxEnd = 0.5*size+thicknessEnd;
	dd4hep::Trd2 body1(0, xC2, rOuter, rOuterEnd, zHalf);
	dd4hep::Trd2 cut(rOuter, rMaxEnd, rOuter, rOuterEnd, zHalf);

        double expandAngle = atan(xC2/(2*zHalf));
        double edge1ToZAngle = atan((rMaxEnd-rOuter)/(2*zHalf));
        double edge2ToZAngle = atan((xC2-rOuterEnd+rOuter)/(2*zHalf));
        double edge2ToXAngle = 90*dd4hep::degree - edge2ToZAngle;
        double bottomAngle = 0.5*(180*dd4hep::degree-(edge2ToZAngle-edge1ToZAngle));
        double rotateAngle = 0.5*(edge1ToZAngle+edge2ToZAngle);
        double edge1ToCAngle = asin(sin(90*dd4hep::degree+edge1ToZAngle)/(xC2/sin(expandAngle))*(rOuter-rOuterEnd));
        double CToEConeAxisAngle = edge1ToCAngle-0.5*(edge2ToZAngle-edge1ToZAngle);
        if(fabs(rotateAngle-(expandAngle-CToEConeAxisAngle))>1e-12){
	  std::cout << "Warning! rotate angle was not calculated rightly. Please check input parameters whether satisfy the Waist case." << std::endl;
        }
	double a1 = rOuter/sin(bottomAngle)*sin(90*dd4hep::degree-edge1ToZAngle);
        double a2 = rOuterEnd/sin(180*dd4hep::degree-bottomAngle)*sin(90*dd4hep::degree-edge2ToZAngle);
        double zC1 = rOuter/sin(edge1ToCAngle)*sin(90*dd4hep::degree+edge1ToZAngle)*cos(CToEConeAxisAngle);
        double zC2 = rOuterEnd/rOuter*zC1;
        double zBottom = a1*tan(bottomAngle);
        double aC1 = a1/zBottom*zC1;
        double aC2 = a1/zBottom*zC2;
        double xC1InECone = zC1*tan(CToEConeAxisAngle);
        double xC2InECone = zC2*tan(CToEConeAxisAngle);
        double bC1 = sqrt(rOuter*rOuter/(1-xC1InECone*xC1InECone/aC1/aC1));
        double bC2 = sqrt(rOuterEnd*rOuterEnd/(1-xC2InECone*xC2InECone/aC2/aC2));
        double b1 = bC1/zC1*zBottom;
        if(fabs(bC1/zC1-bC2/zC2)>1e-12){
	  std::cout << "Warning! bC1/zC1 not equal to bC2/zC2." << std::endl;
        }
        double pzTopCut = 0.5*(a1-a2)*tan(bottomAngle);
	TGeoCone* pCone1 = new TGeoCone(pzTopCut, 0, a1, 0, a2);
	TGeoScale* pScale1 = new TGeoScale(1, b1/a1, 1);
	TGeoScaledShape* pScaledShape1 = new TGeoScaledShape(pCone1,pScale1);
	dd4hep::Solid_type<TGeoScaledShape> side1(pScaledShape1);

	double xshift = 0.5*(rMaxEnd-a2*cos(rotateAngle)-rOuter+a1*cos(bottomAngle-edge2ToXAngle));
        double zshift = 0.5*(a2-a1)*sin(rotateAngle);
	dd4hep::Transform3D unionTransformer1(dd4hep::RotationY(rotateAngle), dd4hep::Position(xshift, 0, zshift));
	dd4hep::Transform3D unionTransformer2(dd4hep::RotationY(-rotateAngle), dd4hep::Position(-xshift, 0, zshift));
	dd4hep::Transform3D sameTransformer(dd4hep::RotationY(0), dd4hep::Position(0, 0, 0));
	dd4hep::UnionSolid tmp1Solid(body1, side1, unionTransformer1);
	dd4hep::UnionSolid tmp2Solid(tmp1Solid, side1, unionTransformer2);
	dd4hep::IntersectionSolid shell(tmp2Solid, cut, sameTransformer);
	dd4hep::Volume shellLog(volName+"Shell", shell, material);
        envelope.placeVolume(shellLog, dd4hep::Position(0, 0, zCenter));
	envelope.placeVolume(shellLog, dd4hep::Transform3D(dd4hep::RotationY(180*dd4hep::degree), dd4hep::Position(0, 0, -zCenter)));
	
        double edge1ToZ = atan((0.5*size-radius)/(2*zHalf));
        double edge2ToZ = atan((xC2-radiusEnd+radius)/(2*zHalf));
        double edge2ToX = 90*dd4hep::degree - edge2ToZ;
        double bottom = 0.5*(180*dd4hep::degree-(edge2ToZ-edge1ToZ));
        double rotate = 0.5*(edge1ToZ+edge2ToZ);
        double edge1ToC = asin(sin(90*dd4hep::degree+edge1ToZ)/(xC2/sin(expandAngle))*(radius-radiusEnd));
        double CToEConeAxis = edge1ToC-0.5*(edge2ToZ-edge1ToZ);
        if(fabs(rotate-(expandAngle-CToEConeAxis))>1e-12){
	  std::cout << "Warning! rotate angle was not calculated rightly. Please check input parameters whether satisfy the Waist case." << std::endl;
        }
	double a1Hole = radius/sin(bottom)*sin(90*dd4hep::degree-edge1ToZ);
        double a2Hole = radiusEnd/sin(180*dd4hep::degree-bottom)*sin(90*dd4hep::degree-edge2ToZ);
        double zC1Hole = radius/sin(edge1ToC)*sin(90*dd4hep::degree+edge1ToZ)*cos(CToEConeAxis);
        double zC2Hole = radiusEnd/radius*zC1Hole;
        double zBottomHole = a1Hole*tan(bottom);
        double aC1Hole = a1Hole/zBottomHole*zC1Hole;
        double aC2Hole = a1Hole/zBottomHole*zC2Hole;
        double xC1InEConeHole = zC1Hole*tan(CToEConeAxis);
        double xC2InEConeHole = zC2Hole*tan(CToEConeAxis);
        double bC1Hole = sqrt(radius*radius/(1-xC1InEConeHole*xC1InEConeHole/aC1Hole/aC1Hole));
        double bC2Hole = sqrt(radiusEnd*radiusEnd/(1-xC2InEConeHole*xC2InEConeHole/aC2Hole/aC2Hole));
        double b1Hole = bC1Hole/zC1Hole*zBottomHole;
        if(fabs(bC1Hole/zC1Hole-bC2Hole/zC2Hole)>1e-12){
	  std::cout << "Warning! bC1/zC1 not equal to bC2/zC2 for Hole." << std::endl;
        }
        double pzTopCutHole = 0.5*(a1Hole-a2Hole)*tan(bottom);
	dd4hep::Trd2 body2(0, xC2, radius, radiusEnd, zHalf);
	dd4hep::Trd2 cut2(radius, 0.5*size, radius, radiusEnd, zHalf);
	TGeoCone* pCone2 = new TGeoCone(pzTopCutHole, 0, a1Hole, 0, a2Hole);
        TGeoScale* pScale2 = new TGeoScale(1, b1Hole/a1Hole, 1);
        TGeoScaledShape* pScaledShape2 = new TGeoScaledShape(pCone2,pScale2);
	dd4hep::Solid_type<TGeoScaledShape> side2(pScaledShape2);

        double xshiftHole = 0.5*(0.5*size-a2Hole*cos(rotate)-radius+a1Hole*cos(bottom-edge2ToX));
        double zshiftHole = 0.5*(a2Hole-a1Hole)*sin(rotate);
	dd4hep::Transform3D unionTransformer3(dd4hep::RotationY(rotate), dd4hep::Position(xshiftHole, 0, zshiftHole));
	dd4hep::Transform3D unionTransformer4(dd4hep::RotationY(-rotate), dd4hep::Position(-xshiftHole, 0, zshiftHole));
	dd4hep::UnionSolid tmp3Solid(body2, side2, unionTransformer3);
	dd4hep::UnionSolid tmp4Solid(tmp3Solid, side2, unionTransformer4);
	dd4hep::IntersectionSolid vacuumPipe(tmp4Solid, cut, sameTransformer);
	dd4hep::Volume pipeLog(volName+"Vacuum", vacuumPipe, beamMaterial);
        shellLog.placeVolume(pipeLog, dd4hep::Position(0, 0, 0));
	
	shellLog.setVisAttributes(theDetector, "TubeVis");
	pipeLog.setVisAttributes(theDetector, "VacVis");
      }
      else if(type == CEPC::kFatWaist){
        double beamAngle = 0.5*angle;
        if(radiusEnd==0) radiusEnd = radius;
        if(size==0) size = (zend*tan(beamAngle)+radiusEnd)*2;
        if(thicknessEnd==0) thicknessEnd = thickness;
        double rOuter = radius+thickness;
        double rOuterEnd = radiusEnd+thicknessEnd;
        double yMaxEnd = 0.5*size+thicknessEnd;
	dd4hep::Transform3D sameTransformer(dd4hep::RotationY(0), dd4hep::Position(0, 0, 0));

	dd4hep::Trd2 body1(rOuter, rOuterEnd, rOuter, yMaxEnd, zHalf);
	dd4hep::ConeSegment cone1(zHalf, 0, rOuter, 0, rOuterEnd, phi0, dPhi);
	dd4hep::IntersectionSolid shell(cone1, body1, sameTransformer);
	dd4hep::Volume shellLog(volName, shell, material);
        envelope.placeVolume(shellLog, dd4hep::Position(0, 0, zCenter));
	envelope.placeVolume(shellLog, dd4hep::Transform3D(dd4hep::RotationY(180*dd4hep::degree), dd4hep::Position(0, 0, -zCenter)));

	dd4hep::Trd2 body2(radius, radiusEnd, radius, 0.5*size, zHalf);
	dd4hep::ConeSegment cone2(zHalf, 0, radius, 0, radius, phi0, dPhi);
	dd4hep::IntersectionSolid vacuumPipe(cone2, body2, sameTransformer);
	dd4hep::Volume pipeLog(volName, vacuumPipe, beamMaterial);
        shellLog.placeVolume(pipeLog, dd4hep::Position(0, 0, 0));

	shellLog.setVisAttributes(theDetector, "TubeVis");
	pipeLog.setVisAttributes(theDetector, "VacVis");
      }
      radius += thickness;
      radiusEnd += thicknessEnd;
    }
    if( type == CEPC::kCenter ) { // store only the central sections !
      ConicalSupportData::Section section ;
      section.rInner = pipeRadius + 0.5*(pipeThicknessRel-pipeThickness) ;
      section.rOuter = section.rInner + pipeThickness;
      section.zPos   = zstart ;

      ConicalSupportData::Section sectionEnd ;
      sectionEnd.rInner = pipeRadiusEnd + 0.5*(pipeThicknessRelEnd-pipeThicknessEnd) ;
      sectionEnd.rOuter = sectionEnd.rInner + pipeThicknessEnd;
      sectionEnd.zPos   = zend ;

      if(beampipeData->sections.size()!=0){
	ConicalSupportData::Section last = beampipeData->sections.back(); 
	if(last.rInner != section.rInner || last.rOuter != section.rOuter){
	  section.zPos = zstart + 1e-9*dd4hep::mm ;
	  beampipeData->sections.push_back( section );
	}
      }
      else beampipeData->sections.push_back( section );
      beampipeData->sections.push_back( sectionEnd ) ;
    }
  }//for all xmlSections
  
  // add a surface just inside the beampipe for tracking:
  double rInner = beampipeData->sections[0].rInner;
  double rOuter = beampipeData->sections[0].rOuter;
  Vector3D oCyl( 0.5*(rInner+rOuter)  , 0. , 0.  ) ;
  VolCylinder pipeSurf( envelope , SurfaceType( SurfaceType::Helper ) ,
			0.5*(rOuter-rInner) , 0.5*(rOuter-rInner), oCyl ) ;
  volSurfaceList( tube )->push_back( pipeSurf ) ;
  
  tube.addExtension< ConicalSupportData >( beampipeData ) ;

  //--------------------------------------
  tube.setVisAttributes( theDetector, x_beampipe.visStr(), envelope );
  
  //debug
  std::cout << "============ConicalSupportData============" << std::endl;
  for(unsigned int i=0;i<beampipeData->sections.size();i++){
    std::cout << std::setw(8) << beampipeData->sections[i].zPos    /dd4hep::mm
              << std::setw(8) << beampipeData->sections[i].rInner  /dd4hep::mm
	      << std::setw(8) << beampipeData->sections[i].rOuter  /dd4hep::mm
	      << std::endl;
  }

  return tube;
}
DECLARE_DETELEMENT(DD4hep_CRDBeamPipe_v01, create_detector)

DD4HEP_INSTANTIATE_SHAPE_HANDLE(TGeoScaledShape);
