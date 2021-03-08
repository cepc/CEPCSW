#ifndef PanUTIL
#define PanUTIL 1

#include <sstream>
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"
namespace PanUtil{
    std::string Convert (float number);
    dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0) ;
    void line_a_b(float x1, float y1, float x2, float y2, float& a, float& b);
    float getPhi(float x, float y);
    int partition(float x, float y);
    int getLayer(float x, float y, float z, std::vector<float>& layers);
    int getLayer_v1(float x, float y, float z, std::vector<float>& layers);
    int getStave(float x, float y);
    double getFieldFromCompact();
    std::vector<double> getTrackingRegionExtent();
}
#endif
