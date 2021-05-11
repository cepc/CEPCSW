/////////////////////////////////////////////////////////////////////////
//  An implementation of genfit AbsMaterialInterface
//
//  Author:
//   Yao ZHANG (zhangyao@ihep.ac.cn)
//
//
/////////////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_GENFITMATERIALINTERFACE_H
#define RECGENFITALG_GENFITMATERIALINTERFACE_H

#include "AbsMaterialInterface.h"

class RKTrackRep;
class TGeoManager;
class TGeoNode;
namespace dd4hep{
    class Detector;
}

/**
 * @brief MaterialInterface implementation for use with ROOT's TGeoManager.
 */
class GenfitMaterialInterface : public genfit::AbsMaterialInterface{
    public:
        GenfitMaterialInterface(const dd4hep::Detector* dd4hepGeo);
        virtual ~GenfitMaterialInterface(){};
        static GenfitMaterialInterface* getInstance(
                const dd4hep::Detector* dd4hepGeo);
        void destruct();

        //Set Detector pointer of dd4hep
        void setDetector(dd4hep::Detector*);

        /** @brief Initialize the navigator at given position and with given
          direction.  Returns true if the volume changed.
          */
        bool initTrack(double posX, double posY, double posZ,
                double dirX, double dirY, double dirZ) override;

        genfit::Material getMaterialParameters() override;


        /** @brief Make a step (following the curvature) until step length
         * sMax or the next boundary is reached.  After making a step to a
         * boundary, the position has to be beyond the boundary, i.e. the
         * current material has to be that beyond the boundary.  The actual
         * step made is returned.
         */
        double findNextBoundary(const genfit::RKTrackRep* rep,
                const genfit::M1x7& state7,
                double sMax,
                bool varField = true) override;

        // ClassDef(GenfitMaterialInterface, 1);

    private:
        static GenfitMaterialInterface* m_instance;
        TGeoManager* getGeoManager();
        TGeoManager* m_geoManager;
        double getSafeDistance();
        double getStep();
        TGeoNode* findNextBoundary(double stepmax, const char* path="",
                bool frombdr=false);
        bool isSameLocation(double posX, double posY, double posZ,
                bool change=false);
        void setCurrentDirection(double nx, double ny, double nz);

};

/** @} */

#endif // GenfitMaterialInterface
