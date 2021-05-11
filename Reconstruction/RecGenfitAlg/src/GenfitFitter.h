//////////////////////////////////////////////////////////////////
///
/// This is an interface of to call genfit fitter
///
/// In this file, including:
///   a genfit fitter class
///   Fitter type is fixed after creation
///   Magnetic field and geometry for material effect are singletons in genfit.
///   They will be set before fitting and before each time call track rep.
///
/// Authors:
///   Yao ZHANG (zhangyao@ihep.ac.cn)
///
//////////////////////////////////////////////////////////////////

#ifndef RECGENFITALG_GENFITFITTER_H
#define RECGENFITALG_GENFITFITTER_H

#include "AbsKalmanFitter.h"
#include <string>

class GenfitField;
class GenfitMaterialInterface;
class GenfitTrack;
class TVector3;
class TGeoManager;
class ISvcLocator;
namespace genfit{
    class AbsKalmanFitter;
    class KalmanFitterRefTrack;
    class DAF;
}
namespace dd4hep{
    class OverlayedField;
    class Detector;
}


/// The GenfitTrack class
class GenfitFitter{
    public:
        /// Type of fitters are :DAFRef,DAF,KalmanFitter,KalmanFitterRefTrack
        GenfitFitter(const char* type="DAFRef",const char* name="GenfitFitter");
        virtual ~GenfitFitter();

        /// Magnetic field and geometry for material effect in genfit
        /// please SET before use !!!!
        void setField(const GenfitField* field);
        /// please SET before use !!!!
        void setGeoMaterial(const dd4hep::Detector* dd4hepGeo);

        /// Main fitting function
        int processTrack(GenfitTrack* track, bool resort=false);

        /// fitting with rep
        int processTrackWithRep(GenfitTrack* track,int repID=0,
                bool resort=false);

        /// Extrapolate the track to the CDC hit
        /// Output: poca pos and dir and poca distance to the hit wire
        /// Input: genfit track, pos and mom, two ends of a wire
        ///        pos, and mom are position & momentum at starting point
        double extrapolateToHit(TVector3& poca, TVector3& pocaDir,
                TVector3& pocaOnWire, double& doca, const GenfitTrack* track,
                TVector3 pos, TVector3 mom, TVector3 end0, TVector3 end1,
                int debug=0, int repID=0, bool stopAtBoundary=false,
                bool calcJacobianNoise=false);

        /// Extrapolate the track to the cyliner at fixed raidus
        /// Output: pos and mom at fixed radius
        /// Input: genfitTrack, radius of cylinder at center of the origin,
        ///        repID, stopAtBoundary and calcAverageState
        double extrapolateToCylinder(TVector3& pos, TVector3& mom,
                GenfitTrack* track, double radius, const TVector3 linePoint,
                const TVector3 lineDirection, int hitID =0, int repID=0,
                bool stopAtBoundary=false, bool calcJacobianNoise=false);

        /// Extrapolate the track to the point
        /// Output: pos and mom of POCA point to point
        /// Input: genfitTrack,point,repID,stopAtBoundary and calcAverageState
        /// repID same with pidType
        double extrapolateToPoint(TVector3& pos, TVector3& mom,
                const GenfitTrack* genfitTrack, const TVector3& point,
                int repID=0, bool stopAtBoundary = false,
                bool calcJacobianNoise = false) const;

        /// setters of fitter properties
        void setFitterType(const char* val);
        void setMinIterations(unsigned int val);
        void setMaxIterations(unsigned int val);
        void setMaxIterationsBetas(double bStart,double bFinal,unsigned int val);
        void setDeltaPval(double val);
        void setRelChi2Change(double val);
        void setBlowUpFactor(double val);
        void setResetOffDiagonals(bool val);
        void setBlowUpMaxVal(double val);
        void setMultipleMeasurementHandling(
                genfit::eMultipleMeasurementHandling val);
        void setMaxFailedHits(int val);
        void setConvergenceDeltaWeight(double val);
        void setAnnealingScheme(double bStart,double bFinal,unsigned int nSteps);
        //TODO chi2cut?
        void setNoEffects(bool val);
        void setEnergyLossBetheBloch(bool val);
        void setNoiseBetheBloch(bool val);
        void setNoiseCoulomb(bool val);
        void setEnergyLossBrems(bool val);
        void setNoiseBrems(bool val);
        void setIgnoreBoundariesBetweenEqualMaterials(bool val);
        void setMscModelName(std::string val);
        void setMaterialDebugLvl(unsigned int val);
        void setDebug(unsigned int val);
        void setHist(unsigned int val);

        /// getters of fitter properties
        std::string getFitterType() const {return m_fitterType;}
        unsigned int getMinIterations() const { return m_minIterations; }
        unsigned int getMaxIterations() const { return m_maxIterations; }
        double getDeltaPval() const { return m_deltaPval; }
        double getRelChi2Change() const { return m_relChi2Change; }
        double getBlowUpFactor() const { return m_blowUpFactor; }
        double getBlowUpMaxVal() const { return m_blowUpMaxVal; }
        bool getResetOffDiagonals() const { return m_resetOffDiagonals; }
        genfit::eMultipleMeasurementHandling getMultipleMeasurementHandling()
            const { return m_multipleMeasurementHandling; }
        int getMaxFailedHits() const { return m_maxFailedHits; }
        double getConvergenceDeltaWeight() const { return m_deltaWeight; }
        float getAnnealingBetaStart() const {return m_annealingBetaStart;}
        float getAnnealingBetaStop() const {return m_annealingBetaStop;}
        bool getNoEffects(){return m_noEffects;}
        bool getEnergyLossBetheBloch(){return m_energyLossBetheBloch;}
        bool getNoiseBetheBloch(){return m_noiseBetheBloch;}
        bool getNoiseCoulomb(){return m_noiseCoulomb;}
        bool getEnergyLossBrems(){return m_energyLossBrems;}
        bool getNoiseBrems(){return m_noiseBrems;}
        bool getIgnoreBoundariesBetweenEqualMaterials()
        {return m_ignoreBoundariesBetweenEqualMaterials;}
        std::string getMscModelName(){return m_mscModelName;}
        int  getDebug() const {return m_debug;}
        int  getHist() const {return m_hist;}

        /// Printer
        void print(const char* name="");
        void initHist(std::string name);
        void writeHist();

        const char* getName(void) const {return m_name;}

    private:
        GenfitFitter(const GenfitFitter&){};
        GenfitFitter& operator=(const GenfitFitter&);

        /// Initialze fitter and setting fitter parameter
        int init(bool deleteOldFitter=false);

        /// Get DAF fitter
        genfit::DAF* getDAF();
        /// Get KalmanFitterRefTrack
        genfit::KalmanFitterRefTrack* getKalRef();

        genfit::AbsKalmanFitter* m_absKalman;/// kalman fitter object

        const GenfitField* m_genfitField;//pointer to genfit field
        GenfitMaterialInterface* m_geoMaterial;//pointer to genfit geo

        ///fitting method: DAFRef,DAF,KalmanFitter,KalmanFitterRefTrack
        std::string m_fitterType;

        const char* m_name;

        /// control parameters of fitter
        unsigned int m_minIterations;     /// minimum number of iterations
        unsigned int m_maxIterations;     /// maximum number of iterations
        double       m_deltaPval;         /// delta pVal that converged
        double       m_relChi2Change;     ///
        double       m_blowUpFactor;      ///
        bool         m_resetOffDiagonals; ///
        double       m_blowUpMaxVal;      ///
        genfit::eMultipleMeasurementHandling m_multipleMeasurementHandling;///
        int          m_maxFailedHits;     ///
        double       m_deltaWeight;       /// delta weight that converged
        double       m_annealingBetaStart;/// start beta for annealing
        double       m_annealingBetaStop; /// stop beta for annealing
        unsigned int m_annealingNSteps;   /// n steps

        /// control parAmeters of material effects
        bool         m_noEffects;
        bool         m_energyLossBetheBloch;
        bool         m_noiseBetheBloch;
        bool         m_noiseCoulomb;
        bool         m_energyLossBrems;
        bool         m_noiseBrems;
        bool         m_ignoreBoundariesBetweenEqualMaterials;
        std::string  m_mscModelName;

        int          m_debug;             /// debug
        int          m_hist;             /// hist
};

#endif
