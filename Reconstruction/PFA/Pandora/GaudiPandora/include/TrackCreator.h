/**
 *  @brief  Header file for the track creator class.
 * 
 *  $Log: $
 */

#ifndef TRACK_CREATOR_H
#define TRACK_CREATOR_H 1




#include "GaudiKernel/ISvcLocator.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackConst.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/ReconstructedParticleConst.h"

#include "Api/PandoraApi.h"
#include "Objects/Helix.h"

namespace gear { class GearMgr; }

class CollectionMaps;

typedef std::vector<const edm4hep::Track *> TrackVector;
typedef std::set<unsigned int> TrackList;
typedef std::map<edm4hep::ConstTrack, int> TrackToPidMap;
/*
inline LCCollectionVec *newTrkCol(const std::string &name, LCEvent *evt , bool isSubset)
{
    LCCollectionVec* col = new LCCollectionVec( LCIO::TRACK ) ;

    LCFlagImpl hitFlag(0) ;
    hitFlag.setBit( LCIO::TRBIT_HITS ) ;
    col->setFlag( hitFlag.getFlag()  ) ;
    evt->addCollection( col , name ) ;
    col->setSubset( isSubset ) ;

    return col ;
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrackCreator class
 */
class TrackCreator
{
public:
    typedef std::vector<double> DoubleVector;
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        StringVector    m_trackCollections;                     ///< The reconstructed track collections
        StringVector    m_kinkVertexCollections;                ///< The kink vertex collections
        StringVector    m_prongVertexCollections;               ///< The prong vertex collections
        StringVector    m_splitVertexCollections;               ///< The split vertex collections
        StringVector    m_v0VertexCollections;                  ///< The v0 vertex collections

        StringVector    m_prongSplitVertexCollections;          ///< Concatenated list of prong and split vertex collections
        int             m_shouldFormTrackRelationships;         ///< Whether to form pandora track relationships using v0 and kink info

        int             m_minTrackHits;                         ///< Track quality cut: the minimum number of track hits
        int             m_minFtdTrackHits;                      ///< Track quality cut: the minimum number of FTD track hits for FTD only tracks
        int             m_maxTrackHits;                         ///< Track quality cut: the maximum number of track hits

        float           m_d0TrackCut;                           ///< Track d0 cut used to determine whether track can be used to form pfo
        float           m_z0TrackCut;                           ///< Track z0 cut used to determine whether track can be used to form pfo

        int             m_usingNonVertexTracks;                 ///< Whether can form pfos from tracks that don't start at vertex
        int             m_usingUnmatchedNonVertexTracks;        ///< Whether can form pfos from unmatched tracks that don't start at vertex
        int             m_usingUnmatchedVertexTracks;           ///< Whether can form pfos from unmatched tracks that start at vertex
        float           m_unmatchedVertexTrackMaxEnergy;        ///< Maximum energy for unmatched vertex track

        float           m_d0UnmatchedVertexTrackCut;            ///< d0 cut used to determine whether unmatched vertex track can form pfo
        float           m_z0UnmatchedVertexTrackCut;            ///< z0 cut used to determine whether unmatched vertex track can form pfo
        float           m_zCutForNonVertexTracks;               ///< Non vtx track z cut to determine whether track can be used to form pfo

        int             m_reachesECalNTpcHits;                  ///< Minimum number of tpc hits to consider track as reaching ecal
        int             m_reachesECalNFtdHits;                  ///< Minimum number of ftd hits to consider track as reaching ecal
        float           m_reachesECalTpcOuterDistance;          ///< Max distance from track to tpc r max to id whether track reaches ecal
        int             m_reachesECalMinFtdLayer;               ///< Min layer in Ftd for tracks to be considered to have reached decal
        float           m_reachesECalTpcZMaxDistance;           ///< Max distance from track to tpc z max to id whether track reaches ecal
        float           m_reachesECalFtdZMaxDistance;           ///< Max distance from track hit to ftd z position to identify ftd hits
        float           m_curvatureToMomentumFactor;            ///< Constant relating track curvature in b field to momentum

        float           m_minTrackECalDistanceFromIp;           ///< Sanity check on separation between ip and track projected ecal position
        float           m_maxTrackSigmaPOverP;                  ///< Track fraction momentum error cut
        float           m_minMomentumForTrackHitChecks;         ///< Min track momentum required to perform final quality checks on number of hits

        float           m_tpcMembraneMaxZ;                      ///< Tpc membrane max z coordinate
        float           m_maxTpcInnerRDistance;                 ///< Track cut on distance from tpc inner r to id whether track can form pfo
        float           m_minTpcHitFractionOfExpected;          ///< Minimum fraction of TPC hits compared to expected
        int             m_minFtdHitsForTpcHitFraction;          ///< Minimum number of FTD hits to ignore TPC hit fraction
        bool            m_use_dd4hep_geo;          ///
        bool            m_debug;
    };

    /**
     *  @brief  Constructor
     * 
     */
     TrackCreator(const Settings &settings, const pandora::Pandora *const pPandora, ISvcLocator* svcloc);

    /**
     *  @brief  Destructor
     */
     ~TrackCreator();

    /**
     *  @brief  Create associations between tracks, V0s, kinks, etc
     * 
     */
    pandora::StatusCode CreateTrackAssociations(const CollectionMaps& collectionMaps);

    const edm4hep::Track* GetTrackAddress(const CollectionMaps& collectionMaps, const edm4hep::ConstTrack& pTrack );
    /**
     *  @brief  Create tracks, insert user code here
     * 
     */
    pandora::StatusCode CreateTracks(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Get the track vector
     * 
     *  @return The track vector
     */
    const TrackVector &GetTrackVector() const;

    /**
     *  @brief  Reset the track creator
     */
    void Reset();

private:
    /**
     *  @brief  Extract kink information from specified lcio collections
     * 
     */
    pandora::StatusCode ExtractKinks(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Extract prong and split information from specified lcio collections
     * 
     *  @param  pLCEvent the lcio event
     */
    //pandora::StatusCode ExtractProngsAndSplits(const EVENT::LCEvent *const pLCEvent);
    pandora::StatusCode ExtractProngsAndSplits(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Extract v0 information from specified lcio collections
     * 
     *  @param  pLCEvent the lcio event
     */
   //pandora::StatusCode ExtractV0s(const EVENT::LCEvent *const pLCEvent);
    pandora::StatusCode ExtractV0s(const CollectionMaps& collectionMaps);

    /**
     *  @brief  Whether the track vertex conflicts with previously provided relationship information
     * 
     *  @param  trackVec the vector of tracks associated with the vertex
     */
    bool IsConflictingRelationship(const edm4hep::ConstReconstructedParticle &Particle) const;

    /**
     *  @brief  Whether a track is a v0 track
     * 
     * 
     */
    bool IsV0(unsigned int pTrack_id) const;

    /**
     *  @brief  Whether a track is a parent track
     * 
     * 
     */
    bool IsParent(unsigned int pTrack_id) const;

    /**
     *  @brief  Whether a track is a daughter track
     * 
     */
    bool IsDaughter(unsigned int pTrack_id) const;

    /**
     *  @brief  Copy track states stored in tracks to pandora track parameters
     * 
     */
    void GetTrackStates(const edm4hep::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Copy track state from track state instance to pandora input track state
     * 
     */
    void CopyTrackState(const edm4hep::TrackState & pTrackState, pandora::InputTrackState &inputTrackState) const;

    /**
     *  @brief  Obtain track time when it reaches ECAL
     * 
     */
    float CalculateTrackTimeAtCalorimeter(const edm4hep::Track *const pTrack) const;

    /**
     *  @brief  Decide whether track reaches the ecal surface
     * 
     */
    void TrackReachesECAL(const edm4hep::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Determine whether a track can be used to form a pfo under the following conditions:
     *          1) if the track proves to be associated with a cluster, OR
     *          2) if the track proves to have no cluster associations
     * 
     */
    void DefineTrackPfoUsage(const edm4hep::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Whether track passes the quality cuts required in order to be used to form a pfo
     * 
     * 
     *  @return boolean
     */
    bool PassesQualityCuts(const edm4hep::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Get number of hits in TPC of a track
     * 
     */
    int GetNTpcHits(const edm4hep::Track *const pTrack) const;

    /**
     *  @brief  Get number of hits in FTD of a track
     * 
     */
    int GetNFtdHits(const edm4hep::Track *const pTrack) const;

    const Settings          m_settings;                     ///< The track creator settings
    const pandora::Pandora *m_pPandora;                     ///< Address of the pandora object to create tracks and track relationships

    float             m_bField;                       ///< The bfield

    float             m_tpcInnerR;                    ///< The tpc inner radius
    float             m_tpcOuterR;                    ///< The tpc outer radius
    unsigned int      m_tpcMaxRow;                    ///< The tpc maximum row number
    float             m_tpcZmax;                      ///< The tpc maximum z coordinate
    float                   m_cosTpc;                       ///< Cos(theta) value at end of tpc

    DoubleVector            m_ftdInnerRadii;                ///< List of ftd inner radii
    DoubleVector            m_ftdOuterRadii;                ///< List of ftd outer radii
    DoubleVector            m_ftdZPositions;                ///< List of ftd z positions
    unsigned int            m_nFtdLayers;                   ///< Number of ftd layers
    float                   m_tanLambdaFtd;                 ///< Tan lambda for first ftd layer

    int               m_eCalBarrelInnerSymmetry;      ///< ECal barrel inner symmetry order
    float             m_eCalBarrelInnerPhi0;          ///< ECal barrel inner phi 0
    float             m_eCalBarrelInnerR;             ///< ECal barrel inner radius
    float             m_eCalEndCapInnerZ;             ///< ECal endcap inner z

    float                   m_minEtdZPosition;              ///< Min etd z position
    float                   m_minSetRadius;                 ///< Min set radius

    TrackVector             m_trackVector;                  ///< The track vector
    TrackList               m_v0TrackList;                  ///< The list of v0 tracks
    TrackList               m_parentTrackList;              ///< The list of parent tracks
    TrackList               m_daughterTrackList;            ///< The list of daughter tracks
    TrackToPidMap           m_trackToPidMap;                ///< The map from track addresses to particle ids, where set by kinks/V0s
    gear::GearMgr* _GEAR;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackVector &TrackCreator::GetTrackVector() const
{
    return m_trackVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackCreator::Reset()
{
    m_trackVector.clear();
    m_v0TrackList.clear();
    m_parentTrackList.clear();
    m_daughterTrackList.clear();
    m_trackToPidMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsV0(unsigned int pTrack_id) const // should check here, if id is correct one to do this
{
    return (m_v0TrackList.end() != m_v0TrackList.find(pTrack_id));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsParent(unsigned int pTrack_id) const
{
    return (m_parentTrackList.end() != m_parentTrackList.find(pTrack_id));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackCreator::IsDaughter(unsigned int pTrack_id) const
{
    return (m_daughterTrackList.end() != m_daughterTrackList.find(pTrack_id));
}

#endif // #ifndef TRACK_CREATOR_H
