/**
 * 
 *  @brief  Header file for the mc particle creator class.
 * 
 *  $Log: $
 */

#ifndef MC_PARTICLE_CREATOR_H
#define MC_PARTICLE_CREATOR_H 1

#include "edm4hep/MCParticle.h"
#include "Api/PandoraApi.h"

#include "CaloHitCreator.h"
#include "TrackCreator.h"
/**
 *  @brief  MCParticleCreator class
 */

class CollectionMaps;

class MCParticleCreator
{
public:
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

        StringVector    m_mcParticleCollections;                ///< The mc particle collections
        StringVector    m_CaloHitRelationCollections;         ///< The SimCaloHit to CaloHit particle relations
        StringVector    m_TrackRelationCollections;           ///< The SimTrackerHit to TrackerHit particle relations
        float           m_bField;                             ///< m_bField
        bool            m_debug;
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~MCParticleCreator();

    /**
     *  @brief  Create MCParticles
     * 
     */    
    pandora::StatusCode CreateMCParticles(const CollectionMaps& collectionMaps ) const;

    /**
     *  @brief  Create Track to mc particle relationships
     *
     */
     pandora::StatusCode CreateTrackToMCParticleRelationships(const CollectionMaps& collectionMaps, const TrackVector &trackVector) const;

     void Reset();
    /**
     *  @brief  Create calo hit to mc particle relationships
     *
     */
      pandora::StatusCode CreateCaloHitToMCParticleRelationships(const CollectionMaps& collectionMaps, const CalorimeterHitVector &calorimeterHitVector) const;

private:
    const Settings          m_settings;                         ///< The mc particle creator settings
    const pandora::Pandora *m_pPandora;                         ///< Address of the pandora object to create the mc particles
    const float             m_bField;                           ///< The bfield
    std::map<unsigned int, const edm4hep::MCParticle*>*  m_id_pMC_map;
};

inline void MCParticleCreator::Reset()
{
    m_id_pMC_map->clear();
}

#endif // #ifndef MC_PARTICLE_CREATOR_H
