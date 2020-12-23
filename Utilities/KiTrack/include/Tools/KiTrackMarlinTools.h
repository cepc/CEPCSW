#ifndef KiTrackMarlinTools_h
#define KiTrackMarlinTools_h

#include <string>
#include <set>
#include <vector>
#include <map>
//#include <sstream>

#include "edm4hep/TrackerHit.h"
#include "edm4hep/Track.h"

#include "ILDImpl/FTDHitSimple.h"
#include "ILDImpl/VXDHitSimple.h"
#include "KiTrack/ITrack.h"

using namespace KiTrack;

namespace KiTrackMarlin{

/** @return information about the contents of the passed CellID0 */ 
std::string getCellID0Info( int cellID0 );

/** @return the layer given by the cellID0 */
int getCellID0Layer( int cellID0 );

/** Generates a root file with a specified name and with a single tree .
 * If a rootfile with the same name already exists, a number will be appended
 * to it, so the old one doesn't get lost.
 * So if I create myRoot.root and it already exists, the old one will be renamed to
 * "myRoot1.root".
 * If "myRoot1.root" already exists, "myRoot2.root" will be created.
 * 
 * So after we did this for some time, we end up with a lot of root files with numbers.
 * The one without a number appended is the newest one.
 * The one with the number 1 appended is the oldest, 2 is newer than 1, 3 is newer than 2 and so on.
 * 
 * @param branchNames A set of all the branchnames that should be created in the tree
 * 
 * @param createNew Whether to create a new root file or create the tree and branches in an existing one
 */
void setUpRootFile( std::string fileNamePath, std::string treeName , std::set<std::string> branchNames = std::set<std::string>() , bool createNew=true );


/** Saves values to a tree in a rootfile.
 * 
 * @param fileNamePath the name (path) of the rootfile
 * 
 * @param treeName the name of the tree
 * 
 * @param map_name_data a map of key = name of what we save , mapped value = value we want to save
 * 
 */
void saveToRoot( std::string fileNamePath, std::string treeName , std::map < std::string , float > map_name_data );

void saveToRoot( std::string rootFileName, std::string treeName , std::vector < std::map < std::string , float > > rootDataVec );




/** method that compares two TrackerHits.
 * 
 * @return true if |a.z| < |b.z| , i.e. true if a is nearer to z=0 than b is
 */
//bool compare_TrackerHit_z( edm4hep::ConstTrackerHit* a, edm4hep::ConstTrackerHit* b );
bool compare_TrackerHit_z( edm4hep::ConstTrackerHit& a, edm4hep::ConstTrackerHit& b );

/** method that compares two TrackerHits.
 * 
 * @return true if |a.R| < |b.R| , i.e. true if a has smaller radius than b 
 *
 * to be used at the VXD-SIT system
 */
bool compare_TrackerHit_R( edm4hep::ConstTrackerHit& a, edm4hep::ConstTrackerHit& b );


FTDHitSimple* createVirtualIPHit( int side , const SectorSystemFTD* sectorSystemFTD );

VXDHitSimple* createVirtualIPHit( const SectorSystemVXD* sectorSystemVXD );


std::string getPositionInfo( edm4hep::ConstTrackerHit hit );

std::string getPositionInfo( IHit* hit );   

std::string getTrackHitInfo( ITrack* track );

std::string getTrackHitInfo( edm4hep::Track* track );

} // end of namespace KiTrackMarlin




#endif


