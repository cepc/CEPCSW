#ifndef UTIL_LCAscHepRdr_H
#define UTIL_LCAscHepRdr_H 1

#include "IMPL/LCCollectionVec.h"
#include <fstream>

namespace UTIL{
  
  /**Basic utility for reading a ASCII HEPEvt file and filling
   * a LCCollectionVec with MCParticles containing the HEPEvt
   * file information.
   * 
   * @author Mora de Freitas
   * @version $Id: 
   */
  class LCAscHepRdr{
    
  public:

    /** Open the HEPEvt input file in the constructer
     */
    LCAscHepRdr(const char* evfile, int fileFormat) ;
    
    /** noop
     */
    ~LCAscHepRdr() ;
    
    /** Read an event and return a LCCollectionVec of MCParticles.
     */
    IMPL::LCCollectionVec * readEvent() ;

  private:
    std::ifstream inputFile;
    int theFileFormat;
    
  }; // class

} // namespace UTIL

#endif /* ifndef UTIL_LCAscHepRdr_H */
