#include "LCAscHepRdr.h"
#include "EVENT/MCParticle.h"
#include "IMPL/MCParticleImpl.h"
#include "lcio.h"
#include "EVENT/LCIO.h"
#include <sstream>
#include "Exceptions.h"

using namespace EVENT ;
using namespace IMPL ;

//
// the time being we leave it as #define, to avoid to depend
// on Mokka business if it could be part of lcio package
//
#define HEPEvt 1
#define hepevt 2

namespace UTIL{
  
  LCAscHepRdr::LCAscHepRdr(const char* evfile, int fileFormat)
    : theFileFormat(fileFormat)
  {
    //
    //   Adapt the Geant4 G4HEPEvtInterface code 
    //
    switch(fileFormat) 
      {
      case HEPEvt :
      case hepevt :
	inputFile.open(evfile);
	if (!inputFile) 
	  {
	    std::stringstream description ; 
	    description << "LCAscHepRdr, no ascii Hep file found: " << evfile << std::ends ;
	    throw IO::IOException( description.str() );
	  }
	break;
      default :
	std::stringstream description ; 
	description << "LCAscHepRdr, bad suffix on file name: " << evfile << std::ends ;
	throw IO::IOException( description.str() );
      }
  }

  LCAscHepRdr::~LCAscHepRdr(){}
  
//
// Read an event and return a LCCollectionVec of MCParticles
//
  IMPL::LCCollectionVec * LCAscHepRdr::readEvent()
  {

    IMPL::LCCollectionVec * mcVec = nullptr;
    double c_light = 299.792;// mm/ns
    //
    //  Read the event, check for errors
    //
    int NHEP;  // number of entries
    int NOUT ;   // number of outgoing particles
    int BRE ;   // beam remnants
    double WEIGHT ;   // weight
    inputFile >> NHEP >> NOUT >> BRE >> WEIGHT; 
    if( inputFile.eof() ) 
      {
	//
	// End of File :: ??? Exception ???
	//   -> FG:   EOF is not an exception as it happens for every file at the end !
	//
	return mcVec;
      }
    
    //
    //  Create a Collection Vector
    //
    mcVec = new IMPL::LCCollectionVec(LCIO::MCPARTICLE);
    MCParticleImpl* p;
    MCParticleImpl* d;
    
    //
    //  Loop over particles
    //
    int ISTHEP;   // status code
    int IDHEP;    // PDG code
    int JMOHEP1;  // first mother
    int JMOHEP2;  // last mother
    int JDAHEP1;  // first daughter
    int JDAHEP2;  // last daughter
    double PHEP1; // px in GeV/c
    double PHEP2; // py in GeV/c
    double PHEP3; // pz in GeV/c
    double PHEP4; // energy in GeV
    double PHEP5; // mass in GeV/c**2
    double VHEP1; // x vertex position in mm
    double VHEP2; // y vertex position in mm
    double VHEP3; // z vertex position in mm
    double VHEP4; // production time in mm/c

    std::vector<int> *daughter1 = new std::vector<int> ();
    std::vector<int> *daughter2 =  new std::vector<int> ();

    for( int IHEP=0; IHEP<NHEP; IHEP++ )
      {
	//if ( theFileFormat == HEPEvt)
	if ( false)
	  inputFile >> ISTHEP >> IDHEP >> JDAHEP1 >> JDAHEP2
		  >> PHEP1 >> PHEP2 >> PHEP3 >> PHEP5;
	else
	  inputFile >> ISTHEP >> IDHEP 
		    >> JMOHEP1 >> JMOHEP2
		    >> JDAHEP1 >> JDAHEP2
		    >> PHEP1 >> PHEP2 >> PHEP3 
		    >> PHEP4 >> PHEP5
		    >> VHEP1 >> VHEP2 >> VHEP3
		    >> VHEP4;

	if(inputFile.eof())
		return nullptr;	
	//
	//  Create a MCParticle and fill it from stdhep info
	//
	MCParticleImpl* mcp = new MCParticleImpl();
	//
	//  PDGID
	//
	mcp->setPDG(IDHEP);
	//
	//  Momentum vector
	//
	float p0[3] = {PHEP1,PHEP2,PHEP3};
	mcp->setMomentum(p0);
	//
	//  Mass
	//
	mcp->setMass(PHEP5);
	//
	//  Vertex 
	// (missing information in HEPEvt files)
	double v0[3] = {0.,0.,0.};
	mcp->setVertex(v0);
	//
	//  Generator status
	//
	mcp->setGeneratorStatus(ISTHEP);
	//
	//  Simulator status 0 until simulator acts on it
	//
	mcp->setSimulatorStatus(0);
	//
	//  Creation time (note the units)
	// (No information in HEPEvt files)
	mcp->setTime(0./c_light);
	//
	//  Add the particle to the collection vector
	//
	mcVec->push_back(mcp);
	//
	// Keep daughters information for later
	//
	daughter1->push_back(JDAHEP1);
	daughter2->push_back(JDAHEP2);

	   // fg: comment out the  mother relationships altogether ....

// //
// // Add the parent information. The implicit assumption here is that
// // no particle is read in before its parents.
// //
// 	   int fp = _reader->mother1(IHEP) - 1;
// 	   int lp = _reader->mother2(IHEP) - 1;
// //
// //  If both first parent and second parent > 0, and second parent >
// //     first parent, assume a range
// //
// 	   if( (fp > -1) && (lp > -1) )
// 	   {
// 	     if(lp >= fp)
// 		 {
// 		   for(int ip=fp;ip<lp+1;ip++)
// 		   {
// 			 p = dynamic_cast<MCParticleImpl*>
// 				 (mcVec->getElementAt(ip));
// 			 mcp->addParent(p);
// 		   }
// 		 }
// //
// //  If first parent < second parent, assume 2 discreet parents
// //
// 		 else
// 		 {
// 		   p = dynamic_cast<MCParticleImpl*>
// 			 (mcVec->getElementAt(fp));
// 		   mcp->addParent(p);
// 		   p = dynamic_cast<MCParticleImpl*>
// 			 (mcVec->getElementAt(lp));
// 		   mcp->addParent(p);
// 		 }
// 	   }
// //
// //  Only 1 parent > 0, set it
// //
// 	   else if(fp > -1)
// 	   {
// 		 p = dynamic_cast<MCParticleImpl*>
// 		   (mcVec->getElementAt(fp));
// 		 mcp->addParent(p);
// 	   }
// 	   else if(lp > -1)
// 	   {
// 		 p = dynamic_cast<MCParticleImpl*>
// 		   (mcVec->getElementAt(lp));
// 		 mcp->addParent(p);
// 	   }

	 }// End loop over particles
//
//  Now make a second loop over the particles, checking the daughter
//  information. This is not always consistent with parent 
//  information, and this utility assumes all parents listed are
//  parents and all daughters listed are daughters
//
    for( int IHEP=0; IHEP<NHEP; IHEP++ )
      {
	//
	//  Get the MCParticle
	//
	MCParticleImpl* mcp = 
	  dynamic_cast<MCParticleImpl*>
	  (mcVec->getElementAt(IHEP));
	//
	//  Get the daughter information, discarding extra information
	//  sometimes stored in daughter variables.
	//

	int fd = daughter1->operator[](IHEP) - 1;
	int ld = daughter2->operator[](IHEP) - 1;

	//
	//  As with the parents, look for range, 2 discreet or 1 discreet 
	//  daughter.
	//
	if( (fd > -1) && (ld > -1) )
	  {
	    if(ld >= fd)
	      {
		for(int id=fd;id<ld+1;id++)
		  {
		    //
		    //  Get the daughter, and see if it already lists this particle as
		    //    a parent.
		    //
		    d = dynamic_cast<MCParticleImpl*>
		      (mcVec->getElementAt(id));
		    int np = d->getParents().size();
		    bool gotit = false;
		    for(int ip=0;ip < np;ip++)
		      {
			p = dynamic_cast<MCParticleImpl*>
			  (d->getParents()[ip]);
			if(p == mcp)gotit = true;
		      }
		    //
		    //  If not already listed, add this particle as a parent
		    //
		    if(!gotit)d->addParent(mcp);
		  }
	      }
	    //
	    //  Same logic, discreet cases
	    //
	    else
	      {
		d = dynamic_cast<MCParticleImpl*>
		  (mcVec->getElementAt(fd));
		int np = d->getParents().size();
		bool gotit = false;
		for(int ip=0;ip < np;ip++)
		  {
		    p = dynamic_cast<MCParticleImpl*>
		      (d->getParents()[ip]);
		    if(p == mcp)gotit = true;
		  }
		if(!gotit)d->addParent(mcp);
		d = dynamic_cast<MCParticleImpl*>
		  (mcVec->getElementAt(ld));
		np = d->getParents().size();
		gotit = false;
		for(int ip=0;ip < np;ip++)
		  {
		    p = dynamic_cast<MCParticleImpl*>
		      (d->getParents()[ip]);
		    if(p == mcp)gotit = true;
		  }
		if(!gotit)d->addParent(mcp);
	      }
	  }
	else if(fd > -1)
	  {
	    d = dynamic_cast<MCParticleImpl*>
	      (mcVec->getElementAt(fd));
	    int np = d->getParents().size();
	    bool gotit = false;
	    for(int ip=0;ip < np;ip++)
	      {
		p = dynamic_cast<MCParticleImpl*>
		  (d->getParents()[ip]);
		if(p == mcp)gotit = true;
	      }
	    if(!gotit)d->addParent(mcp);
	  }
	else if(ld > -1)
	  {
	    d = dynamic_cast<MCParticleImpl*>
	      (mcVec->getElementAt(ld));
	    int np = d->getParents().size();
	    bool gotit = false;
	    for(int ip=0;ip < np;ip++)
	      {
		p = dynamic_cast<MCParticleImpl*>
		  (d->getParents()[ip]);
		if(p == mcp)gotit = true;
	      }
	    if(!gotit)d->addParent(mcp);
	  }
      }// End second loop over particles
    //
    //  Return the collection
    //
    return mcVec;    
  }
  
  
} // namespace UTIL
