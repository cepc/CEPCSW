#include "MarlinEvtSeeder.h"
#include "jenkinsHash.h"
#include "GaudiKernel/Algorithm.h"

#include <stdlib.h>
#include <limits>
#include <algorithm>

DECLARE_COMPONENT(MarlinEvtSeeder)

MarlinEvtSeeder::MarlinEvtSeeder(const std::string& name, ISvcLocator* svc) 
  : base_class(name, svc),
    _evtNo(-1),
    _runNo(-1)
{
} 

void MarlinEvtSeeder::registerAlg ( Algorithm* alg )
{
  srand( _globalSeed );
  debug() << "initialised with global seed " << _globalSeed << endmsg; 

  _vAlgSeed.push_back( std::make_pair( alg, rand() ) );
  debug() << alg->name() << "registered for random seed service. Allocated "
          <<  _vAlgSeed.back().second << " as initial seed."
          << endmsg; 
}

void MarlinEvtSeeder::refreshSeeds()
{
  // get hashed seed using jenkins_hash
  unsigned int seed = 0 ; // initial state
  unsigned int eventNumber = _evtNo;  //evt->getEventNumber() ;
  unsigned int runNumber = _runNo;  //evt->getRunNumber() ;

  unsigned char * c = (unsigned char *) &eventNumber ;
  seed = jenkins_hash( c, sizeof eventNumber, seed) ;

  c = (unsigned char *) &runNumber ;
  seed = jenkins_hash( c, sizeof runNumber, seed) ;

  int _global_seed = _globalSeed;
  c = (unsigned char *) &_global_seed ;
  seed = jenkins_hash( c, sizeof _global_seed, seed) ;

  // set the seed for rand() for this event
  if ( seed == 1 ) seed = 123456789 ; // can't used a seed value of 1 as srand(1) sets rand() back the state of the last call to srand( seed ).

  debug() << "MarlinEvtSeeder: Refresh Seeds using " << seed << " as seed for srand( seed )" << endmsg; 
  srand( seed );

  // fill vector with seeds for each registered processor using rand() 
  for( auto& iAlgSeed : _vAlgSeed ) 
  {
    iAlgSeed.second = rand();
  }
}

unsigned int MarlinEvtSeeder::getSeed(Algorithm* alg, int eventNumber, int runNumber)
{
  if ( eventNumber != _evtNo || runNumber != _runNo ) {
    _evtNo = eventNumber;
    _runNo = runNumber;
    refreshSeeds();
  }

  typedef std::pair<Algorithm*, unsigned int> Pair;

  auto it = find_if( _vAlgSeed.begin(), _vAlgSeed.end(), [&](Pair const& pair){ return pair.first == alg;  }  );

  return ( it != _vAlgSeed.end() ? it->second : throw  ) ;
}

StatusCode MarlinEvtSeeder::initialize()
{
  return StatusCode::SUCCESS;
}

StatusCode MarlinEvtSeeder::finalize()
{
  return StatusCode::SUCCESS;
}
