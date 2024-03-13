#ifndef _CRD_CALOHIT_
#define _CRD_CALOHIT_

#include <DD4hep/Objects.h>
#include "TVector3.h"

class CaloHit{

public:
  CaloHit(unsigned long long _cellID, int _system, int _module, int _layer, int _stave, int _tower, int _slice, TVector3 _pos, double _E, int _step_numbers)
  : cellID(_cellID), system(_system), module(_module), layer(_layer), stave(_stave), tower(_tower), slice(_slice), position(_pos), E(_E), step_numbers(_step_numbers) {}; 
  CaloHit() {};

  inline bool operator == (const CaloHit &x) const{
    return ( (cellID == x.cellID) && getEnergy()==x.getEnergy() );
  }
  unsigned long long getcellID() const { return cellID; }
  int getSystem() const { return system; }
  int getModule() const { return module; }
  int getStave()  const { return stave;  }
  int getLayer()  const { return layer;  }
  int getTower()  const { return tower;  }
  int getSlice()  const { return slice;  }

  TVector3 getPosition() const { return position; }
  double getEnergy() const { return E; }
  int getStep_numbers() const { return step_numbers; }

  void setcellID(unsigned long long _cellid) { cellID = _cellid; }
  void setcellID(int _system, int _module, int _layer, int _stave, int _tower, int _slice) { system=_system; module=_module; stave=_stave; layer=_layer; tower=_tower; slice=_slice;}
  void setPosition( TVector3 posv3) { position.SetXYZ( posv3.x(), posv3.y(), posv3.z() ); }
  void setE(double _E) { E=_E; }
  void setStep_numbers(int _step_numbers) { step_numbers=_step_numbers; }


private:
	unsigned long long cellID;
	int system;
	int module;
	int stave;
	int layer;
	int tower;
	int slice;
	TVector3 position;
	double E;
  int step_numbers;

};
  
#endif
