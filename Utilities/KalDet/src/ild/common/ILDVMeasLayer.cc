
#include "ILDVMeasLayer.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

// #include "streamlog/streamlog.h"

Bool_t   ILDVMeasLayer::kActive = kTRUE;
Bool_t   ILDVMeasLayer::kDummy = kFALSE;

//ClassImp(ILDVMeasLayer)

ILDVMeasLayer::ILDVMeasLayer(TMaterial &min,
                             TMaterial &mout,
                             Double_t   Bz,
                             Bool_t     is_active,
                             int        cellID ,
                             const Char_t    *name)  
: TVMeasLayer(min, mout, is_active, name),
_Bz(Bz),
_isMultiLayer(false)
{
  _cellIDs.push_back(cellID);

  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  encoder.setValue(cellID);
  encoder[lcio::ILDCellID0::module] = 0;
  encoder[lcio::ILDCellID0::sensor] = 0;

  _layerID = encoder.lowWord();
  
}


ILDVMeasLayer::ILDVMeasLayer(TMaterial &min,
                             TMaterial &mout,
                             Double_t  Bz,
                             const std::vector<int>& cellIDs,
                             Bool_t    is_active,
                             const Char_t    *name)
: TVMeasLayer(min, mout, is_active, name),
_Bz(Bz),
_cellIDs(cellIDs),
_isMultiLayer(true)
{
  
  if (cellIDs.size() == 0 ) {
    // streamlog_out(ERROR) << __FILE__ << " line " << __LINE__ << " size of cellIDs == 0" << std::endl;
  }

  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  encoder.setValue(cellIDs.at(0));
  encoder[lcio::ILDCellID0::module] = 0;
  encoder[lcio::ILDCellID0::sensor] = 0;

  _layerID = encoder.lowWord();
  
}



