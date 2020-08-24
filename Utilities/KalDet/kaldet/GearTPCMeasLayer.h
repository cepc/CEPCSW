#ifndef GEARTPC_MEASLAYER_H
#define GEARTPC_MEASLAYER_H

#include "kaltest/TVMeasLayer.h"
#include <set>

namespace kaldet
{

  class GearTPCHit;

  /**
   * The GearTPCMeasLayer class introduces the z-dependent resolutions sigmaX and sigmaZ
   * as well as Gear modules and rows which correspond to this layer.
   *
   * If the layer is defined as a perfect layer this means all modules are perfectly alligned 
   * and more than one module/row can be assigned to this layer. You can add them using AddModuleRow.
   * The perfect layer should contain all the moduleRows on it, so it is guaranteed that the
   * user can access all neighbouring modules this way. 
   *
   * If the layer is not defined as perfect (default) there can only be one module on this layer.
   * Calling AddModuleRow will throw an exception. This is the default behaviour because Gear does
   * not guarantee that the modules are alligned. Displaced modules do not make up a perfect
   * cylinder / plane and have to be treated as separate segments. Finding a neighbouring module/row
   * is not trivial and has to be left to the user or a future Gear version.
   */
		 
  class GearTPCMeasLayer 
    : public TVMeasLayer
  {

  public:
    /** The constructor.
     *  The materials and the type (active or passive) are passed on to the
     *  TVMeasLayer. sigmaX0 [mm] is the constant part of sigmaX, sigmaX1 [mm/sqrt(mm)] 
     *  the z-dependent part, accordingly for sigmaZ.
     *
     *  Module and row have to be specified. They will be added as the first
     *  module/row pair of this measurement layer. 
     *  For a perfect layer modules can be added with AddModuleRow.
     *
     *  Note: This class cannot be instantiated because the parent's geometry dependent
     *  purely virtual
     *  functions like XvToMv are not implemented. This will happen in the cylindrical or planar
     *  implementations.
     *
     *  For inactive layers you will usually leave the sigmas at 0, they have no useful meaning in 
     *  this case.
     */
    GearTPCMeasLayer(TMaterial &min,
		     TMaterial &mout,
		     Int_t      module,
		     Int_t      row,
		     Bool_t     isPerfect,
		     Bool_t     isActive,
		     Double_t   sigmaX0 = 0., //< the constant part of sigmaX
		     Double_t   sigmaX1 = 0., //< the z-dependent part of sigmaX
		     Double_t   sigmaZ0 = 0. , //< the constant part of sigmaZ
		     Double_t   sigmaZ1 = 0.); //< the z-dependent part of sigmaZ
    
    /// The destructor
    virtual  ~GearTPCMeasLayer();

   /**
   * A perfect measurement layer contains all the modules with rows (row segments)
   * that make up the layer.
   */
    virtual std::set< std::pair <int, int> > const & GetModuleRows() const;
    
    /**
     * Add another row on another module which lies on the same cylinder.
     */
    virtual void AddModuleRow(int module, int row);
    
    /**
     * Get the measurement vector (mv) for this layer from a space point (xv)
     */
    virtual TKalMatrix XvToMv    (const TVector3   &xv)   const = 0;

    /**
     * Get the z-depenent resolution in the readout plane 
     * (usually x or r\f$\phi\f$).
     */
    virtual Double_t GetSigmaX(Double_t z) const;

     /**
     * Get the z-depenent resolution in z (drift direction).
     */
   virtual Double_t GetSigmaZ(Double_t z) const;
    
  
    /**
     * Get the flag whether the layer is declared as perfect.
     */
    virtual Bool_t IsPerfect() const;

    /**
     * A virtual function to create the appropriate hit. Depending on the implementation
     * (cylindrical or straight measurement layer) you get the appropriate implementation 
     * of GearTPCHit.
     * It creates a new hit on the heap and hands over the ownership.
     */
    virtual GearTPCHit * createHit(Double_t * meas,
				   Double_t * dmeas,
				   void * hitPointer, 
				   Double_t bField,
				   Double_t vDrift,
				   Int_t           m = kMdim) const = 0;

  protected:
    Double_t fSigmaX0;  // xy resolution
    Double_t fSigmaX1;  // xy resolution
    Double_t fSigmaZ0;  // z  resolution
    Double_t fSigmaZ1;  // z  resolution

    /// A set to hold all the module/row combinations associated to this layer
    std::set< std::pair<int, int> > fModuleRows;

    Bool_t fIsPerfect;
  };

}// namespace kaldet
#endif // GEARTPC_MEASLAYER_H
