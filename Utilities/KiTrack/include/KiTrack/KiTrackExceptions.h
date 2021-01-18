#ifndef KiTrackExceptions_h
#define KiTrackExceptions_h

#include <string>
#include <exception> 

//Exceptions for the KiTrack namespace

namespace KiTrack {

  /**Base exception class for KiTrack - all other exceptions extend this.
   * @author R. Glattauer, HEPHY
   * 
   */

  class KiTrackException : public std::exception {

    
  protected:
    std::string message{} ;
    
    KiTrackException(){  /*no_op*/ ; } 
    
  public: 
     virtual ~KiTrackException()  { /*no_op*/; } 
    
    KiTrackException( const std::string& text ){
      message = "KiTrack::Exception: " + text ;
    }

    virtual const char* what() const  noexcept { return  message.c_str() ; } 

  };

  
  
  /**Out of range exception, used when the user tries to access layers oder sensors, that are not implemented
   * @author R. Glattauer, HEPHY
   */
  class OutOfRange : public KiTrackException{
    
  protected:
    OutOfRange() {  /*no_op*/ ; } 
  public: 
    virtual ~OutOfRange()  { /*no_op*/; } 

    OutOfRange( std::string text ){
      message = "KiTrack::OutOfRange: " + text ;
    }
  }; 
  
  
  /**Invalid Parameter exception.
   * @author R. Glattauer, HEPHY
   */
  class InvalidParameter : public KiTrackException{
     
  protected:
     InvalidParameter() {  /*no_op*/ ; } 
  public: 
     virtual ~InvalidParameter()  { /*no_op*/; } 
     
     InvalidParameter( std::string text ){
        message = "KiTrack::InvalidParameter: " + text ;
     }
  }; 
  
  
  /**Wrong segment length exception.
   * @author R. Glattauer, HEPHY
   */
  class BadSegmentLength : public KiTrackException{
     
  protected:
     BadSegmentLength() {  /*no_op*/ ; } 
  public: 
     virtual ~BadSegmentLength()  { /*no_op*/; } 
     
     BadSegmentLength( std::string text ){
        message = "KiTrack::BadSegmentLength: " + text ;
     }
  }; 

  
  
  /**Unknown criterion exception.
   * @author R. Glattauer, HEPHY
   */
  class UnknownCriterion : public KiTrackException{
     
  protected:
     UnknownCriterion() {  /*no_op*/ ; } 
  public: 
     virtual ~UnknownCriterion()  { /*no_op*/; } 
     
     UnknownCriterion( std::string text ){
        message = "KiTrack::UnknownCriterion: " + text ;
     }
  }; 

} 

#endif


