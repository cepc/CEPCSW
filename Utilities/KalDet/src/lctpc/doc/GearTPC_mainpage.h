/**
@mainpage

GearTPCKalDetector is a KalDet implementation for the LCTPC collaboration 
which is flexible enough to handle multiple detector geometries 
without having to adapt the code. It instantiates the full KalDet geometry only from information from the Gear file.
The design is flexible enough that any geometry which can be described with gear
can also be treated by the Kalman detector. Currently not all parts are
implemented, for instance straight pad rows / Kalman layers. In this
case an error message is printed and a gear::NotImplementedException is thrown.

For backward compatibility the EXTPC interface has been reimplemented using 
the GearTPC classes. Once all dependencies on EXTPC have been removed from
MarlinTPC these compatibility classes will be removed from KalDet.
 */
