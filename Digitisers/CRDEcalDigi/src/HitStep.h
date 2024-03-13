#ifndef HIT_STEP_H
#define HIT_STEP_H

class HitStep{

public: 
  HitStep (double _Q, double _T): Q(_Q), T(_T) {};
  //HitStep (double _Q, double _T) { Q=_Q; T=_T; };
  HitStep() {};

  void setQ(double _Q) { Q =_Q; }
  void setT(double _T) { T =_T; }

  double getQ() const { return Q; }
  double getT() const { return T; }
  inline bool operator < (const HitStep &x) const { 
    return T <x.T ;
  }

private: 
  double Q;
  double T;

};

#endif
