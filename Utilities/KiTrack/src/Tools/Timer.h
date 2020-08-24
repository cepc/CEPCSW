#ifndef dataharvester_Timer_H_
#define dataharvester_Timer_H_

namespace KiTrackMarlin {
class Timer {
public:
  /**
   *  A fast, precise timer
   *  Works on linux only
   * 
   * Author: Wolfgang Waltenberger
   */
  static void start_counter();

  /// Return what the harvester thinks is the CPU frequency.
  static double cpuMHz();
  static double lap(); //< lapsed time, in seconds.
  static double ticks(); //< lapsed time, in clock ticks.
};
}

#endif // dataharvester_Timer_H_
