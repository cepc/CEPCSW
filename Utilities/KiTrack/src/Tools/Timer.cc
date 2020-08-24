#include "Tools/Timer.h"

#include <string>
#include <fstream>
#include <iterator>
#include <cstdlib>
#include <iostream>

#define FINE_GRAINED_LINUX_TIMING

#include <ctime>
#include <sys/times.h>
#include <unistd.h>

using namespace std;
using namespace KiTrackMarlin;

namespace {
  #ifdef FINE_GRAINED_LINUX_TIMING
  static unsigned cyc_hi = 0;
  static unsigned cyc_lo = 0;

  void access_counter(unsigned *hi, unsigned *lo)
  {
    asm("rdtsc; movl %%edx,%0; movl %%eax,%1"
        : "=r" (*hi), "=r" (*lo)
        :
        : "%edx", "%eax");
  }

  double get_counter()
  {
    unsigned ncyc_hi, ncyc_lo;
    unsigned hi, lo, borrow;
    access_counter(&ncyc_hi, &ncyc_lo);
    lo = ncyc_lo - cyc_lo;
    borrow = lo > ncyc_lo;
    hi = ncyc_hi - cyc_hi - borrow;
    return (double) hi * (1 << 30) * 4 + lo;
  }
  #else
  clock_t before;
  #endif
}

double Timer::cpuMHz()
{
  static double ret=-1.;
  if ( ret > -.5 ) return ret;
  string input="Unknow CPU";
  {
    ifstream cpuinfo("/proc/cpuinfo");
    if ( cpuinfo.is_open() )
    {
      cpuinfo.unsetf( ios::skipws );
      istream_iterator<char> sbegin(cpuinfo),send;
      copy(sbegin,send,inserter(input,input.end()));
      cpuinfo.close();
    };
  }

  size_t i = input.find("model name");
  if (i!=string::npos )
  {
    i = input.find("@ ",i);
    if ( i != string::npos )
    {
      size_t j = input.find("GHz",i);
      //cout << "[debug] i=" << i << " j=" << j << " string=" << input.substr(i+1,j-i-1) << "X" << endl;
      return 1000.*atof(input.substr(i+1,j-i-1).c_str() );
    };

  }

  i = input.find("cpu MHz");
  if (i==string::npos)
  {
    cout << "[Timer] /proc/cpuinfo does not contain cpu speed..." << endl;
    ret=0.;
  };
  i = input.find(":",i);
  ret=atof(input.substr(i+1,input.find("/n",i)-i).c_str());
  return ret;
}

void Timer::start_counter()
{
  #ifdef FINE_GRAINED_LINUX_TIMING
  access_counter(&cyc_hi, &cyc_lo);
  #else
  tms buf;
  times( &buf );     
  before = buf.tms_utime;
  #endif
}

double Timer::lap()
{
  #ifdef FINE_GRAINED_LINUX_TIMING
  return get_counter() / cpuMHz() / 1.e6;
  #else
  tms buf;
  times ( &buf );
  clock_t after = buf.tms_utime;
  static float cps=sysconf(_SC_CLK_TCK);
  double t = ( after - before )  / cps;    
  return t;
  #endif
}

double Timer::ticks()
{
  #ifdef FINE_GRAINED_LINUX_TIMING
  return get_counter();
  #else
  tms buf;
  times ( &buf );
  clock_t after = buf.tms_utime;
  // static float cps=sysconf(_SC_CLK_TCK);
  double t = ( after - before ); 
  return t;
  #endif
}

#ifdef FINE_GRAINED_LINUX_TIMING
#undef FINE_GRAINED_LINUX_TIMING
#endif
