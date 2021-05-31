
#include "GenericBFieldMapBrBz.h"

#include "FieldMapFileProvider.h"

#include <cmath>

#include "DD4hep/DD4hepUnits.h"

GenericBFieldMapBrBz::GenericBFieldMapBrBz()
    : m_provider(nullptr) {
    type = dd4hep::CartesianField::MAGNETIC;

}

void GenericBFieldMapBrBz::fieldComponents(const double* pos, double* field) {
    if (!m_provider) {
        std::string error_msg = "[ERROR] GenericBFieldMapBrBz: No provider! ";
        throw std::runtime_error(error_msg); 
    }

    // convert pos to r/z
    double x = pos[0] / dd4hep::m; // convert to meter
    double y = pos[1] / dd4hep::m;
    double z = pos[2] / dd4hep::m;
    double r = sqrt(x*x+y*y);
    double phi = atan2(y, x);

    // std::cout << " r: " << r
    //           << " x: " << x
    //           << " y: " << y
    //           << " z: " << z
    //           << " field Bx: " << field[0]
    //           << " field By: " << field[1]
    //           << " field Bz: " << field[2]
    //           << std::endl;

    // get the point r0/z0, r1/z1

    // -- normalized r/z at the bin
    double rn = 0.0;
    double zn = 0.0;

    int ir0 = 0;
    int iz0 = 0;

    ir0 = m_provider->rBinIdx(r, rn);
    iz0 = m_provider->zBinIdx(z, zn);

    // not in the valid return
    if (ir0 < 0 || iz0 < 0) {
        // std::cout << "SKIP due to "
        //           << " r: " << r
        //           << " x: " << x
        //           << " y: " << y
        //           << " z: " << z
        //           << " ir0: " << ir0
        //           << " iz0: " << iz0
        //           << " rn: " << rn
        //           << " zn: " << zn
        //           << std::endl;
        return;
    }

    // std::cout << " r: " << r
    //           << " x: " << x
    //           << " y: " << y
    //           << " z: " << z
    //           << " ir0: " << ir0
    //           << " iz0: " << iz0
    //           << " rn: " << rn
    //           << " zn: " << zn
    //           << std::endl;

    int ir1 = ir0 + 1;
    int iz1 = iz0 + 1;

    // get the Br/Bz at four points
    double Br_r0z0 = 0.0;    double Bz_r0z0 = 0.0;
    double Br_r1z0 = 0.0;    double Bz_r1z0 = 0.0;
    double Br_r0z1 = 0.0;    double Bz_r0z1 = 0.0;
    double Br_r1z1 = 0.0;    double Bz_r1z1 = 0.0;

    // calculate the field at r/z
    m_provider->access(ir0, iz0, Br_r0z0, Bz_r0z0);
    m_provider->access(ir1, iz0, Br_r1z0, Bz_r1z0);
    m_provider->access(ir0, iz1, Br_r0z1, Bz_r0z1);
    m_provider->access(ir1, iz1, Br_r1z1, Bz_r1z1);

    double Br = (1.0 - rn) * (1.0 - zn) * Br_r0z0 
              +        rn  * (1.0 - zn) * Br_r1z0
              + (1.0 - rn) *        zn  * Br_r0z1
              +        rn  *        zn  * Br_r1z1;

    double Bz = (1.0 - rn) * (1.0 - zn) * Bz_r0z0 
              +        rn  * (1.0 - zn) * Bz_r1z0
              + (1.0 - rn) *        zn  * Bz_r0z1
              +        rn  *        zn  * Bz_r1z1;

    // update the global field
    field[0] += Br*cos(phi);
    field[1] += Br*sin(phi);
    field[2] += Bz;

    return;
}

void GenericBFieldMapBrBz::init_provider(const std::string& provider, const std::string& url) {
    if (provider == "file") {
        std::cout << "Initialize provider with file. " << std::endl;
        m_provider = new FieldMapFileProvider(url);
    } else if (provider == "db") {
        std::cout << "Initialize provider with file. " << std::endl;
    } else {
        std::string error_msg = "[ERROR] GenericBFieldMapBrBz: Unknown provider: " + provider;
        throw std::runtime_error(error_msg); 
    }
}
