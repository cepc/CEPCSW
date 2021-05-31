#ifndef GenericBFieldMapBrBz_h
#define GenericBFieldMapBrBz_h

/*
 * GenericBFieldMapBrBz is an extension of Cartesian Field in DD4hep.
 * - It enables the DD4hep to access the Magnetic Service from Gaudi. 
 * - It also enables the calculation of Br/Bz at position X. 
 * - It will get the map from an abstract class IFieldMapProvider.
 *
 * -- Tao Lin <lintao AT ihep.ac.cn>
 */

#include <DD4hep/FieldTypes.h>

#include "IFieldMapProvider.h"


class GenericBFieldMapBrBz: public dd4hep::CartesianField::Object {
public:

    GenericBFieldMapBrBz();

    virtual void fieldComponents(const double* pos, double* field);

public:
    // following are interfaces to configure this field map
    void init_provider(const std::string& provider, const std::string& url);

private:

    IFieldMapProvider* m_provider;
};

#endif

