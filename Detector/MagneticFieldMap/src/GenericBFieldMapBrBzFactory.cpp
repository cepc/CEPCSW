/*
 * In this file, the xml is parsed and the GenericBFieldMapBrBz object is created and configured.
 *
 * -- Tao Lin <lintao AT ihep.ac.cn>
 */
#include "GenericBFieldMapBrBz.h"


#include <DD4hep/Version.h>
#if DD4HEP_VERSION_GE(0,24)
#include <DD4hep/detail/Handle.inl>
#else
#include <DD4hep/Handle.inl>
#endif

#include <DD4hep/FieldTypes.h>


#include <DD4hep/DetFactoryHelper.h>
#include <XML/Utilities.h>


static dd4hep::Ref_t create_GenericBFieldMapBrBz(dd4hep::Detector& , 
                                                 dd4hep::xml::Handle_t handle ) {

    // 1. retrieve the parameters from xml
    dd4hep::xml::Component xmlParameter(handle);

    // 2. create the CartesianField
    dd4hep::CartesianField obj;
    GenericBFieldMapBrBz* ptr = new GenericBFieldMapBrBz();

    obj.assign(ptr, xmlParameter.nameStr(), xmlParameter.typeStr());

    return obj;

}

DECLARE_XMLELEMENT(GenericBFieldMapBrBz,create_GenericBFieldMapBrBz)
