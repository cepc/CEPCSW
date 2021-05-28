/*
 * In this file, the xml is parsed and the GenericBFieldMapBrBz object is created and configured.
 * 
 * The properties for the GenericBFieldMapBrBz
 * - provider (attribute)
 *   - [file, db]
 * - source (tag)
 *   - the attributes include: 
 *     - url
 *       - file path for the 'file' mode.
 *       - DB instance ... for the 'db' mode.
 * - rhoMin, rhoMax, zMin, zMax
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

    // - provider
    bool hasProvider = xmlParameter.hasAttr(_Unicode(provider));
    if (!hasProvider) {
        std::string error_msg = "[ERROR] GenericBFieldMapBrBz: Must specify the 'provider'. ";
        throw std::runtime_error(error_msg);
    }

    std::string provider = xmlParameter.attr<std::string>(_Unicode(provider));

    // - source
    bool hasSource = xmlParameter.hasChild(_Unicode(source));
    if (!hasSource) {
        std::string error_msg = "[ERROR] GenericBFieldMapBrBz: Must specify the 'source' tag. ";
        throw std::runtime_error(error_msg);

    }

    dd4hep::xml::Component source(xmlParameter.child(_Unicode(source)));
    
    // get the url
    bool hasUrl = source.hasAttr(_Unicode(url));
    if (!hasUrl) {
        std::string error_msg = "[ERROR] GenericBFieldMapBrBz: Must specify the 'url' in 'source' tag. ";
        throw std::runtime_error(error_msg);
    }

    std::string url = source.attr<std::string>(_Unicode(url));

    // 2. create the CartesianField
    dd4hep::CartesianField obj;
    GenericBFieldMapBrBz* ptr = new GenericBFieldMapBrBz();

    ptr->init_provider(provider, url);

    obj.assign(ptr, xmlParameter.nameStr(), xmlParameter.typeStr());

    return obj;

}

DECLARE_XMLELEMENT(GenericBFieldMapBrBz,create_GenericBFieldMapBrBz)
