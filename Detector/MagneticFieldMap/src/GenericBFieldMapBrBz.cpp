
#include "GenericBFieldMapBrBz.h"

#include "FieldMapFileProvider.h"

GenericBFieldMapBrBz::GenericBFieldMapBrBz()
    : m_provider(nullptr) {
    type = dd4hep::CartesianField::MAGNETIC;

}

void GenericBFieldMapBrBz::fieldComponents(const double* pos, double* field) {
    double curfield[3] = {0.0, 0.0, 0.0};

    field[0] += curfield[0];
    field[1] += curfield[1];
    field[1] += curfield[2];

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
