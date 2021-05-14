
#include "GenericBFieldMapBrBz.h"

GenericBFieldMapBrBz::GenericBFieldMapBrBz()
    : m_provider(nullptr) {

}

void GenericBFieldMapBrBz::fieldComponents(const double* pos, double* field) {
    double curfield[3] = {0.0, 0.0, 0.0};

    field[0] += curfield[0];
    field[1] += curfield[1];
    field[1] += curfield[2];

    return;
}
