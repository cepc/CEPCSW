#ifndef FieldMapFileProvider_h
#define FieldMapFileProvider_h

#include "IFieldMapProvider.h"

#include <string>
#include <vector>

class FieldMapFileProvider: public IFieldMapProvider {

public:
    FieldMapFileProvider(const std::string& url_);

    // Meta data about the map
    virtual int rBinIdx(double r, double& rn);
    virtual int zBinIdx(double z, double& zn);

    // The Br and Bz
    virtual void access(int rbin, int zbin, double& Br, double& Bz);

private:
    void init();

    bool loadBr(const std::string& fn);
    bool loadBz(const std::string& fn);

    bool loadCSV(const std::string& fn, std::vector<double>& data, int& ncol, int& nrow);
private:
    std::string m_url; // the url passed from GenericBFieldMapBrBz

    std::string m_filename_Br;
    std::string m_filename_Bz;

    // the data will include the left col (z) and top row (r)
    std::vector<double> Brdata;
    std::vector<double> Bzdata;

    // in this impl, the rBin/zBin should be consistent between Br and Bz.
    int nr; // include the two endpoints, so nr = nrBin + 1
    int nz;

    double rBinMin;
    double zBinMin;

    double rBinMax;
    double zBinMax;

    double drBin;
    double dzBin;
};

#endif
