#include "FieldMapFileProvider.h"

#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

FieldMapFileProvider::FieldMapFileProvider(const std::string& url_)
    : m_url(url_), nr(-1), nz(-1), 
      rBinMin(-1), zBinMin(-1), 
      rBinMax(-1), zBinMax(-1), 
      drBin(-1), dzBin(-1) {
    init();
}

int FieldMapFileProvider::rBinIdx(double r, double& rn) {
    // 
    // | --- | --- | --- |
    // ^        ^        ^
    // rmin     |        rmax
    //          r
    //       | 
    //    return


    double absr = std::fabs(r);

    int idx = -1;

    // std::cout << "FieldMapFileProvider::rBinIdx: "
    //           << " r: " << r
    //           << " rBinMin: " << rBinMin
    //           << " rBinMax: " << rBinMax
    //           << std::endl;

    if ( rBinMin <= absr && absr < rBinMax) {
        idx = (absr - rBinMin) / drBin;
        double r0 = rBinMin + idx*drBin;
        rn = (absr - r0)/drBin;
    }
    
    return idx;
}

int FieldMapFileProvider::zBinIdx(double z, double& zn) {
    double absz = std::fabs(z);

    int idx = -1;

    if ( zBinMin <= absz && absz < zBinMax) {
        idx = (absz - zBinMin) / dzBin;
        double z0 = zBinMin + idx*dzBin;
        zn = (absz - z0)/dzBin;
    }
    
    return idx;
}

void FieldMapFileProvider::access(int rbin, int zbin, double& Br, double& Bz) {
    // the valid bin should between [0, n)
    // if the point is not in the valid region, return 0
    if ((rbin < 0 || rbin >= nr) || (zbin < 0 || zbin >= nr)) {
        Br = 0;
        Bz = 0;
        return;
    }

    // convert to the internal table (with left col and top row)
    int ridx = rbin+1;
    int zidx = zbin+1;

    int globalidx = ridx + zidx*(nr+1);

    Br = Brdata[globalidx];
    Bz = Bzdata[globalidx];
}

// ======================
// Below are private impl
// ======================

void FieldMapFileProvider::init() {
    // parse the url
    //   Br=/tmp/lint/CEPCSW/Br.csv;Bz=/tmp/lint/CEPCSW/Bz.csv

    std::map<std::string, std::string> url_parsed;
    
    size_t idx_begin = 0;
    size_t idx_end = m_url.find(";");

    while (true) {
        std::string keyval = m_url.substr(idx_begin, idx_end-idx_begin);

        if (keyval.size()) {
            std::cout << "---> keyval: " << keyval << std::endl;
            // separate by '='
            size_t idx_sep = keyval.find("=");
            if (idx_sep == std::string::npos) {
                std::string error_msg = "[ERROR] FieldMapFileProvider: Pleaes specify label=path. ";
                throw std::runtime_error(error_msg);

            }
            std::string key = keyval.substr(0, idx_sep);
            std::string val = keyval.substr(idx_sep+1);
            std::cout << "-----> key: " << key << std::endl;
            std::cout << "-----> val: " << val << std::endl;
            // insert into map

            // duplicated
            if (url_parsed.count(key)) {
                std::string error_msg = "[ERROR] FieldMapFileProvider: duplicated key '" + key + "'";
                throw std::runtime_error(error_msg);
            }

            url_parsed[key] = val;
        }

        if (idx_end == std::string::npos) {
            break;
        }

        idx_begin = idx_end+1;
        idx_end = m_url.find(';', idx_begin);
    }

    // access the map 'url_parsed'
    std::string key;

    // load Br
    key = "Br";
    if (url_parsed.count(key) == 0) {
        std::string error_msg = "[ERROR] FieldMapFileProvider: missing key '" + key + "'";
        throw std::runtime_error(error_msg);
    }
    loadBr(url_parsed[key]);

    // load Bz
    key = "Bz";
    if (url_parsed.count(key) == 0) {
        std::string error_msg = "[ERROR] FieldMapFileProvider: missing key '" + key + "'";
        throw std::runtime_error(error_msg);
    }
    loadBz(url_parsed[key]);

}



bool FieldMapFileProvider::loadBr(const std::string& fn) {

    //  ----> r
    // |
    // v
    // z
    // the shape
    int ncol = -1;
    int nrow = -1;

    // load data
    bool st = loadCSV(fn, Brdata, ncol, nrow);

    // update the metadata
    // for (int i = 0; i < ncol; ++i) {
    //     std::cout << "i: " << i << " data[i]: " << data[i] << std::endl;
    // }
    nr = ncol-1; // skip the top row and left col
    nz = nrow-1; 

    rBinMin = Brdata[1];
    zBinMin = Brdata[ncol];

    rBinMax = Brdata[ncol-1];
    zBinMax = Brdata[ncol*(nrow-1)];

    drBin = (rBinMax-rBinMin) / (nr-1);
    dzBin = (zBinMax-zBinMin) / (nz-1);

    std::cout << "nr: " << nr << std::endl;
    std::cout << "nz: " << nz << std::endl;
    std::cout << "rBinMin: " << rBinMin << std::endl;
    std::cout << "zBinMin: " << zBinMin << std::endl;
    std::cout << "rBinMax: " << rBinMax << std::endl;
    std::cout << "zBinMax: " << zBinMax << std::endl;
    std::cout << "drBin: " << drBin << std::endl;
    std::cout << "dzBin: " << dzBin << std::endl;



    return true;
}

bool FieldMapFileProvider::loadBz(const std::string& fn) {

    // the shape
    int ncol = -1;
    int nrow = -1;
    bool st = loadCSV(fn, Bzdata, ncol, nrow);

    return true;
}

bool FieldMapFileProvider::loadCSV(const std::string& fn, 
                                   std::vector<double>& data, 
                                   int& ncol, int& nrow) {

    std::ifstream input(fn);
    std::string tmpline;

    ncol = 0;
    nrow = 0;


    while(std::getline(input, tmpline)) {
        std::cout << "------> " << tmpline << std::endl;
        ++nrow;

        // split by ,
        int col_per_line = 0;
        
        size_t idx_begin = 0;
        size_t idx_end = tmpline.find(",");

        while (true) {
            std::cout << "----------------> idx: "
                      << idx_begin << "->" << idx_end << std::endl;
            std::string val = tmpline.substr(idx_begin, idx_end-idx_begin);
            double v = 0.0;
            // parse val
            std::stringstream ss;
            ss << val;
            ss >> v;
            // if (not ss.good()) {
            //     v = 0;
            //     ss.clear();
            // }

            std::cout << "----------------> v: " << v << std::endl;
            // if (val.size()==0 || ) { // if it is empty or not number, use 0
            //     data.push_back(0.0);
            // } else {
            // }
            data.push_back(v);

            ++col_per_line;

            // break
            if (idx_end == std::string::npos) {
                break;
            }

            // goto next val
            idx_begin = idx_end+1;
            idx_end = tmpline.find(',', idx_begin);

        }

        // if it is the first line, we need to store the column sizes
        if (ncol == 0) {
            ncol = col_per_line;
        }

        std::cout << "--> ncol: " << ncol
                  << " col_per_line: " << col_per_line
                  << std::endl;

        if (ncol != col_per_line) {
            std::string error_msg = "[ERROR] FieldMapFileProvider: Mismatch columns. ";
            throw std::runtime_error(error_msg);
        }

    }

    std::cout << "SUMMARY: "
              << " size of the data: " << data.size()
              << " shape (ncol: " << ncol
              << " , nrow: " << nrow
              << ")" << std::endl;

    return true;
}
