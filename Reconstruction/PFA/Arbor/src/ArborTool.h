#ifndef ARBORTOOL_H_
#define ARBORTOOL_H_

#include "TVector3.h"
#include "TMatrixF.h"
#include <iostream>
#include <vector>
#include <string>


class QuickUnion{
  std::vector<unsigned> _id;
  std::vector<unsigned> _size;
  int _count;

  public:
    QuickUnion(const unsigned NBranches) {
      _count = NBranches;
      _id.resize(NBranches);
      _size.resize(NBranches);
      for( unsigned i = 0; i < NBranches; ++i ) {
        _id[i] = i;
        _size[i] = 1;
      }
    }

    int count() const { return _count; }

    unsigned find(unsigned p) {
      while( p != _id[p] ) {
        _id[p] = _id[_id[p]];
        p = _id[p];
      }
      return p;
    }

    bool connected(unsigned p, unsigned q) { return find(p) == find(q); }

    void unite(unsigned p, unsigned q) {
      unsigned rootP = find(p);
      unsigned rootQ = find(q);
      _id[p] = q;

      if(_size[rootP] < _size[rootQ] ) {
        _id[rootP] = rootQ; _size[rootQ] += _size[rootP];
      } else {
        _id[rootQ] = rootP; _size[rootP] += _size[rootQ];
      }
      --_count;
    }
  };

typedef std::vector< std::vector<int> > branchcoll;
typedef std::vector<int> branch;
typedef std::vector< std::pair<int, int> > linkcoll;


TMatrixF MatrixSummarize( TMatrixF inputMatrix );		//ArborCoreNeed

std::vector<int>SortMeasure( std::vector<float> Measure, int ControlOrderFlag );	//ArborCoreNeed

branchcoll ArborBranchMerge(branchcoll inputbranches, TMatrixF inputMatrix);		//ArborCoreNeed

#endif //
