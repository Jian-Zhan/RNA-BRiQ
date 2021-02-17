/*
 * KDTree.h
 *
 *  Created on: Oct 30, 2018
 *      Author: s2982206
 */

#ifndef MATH_KDTREE_H_
#define MATH_KDTREE_H_

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

namespace NSPmath {

/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 *  https://rosettacode.org/wiki/K-d_tree
 * It is a reimplementation of the C code using C++.
 * It also includes a few more queries than the original
 *
 */

using namespace std;
//using point_t = std::vector<double>;
//using indexArr = std::vector<int>;

using pointIndex = typename std::pair<double*, int >;

class KDNode {
   public:
    //using KDNodePtr = std::shared_ptr<KDNode>;
    int index;
    int dim;
    double* x;
    KDNode* left;
    KDNode* right;
    int level; //highest sd dimension for partition

    // initializer
    KDNode();
    KDNode(double* x, int index, int dim, KDNode* left, KDNode* right);
    KDNode(const pointIndex &pi, int dim, KDNode* left, KDNode* right, int level);
    ~KDNode();

    // getter
    double coord(const int &);

    // conversions
    explicit operator bool();
    explicit operator int();
    explicit operator pointIndex();
};

//using KDNodePtr = std::shared_ptr< KDNode >;

//KDNodePtr NewKDNodePtr();

inline double dist(double* a, double* b, int dim);


// Need for sorting
class comparer {
   public:
    int idx;
    explicit comparer(int idx_);
    inline bool compare_idx(
    		const pointIndex &a,  //
			const pointIndex &b   //
    );
};

using pointIndexArr = typename std::vector< pointIndex >;

inline void sort_on_idx(const pointIndexArr::iterator &,  //
                        const pointIndexArr::iterator &,  //
                        int idx);

using pointVec = std::vector<double*>;

class KDTree {
    KDNode* root;
    KDNode* leaf;
    int dim;
    vector<KDNode*> nodeList;

    int selectLevel(const pointIndexArr::iterator &begin,
    		        const pointIndexArr::iterator &end);

    KDNode* make_tree(const pointIndexArr::iterator &begin,  //
                        const pointIndexArr::iterator &end,    //
                        int length);

   public:
    	KDTree(){this->root = NULL; this->leaf = NULL; this->dim = 0;}
    	KDTree(pointVec& point_array, int dim);

   public:
    KDNode* nearest_(           //
        KDNode* branch,  //
		double* pt,        //
        KDNode* best,    //
        double best_dist   //
    );

    // default caller

      KDNode* nearest_(double*);
      int nearest_index(double*);
      void printTree();
      void printTree(KDNode* root);
      ~KDTree();

};


inline double dist(double* a, double* b, int dim) {
    double distc = 0;
    for (int i = 0; i < dim; i++) {
        double di = a[i] - b[i];
        distc += di * di;
    }
    return distc;
}


} /* namespace NSPmodel */

#endif /* MATH_KDTREE_H_ */
