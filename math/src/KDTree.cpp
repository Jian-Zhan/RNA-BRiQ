/*
 * KDTree.cpp
 *
 *  Created on: Oct 30, 2018
 *      Author: s2982206
 */

#include "math/KDTree.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>


namespace NSPmath {


KDNode::KDNode(){
	this->index = 0;
	this->dim = 0;
	this->left = NULL;
	this->right = NULL;
	this->level = 0;
	this->x = new double[0];
}

KDNode::KDNode(double* x, int index, int dim, KDNode* left, KDNode* right) {
    this->index = index;
    this->dim = dim;
    this->x = new double[dim];
    for(int i=0;i<dim;i++) {
    	this->x[i] = x[i];
    }
    this->left = left;
    this->right = right;
    this->level = 0;
}

KDNode::KDNode(const pointIndex &pi, int dim, KDNode* left, KDNode* right, int level) {
    this->index = pi.second;
    this->dim = dim;
    this->x = new double[dim];
    for(int i=0;i<dim;i++)
    	this->x[i] = pi.first[i];
    this->left = left;
    this->right = right;
    this->level = level;
}

KDNode::~KDNode(){
	delete x;
}

double KDNode::coord(const int &idx) { return x[idx]; }
KDNode::operator bool() { return (this->dim == 0); }
KDNode::operator int() { return index; }
KDNode::operator pointIndex() { return pointIndex(x, index); }



comparer::comparer(int idx_) : idx{idx_} {};

inline bool comparer::compare_idx(const pointIndex &a,  //
                                  const pointIndex &b   //
) {
    return (a.first[idx] < b.first[idx]);  //
}

inline void sort_on_idx(const pointIndexArr::iterator &begin,  //
                        const pointIndexArr::iterator &end,    //
                        int idx) {
    comparer comp(idx);
    comp.idx = idx;

    using std::placeholders::_1;
    using std::placeholders::_2;

    std::sort(begin, end, std::bind(&comparer::compare_idx, comp, _1, _2));
}

using pointVec = std::vector<double*>;

int KDTree::selectLevel(const pointIndexArr::iterator &begin,
        const pointIndexArr::iterator &end){
	if(begin == end){
		return 0; //empty tree
	}
	double mean[dim];
	double sd[dim];
	double maxSD = 0;
	int highestSDLevel = 0;
	for(int i=0;i<dim;i++){
		double tot = 0;
		int n = 0;
		for(auto it=begin;it<end;++it){
			tot += it->first[i];
			n ++;
		}
		mean[i] = tot/n;
		double sdTot = 0;
		for(auto it=begin;it<end;++it){
			sdTot += (it->first[i] - mean[i])*(it->first[i] - mean[i]);
		}
		sd[i] = sdTot/n;
		if(sd[i] > maxSD){
			maxSD = sd[i];
			highestSDLevel = i;
		}
	}
	return highestSDLevel;
}

KDNode* KDTree::make_tree(const pointIndexArr::iterator &begin,  //
                            const pointIndexArr::iterator &end,    //
                            int length) {
    if (begin == end) {
        return NULL;  // empty tree
    }

    int level = selectLevel(begin, end);

    if (length > 1) {
        sort_on_idx(begin, end, level);
    }

    auto middle = begin + (length / 2);
    auto l_begin = begin;
    auto l_end = middle;
    auto r_begin = middle + 1;
    auto r_end = end;

    int l_len = length / 2;
    int r_len = length - l_len - 1;

    KDNode* left;
    if (l_len > 0 && dim > 0) {
        left = make_tree(l_begin, l_end, l_len);
    } else {
        left = leaf;
    }
    KDNode* right;
    if (r_len > 0 && dim > 0) {
        right = make_tree(r_begin, r_end, r_len);
    } else {
        right = leaf;
    }

    // KDNode result = KDNode();
    KDNode* node = new KDNode(*middle, dim, left, right, level);
    nodeList.push_back(node);
    return node;
}

KDTree::KDTree(pointVec& point_array, int dim) {
    leaf = NULL;
    // iterators
    this->dim = dim;

    pointIndexArr arr;
    for(int i=0;i<point_array.size();i++){
    	arr.push_back(pointIndex(point_array[i], i));
    }
    auto begin = arr.begin();
    auto end = arr.end();
    int length = arr.size();
    root = KDTree::make_tree(begin, end, length);
}

KDNode* KDTree::nearest_(   //
    KDNode* branch,  //
	double* pt,        //
    KDNode* best,    //
    double best_dist   //
) {
	//return nodeList[0];
	if (branch==NULL) {
	        return NULL;  // basically, null
	}


	double d, dx, dx2;
    int dim = branch->dim;
    int level = branch->level;

    d = dist(branch->x, pt, dim);
    dx = branch->x[level] - pt[level];
    dx2 = dx * dx;

    KDNode* best_l = best;
    double best_dist_l = best_dist;

    if (d < best_dist) {
        best_dist_l = d;
        best_l = branch;
    }

    KDNode* section;
    KDNode* other;

    // select which branch makes sense to check
    if (dx > 0) {
        section = branch->left;
        other = branch->right;
    } else {
        section = branch->right;
        other = branch->left;
    }

    // keep nearest neighbor from further down the tree
    KDNode* further = nearest_(section, pt, best_l, best_dist_l);
    if (further!=NULL) {
        double dl = dist(further->x, pt,dim);
        if (dl < best_dist_l) {
            best_dist_l = dl;
            best_l = further;
        }
        // only check the other branch if it makes sense to do so
        if (dx2 < best_dist_l) {
            further = nearest_(other, pt, best_l, best_dist_l);
            if (further!=NULL) {
                dl = dist(further->x, pt, dim);
                if (dl < best_dist_l) {
                    best_dist_l = dl;
                    best_l = further;
                }
            }
        }
    }
    return best_l;
};


// default caller


KDNode* KDTree::nearest_(double* pt) {
    // KDNodePtr best = branch;
	double branch_dist = dist(root->x, pt, root->dim);
    return nearest_(root,          // beginning of tree
                    pt,            // point we are querying
                    root,          // best is the root
                    branch_dist);  // best_dist = branch_dist
};




int KDTree::nearest_index(double* pt) {
    return nearest_(pt)->index;
};


void KDTree::printTree(){
	printTree(this->root);
}

void KDTree::printTree(KDNode* root){
	if(root == NULL) return;
	int rootIndex, leftIndex, rightIndex;
	if(root->left == NULL)
		leftIndex = -1;
	else
		leftIndex = root->left->index;

	if(root->right == NULL)
		rightIndex = -1;
	else
		rightIndex = root->right->index;

	rootIndex = root->index;

	printf("%d %d %d",rootIndex, leftIndex, rightIndex);
	for(int i=0;i<dim;i++) {
		printf(" %4.2f",root->x[i]);
	}
	printf("\n");
	printTree(root->left);
	printTree(root->right);
}

KDTree::~KDTree(){
	for(int i=0;i<nodeList.size();i++){
		delete nodeList[i];
	}
}


} /* namespace NSPmodel */
