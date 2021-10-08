// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <stdexcept>
#include <queue>
#include "Point.hpp"

template<size_t N, typename ElemType>
class KDTree {
public:
    typedef std::pair<Point<N>, ElemType> value_type;

    KDTree();

    ~KDTree();

    KDTree(const KDTree &rhs);

    KDTree &operator=(const KDTree &rhs);

    size_t dimension() const;

    size_t size() const;

    bool empty() const;

    bool contains(const Point<N> &pt) const;

    void insert(const Point<N> &pt, const ElemType &value);

    ElemType &operator[](const Point<N> &pt);

    ElemType &at(const Point<N> &pt);

    const ElemType &at(const Point<N> &pt) const;

    ElemType knn_value(const Point<N> &key, size_t k) const;

    std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

private:
    struct KDTreeNode {
        Point<N> coords;
        ElemType val;
        KDTreeNode *childrens[2];

        explicit KDTreeNode(const value_type &value) {
            coords = std::get<0>(value);
            val = std::get<1>(value);
            childrens[0] = nullptr;
            childrens[1] = nullptr;
        }
    };

    size_t dimension_;
    size_t size_;
    KDTreeNode *root = nullptr;

    bool find(Point<N> x, KDTreeNode **&p) const;
};

template<size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    // TODO(me): Fill this in.
    dimension_ = N;
    size_ = 0;
    root = nullptr;
}

template<size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    std::queue<KDTreeNode *> queueNodes;
    if (root) {
        queueNodes.push(root);
        while (queueNodes.size()) {
            KDTreeNode *x = queueNodes.front();
            if(x->childrens[0]){
                queueNodes.push(x->childrens[0]);
            }
            if(x->childrens[1]){
                queueNodes.push(x->childrens[1]);
            }
            delete x;
            queueNodes.pop();
        }
    }
}

template<size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
    dimension_ = rhs.dimension_;
    size_ = 0;
    root = nullptr;
    std::queue<KDTreeNode *> queueNodes;
    if (rhs.root) {
        queueNodes.push(rhs.root);
        while (queueNodes.size()) {
            KDTreeNode *x = queueNodes.front();
            insert(x->coords,x->val);
            if(x->childrens[0]){
                queueNodes.push(x->childrens[0]);
            }
            if(x->childrens[1]){
                queueNodes.push(x->childrens[1]);
            }
            queueNodes.pop();
        }
    }
}

template<size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
    dimension_ = rhs.dimension_;
    size_ = 0;
    root = nullptr;
    std::queue<KDTreeNode *> queueNodes;
    if (rhs.root) {
        queueNodes.push(rhs.root);
        while (queueNodes.size()) {
            KDTreeNode *x = queueNodes.front();
            insert(x->coords,x->val);
            if(x->childrens[0]){
                queueNodes.push(x->childrens[0]);
            }
            if(x->childrens[1]){
                queueNodes.push(x->childrens[1]);
            }
            queueNodes.pop();
        }
    }
    return *this;
}

template<size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return dimension_;
}

template<size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return size_;
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    return (!size_);
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(Point<N> x, KDTreeNode **&p) const {
    size_t co = 0;
    for (p = const_cast<KDTreeNode **>(&root);
         *p && !((*p)->coords == x); p = &((*p)->childrens[x[co] > (*p)->coords[co++]])) {
        if (co >= N) {
            co = 0;
        }
    }
    return (*p) != 0;
}

template<size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
    KDTreeNode **p;
    return find(pt, p);
}

template<size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
    KDTreeNode **p;
    if (find(pt, p)) {
        (*p)->val = value;
    } else {
        value_type x(pt, value);
        *p = new KDTreeNode(x);
        size_++;
    }
}

template<size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
    KDTreeNode **p;
    if (find(pt, p)) {
        return (*p)->val;
    } else {
        insert(pt, 0);
        KDTreeNode **q;
        find(pt, q);
        return (*q)->val;
    }
}

template<size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
    KDTreeNode **p;
    if (find(pt, p)) {
        return (*p)->val;
    } else {
        throw std::out_of_range("Error");
    }

}

template<size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
    KDTreeNode **p;
    if (find(pt, p)) {
        return (*p)->val;
    } else {
        throw std::out_of_range("Error");
    }
}

template<size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
    // TODO(me): Fill this in.
    ElemType new_element;
    return new_element;
}

template<size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
    // TODO(me): Fill this in.
    std::vector<ElemType> values;
    return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
