
#ifndef KDTREE_NEIGHBORSQUEUE_HPP
#define KDTREE_NEIGHBORSQUEUE_HPP

#include <queue>
#include <map>
#include <vector>

template<typename ElemType>
class NeighborsQueue {
public:
    explicit NeighborsQueue(size_t k);

    void enqueue(ElemType val, double distance);

    ElemType value();

    std::vector<ElemType> neighbors();

    bool complete();

private:
    struct QueueNode {
        ElemType val;
        double distance;

        QueueNode(ElemType val, double distance) {
            this->val = val;
            this->distance = distance;
        }

        bool operator<(const QueueNode &u) const {
            return distance < u.distance;
        }
    };

    size_t k;
    std::priority_queue<QueueNode> Heap;

};

template<typename ElemType>
NeighborsQueue<ElemType>::NeighborsQueue(size_t k) {
    this->k = k;
}

template<typename ElemType>
void NeighborsQueue<ElemType>::enqueue(ElemType val, double distance) {
    QueueNode node(val, distance);
    Heap.push(node);
    if (Heap.size() == k + 1) {
        Heap.pop();
    }
}

template<typename ElemType>
ElemType NeighborsQueue<ElemType>::value() {
    std::map<ElemType, int> counting;
    std::vector<ElemType> vectorN = neighbors();
    for (ElemType i: vectorN) {
        counting[i] = 0;
    }
    for (ElemType i: vectorN) {
        counting[i] += 1;
    }
    int mayor = 0;
    auto iterMax = counting.begin();
    for (auto iter = counting.begin(); iter != counting.end(); iter++) {
        if (iter->second > mayor) {
            mayor = iter->second;
            iterMax = iter;
        }
    }
    return iterMax->first;
}

template<typename ElemType>
std::vector<ElemType> NeighborsQueue<ElemType>::neighbors() {
    std::vector<ElemType> neighbors;
    while (Heap.size()) {
        neighbors.push_back(Heap.top().val);
        Heap.pop();
    }
    return neighbors;
}

template<typename ElemType>
bool NeighborsQueue<ElemType>::complete() {
    return (Heap.size() == k);
}

#endif //KDTREE_NEIGHBORSQUEUE_HPP
