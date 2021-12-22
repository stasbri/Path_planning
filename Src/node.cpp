#include "node.h"


Node::Node(std::pair<int, int> p){
    i = p.first;
    j = p.second;
    f = 0;
    g = 0;
    h = 0;
    parent = std::make_pair(-1, -1);
}


