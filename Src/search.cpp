#include "search.h"
#include "node.cpp"

Search::Search()
{
//set defaults here
}

Search::~Search() {}

double min(double a, double b){
    if (a < b){
        return a;
    }
    return b;
}


void Search::heuristic(Node &n, Node end, const EnvironmentOptions &options) {
    if (options.metrictype == CN_SP_MT_EUCL){
        double x = n.i - end.i;
        double y = n.j - end.j;
        n.h = std::sqrt(x * x + y * y);
    }
    else if (options.metrictype == CN_SP_MT_DIAG){
        double dx = abs(n.i - end.i);
        double dy = abs(n.j - end.j);
        n.h = dx + dy + (std::sqrt(2) - 2) * min(dx, dy);
    }
    else if (options.metrictype == CN_SP_MT_CHEB){
        double dx = abs(n.i - end.i);
        double dy = abs(n.j - end.j);
        n.h = dx + dy - min(dx, dy);
    }
    else if (options.metrictype == CN_SP_MT_MANH){
        n.h = abs(n.i - end.i) + abs(n.j - end.j);
    }
    n.h = 0;
    n.f = n.h * 1.001 + n.g;
    //std::cout << "f: " << n.f << '\n';


}

std::vector<std::pair<std::pair<int, int>, double>> Search::Neighs(int type, Node n,const Map &map, const EnvironmentOptions &options){
    std::vector<std::pair<std::pair<int, int>, double>> res;
    std::pair<int, int> p1(n.i, n.j + 1);
    std::pair<int, int> p2(n.i, n.j - 1);
    std::pair<int, int> p3(n.i + 1, n.j);
    std::pair<int, int> p4(n.i - 1, n.j);
    std::pair<int, int> p5(n.i + 1, n.j + 1);
    std::pair<int, int> p6(n.i + 1, n.j - 1);
    std::pair<int, int> p7(n.i - 1, n.j - 1);
    std::pair<int, int> p8(n.i - 1, n.j + 1);
    if(!map.CellIsTraversable(p1)){
        //std::cout << "unpassable: " << p1.first << ' '<< p1.second << '\n';
    }
    if(!map.CellIsTraversable(p2)){
        //std::cout << "unpassable: " << p2.first << ' '<< p2.second << '\n';
    }
    if(!map.CellIsTraversable(p3)){
        //std::cout << "unpassable: " << p3.first << ' '<< p3.second << '\n';
    }
    if(!map.CellIsTraversable(p4)){
        //std::cout << "unpassable: " << p4.first << ' '<< p4.second << '\n';
    }
    if(!map.CellIsTraversable(p5)){
        //std::cout << "unpassable: " << p5.first << ' '<< p5.second << '\n';
    }
    if(!map.CellIsTraversable(p6)){
        //std::cout << "unpassable: " << p6.first << ' '<< p6.second << '\n';
    }
    if(!map.CellIsTraversable(p7)){
        //std::cout << "unpassable: " << p7.first << ' '<< p7.second << '\n';
    }
    if(!map.CellIsTraversable(p8)){
        //std::cout << "unpassable: " << p8.first << ' '<< p8.second << '\n';
    }
    if (map.CellIsTraversable(p1)){
    res.push_back(std::make_pair(p1, 1));}
    if (map.CellIsTraversable(p2)){
    res.push_back(std::make_pair(p2, 1));}
    if (map.CellIsTraversable(p3)){
    res.push_back(std::make_pair(p3, 1));}
    if (map.CellIsTraversable(p4)){
    res.push_back(std::make_pair(p4, 1));}
    if (options.allowdiagonal){
            if (options.allowsqueeze){
                if (map.CellIsTraversable(p5)){
                res.push_back(std::make_pair(p5, std::sqrt(2)));}
                if (map.CellIsTraversable(p6)){
                res.push_back(std::make_pair(p6, std::sqrt(2)));}
                if (map.CellIsTraversable(p7)){
                res.push_back(std::make_pair(p7, std::sqrt(2)));}
                if (map.CellIsTraversable(p8)){
                res.push_back(std::make_pair(p8, std::sqrt(2)));}
            }
            if (options.cutcorners && !options.allowsqueeze){
                if((map.CellIsTraversable(p1) || map.CellIsTraversable(p3)) && map.CellIsTraversable(p5)){
                    res.push_back(std::make_pair(p5, std::sqrt(2)));}
                if((map.CellIsTraversable(p1) || map.CellIsTraversable(p4)) && map.CellIsTraversable(p8)){
                    res.push_back(std::make_pair(p8, std::sqrt(2)));}
                if((map.CellIsTraversable(p2) || map.CellIsTraversable(p3)) && map.CellIsTraversable(p6)){
                    res.push_back(std::make_pair(p6, std::sqrt(2)));}
                if((map.CellIsTraversable(p2) || map.CellIsTraversable(p4)) && map.CellIsTraversable(p7)){
                    res.push_back(std::make_pair(p7, std::sqrt(2)));}

            }
            if (!options.cutcorners && !options.allowsqueeze){
                if((map.CellIsTraversable(p1) && map.CellIsTraversable(p3)) && map.CellIsTraversable(p5)){
                    res.push_back(std::make_pair(p5, std::sqrt(2)));}
                if((map.CellIsTraversable(p1) && map.CellIsTraversable(p4)) && map.CellIsTraversable(p8)){
                    res.push_back(std::make_pair(p8, std::sqrt(2)));}
                if((map.CellIsTraversable(p2) && map.CellIsTraversable(p3)) && map.CellIsTraversable(p6)){
                    res.push_back(std::make_pair(p6, std::sqrt(2)));}
                if((map.CellIsTraversable(p2) && map.CellIsTraversable(p4)) && map.CellIsTraversable(p7)){
                    res.push_back(std::make_pair(p7, std::sqrt(2)));}

            }
        }
    return res;

}

void Search::makePrimaryPath(Node curNode, std::map<std::pair<int, int>, Node> &m) {
    while (curNode.parent != std::make_pair(-1, -1)) {
        lppath.push_front(curNode);
        curNode = m[curNode.parent];
    }
}

SearchResult Search::startSearch(ILogger *Logger, const Map &map, const EnvironmentOptions &options)
{
    bool pathfound = false;
    auto t1 = std::chrono::high_resolution_clock::now();
    Node start_node(map.getStart());
    Node cur_node(map.getStart());
    Node goal_node(map.getEnd());
    cur_node.g = 0;
    heuristic(cur_node, goal_node, options);
    nodes[map.getStart()] = cur_node;
    bool done = false;
    std::pair<int, int> min_node;
    double min_f_value;
    std::vector<std::pair<std::pair<int, int>, double>> neighs;
    std::pair<int, int> new_pair;
    std::pair<int, int> cur_pair;
    int x, y;
    bool no_way_to_go = true;
    int number_of_steps = 0;
    while (!done){
        min_f_value = INFINITY;
        for (const auto& [key, value] : nodes){
            if (value.open){
                no_way_to_go = false;
                if (value.f <= min_f_value){
                    min_f_value = value.f;
                    min_node = key;
                    cur_node = value;
                }
            }
        }
        if (no_way_to_go){
            done = true;
        }
        no_way_to_go = true;
        cur_pair = min_node;

        neighs = Neighs(0, nodes[min_node],map, options);
        for (auto & neigh : neighs){
            x = neigh.first.first;
            y = neigh.first.second;
            if(!map.CellIsTraversable(x, y)){
                //std::cout << "unpassable: " << x << ' '<< y << '\n';
            }
            if (map.CellIsTraversable(x, y)){
                new_pair = std::make_pair(x, y);
                auto search = nodes.find(new_pair);
                if (search != nodes.end()){
                    if (nodes[new_pair].open){
                        if (nodes[new_pair].g > cur_node.g + neigh.second){
                            nodes[new_pair].g = cur_node.g + neigh.second;
                            //std::cout << "changing g, neigh.second = " << neigh.second <<'\n';
                            heuristic(nodes[new_pair], goal_node, options);
                            nodes[new_pair].parent = cur_pair;
                        }

                    }
                }
                else{
                    nodes[new_pair] = Node(new_pair);
                    nodes[new_pair].g = cur_node.g + neigh.second;
                    nodes[new_pair].parent = cur_pair;
                    //std::cout << "creating g, neigh.second = " << neigh.second <<'\n';


                    heuristic(nodes[new_pair], goal_node, options);
                }
            }

        }
        nodes[cur_pair].open = false;
        if (cur_pair == map.getEnd()){
            pathfound = true;
            done = true;
        }
        number_of_steps ++;
        //std::cout << "begin, end: ";
        //std::cout << map.getStart().first <<' ' << map.getStart().second << ' ' << '#' << ' ' << map.getEnd().first << ' '<< map.getEnd().second << '\n';
        //std::cout << "height, width: ";
        //std::cout << map.getMapHeight() << ' ';
        //std::cout << map.getMapWidth() << '\n';
        //std::cout << "cur: ";
        //std::cout << cur_pair.first <<' '<<cur_pair.second << '\n';
        //std::cout << "number of steps ";
        //std::cout << number_of_steps << '\n';
        //std::cout << "f: " << nodes[cur_pair].f << '\n' << "g: " << nodes[cur_pair].g << '\n' << "h: " << nodes[cur_pair].h << '\n';
        //std::cout << "distance: ";
        //std::cout << nodes[cur_pair].f << '\n';



    }
    //need to implement
    sresult.pathfound = pathfound;
    sresult.numberofsteps = number_of_steps;
    sresult.nodescreated = nodes.size();

    if (pathfound){
        makePrimaryPath(cur_node, nodes);
    }
    sresult.pathlength = cur_node.g;
    //std::cout << "pathlenght = " << cur_node.g << '\n' << "map = " << map.getCellSize() << '\n';
    auto t2 = std::chrono::high_resolution_clock::now();
    sresult.time = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
//    sresult.hppath = &hppath; //Here is a constant pointer
    //sresult.lppath = &lppath;*/
    return sresult;
}



/*void Search::makePrimaryPath(Node curNode)
{
    //need to implement
}*/

/*void Search::makeSecondaryPath()
{
    //need to implement
}*/
