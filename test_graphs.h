//
// Created by loryt on 8/21/2020.
//

#ifndef BGL_TEST_GRAPHS_H
#define BGL_TEST_GRAPHS_H

#include "graph.h"

labeled_graph_type test_graph_1()
{
    cout << "Creating test graph 1 (from notebook)" << endl;

    char name[] = "ABCD";

    labeled_graph_type g;
    add_edge(0, 1, g);
    add_edge(0, 2, g);
    add_edge(0, 3, g);
    add_edge(2, 3, g);

    /*
    for (auto v : make_iterator_range(vertices(g)))
        g[v].label = name[rand() % num_vertices(g)]; // random labels
    */

    for (size_t i = 0; i < num_vertices(g); i++)
        g[i].label = name[i];

    g[graph_bundle].m_graph_index = 1;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_2()
{
    cout << "Creating test graph 2 (from notebook)" << endl;

    char name[] = "ABCDE";

    labeled_graph_type g;
    add_edge(2, 0, g);
    add_edge(0, 3, g);
    add_edge(1, 3, g);
    add_edge(4, 3, g);
    add_edge(4, 3, g);

    /*
    for (auto v : make_iterator_range(vertices(g)))
        g[v].label = name[rand() % num_vertices(g)]; // random labels
    */

    for (size_t i = 0; i < num_vertices(g); i++)
        g[i].label = name[i];

    g[graph_bundle].m_graph_index = 2;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_3()
{
    cout << "Creating test graph 1 (from kernel survey p.20 f.5)" << endl;
    labeled_graph_type g;

    vertex_descriptor a,b,c;
    edge_descriptor e;
    bool success;

    a = add_vertex(g);
    b = add_vertex(g);
    c = add_vertex(g);

    g[a].label = "a";
    g[b].label = "b";
    g[c].label = "a";

    tie(e, success) = add_edge(a, b, g);
    assert(success);
    tie(e, success) = add_edge(b, c, g);
    assert(success);

    g[graph_bundle].m_graph_index = 1;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_4()
{
    cout << "Creating test graph 4 (from notebook)" << endl;
    string name = "ABCDEF";
    enum {A, B, C, D, E, F, n};
    labeled_graph_type g(n);

    for (int i = 0; i < n; i++)
        g[i].label = name[i];

    boost::add_edge(A, B, g);
    boost::add_edge(A, C, g);
    boost::add_edge(A, D, g);
    boost::add_edge(F, E, g);
    boost::add_edge(E, C, g);
    boost::add_edge(C, F, g);


    g[graph_bundle].m_graph_index = 4;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_5()
{
    cout << "Creating test graph 2 (from kernel survey p.20 f.5)" << endl;
    labeled_graph_type g;

    vertex_descriptor a, b1, b2, c;
    edge_descriptor e;
    bool success;

    a = add_vertex(g);
    b1 = add_vertex(g);
    b2 = add_vertex(g);
    c = add_vertex(g);

    g[a].label = "a";
    g[b1].label = "b";
    g[b2].label = "b";
    g[c].label = "c";

    tie(e, success) = add_edge(a, b1, g);
    assert(success);
    tie(e, success) = add_edge(b1, c, g);
    assert(success);
    tie(e, success) = add_edge(b1, b2, g);
    assert(success);

    g[graph_bundle].m_graph_index = 2;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_6()
{
    cout << "Creating test graph 3 (from notebook)" << endl;

    labeled_graph_type g;
    add_edge(0, 1, g);
    add_edge(0, 2, g);
    add_edge(0, 3, g);
    add_edge(1, 3, g);

    g[0].label = "A";
    g[1].label = "B";
    g[2].label = "C";
    g[3].label = "D";

    g[graph_bundle].m_graph_index = 3;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_7()
{
    cout << "Creating test graph 1 (from WL-graph-kernel p.10f.2)" << endl;

    labeled_graph_type g(6);
    add_edge(0, 1, g);
    add_edge(0, 2, g);
    add_edge(0, 3, g);
    add_edge(1, 3, g);
    add_edge(2, 3, g);
    add_edge(2, 4, g);
    add_edge(2, 5, g);
    //*/
    g[0].label = "5";
    g[1].label = "2";
    g[2].label = "4";
    g[3].label = "3";
    g[4].label = "1";
    g[5].label = "1";
    /*/
    g[0].label = "green";
    g[1].label = "pink";
    g[2].label = "grey";
    g[3].label = "yellow";
    g[4].label = "red";
    g[5].label = "red";
    **/

    g[graph_bundle].m_graph_index = 1;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_8()
{
    cout << "Creating test graph 2 (from WL-graph-kernel p.10f.2)" << endl;

    labeled_graph_type g(6);
    add_edge(0, 1, g);
    add_edge(0, 2, g);
    add_edge(1, 2, g);
    add_edge(1, 3, g);
    add_edge(2, 3, g);
    add_edge(2, 4, g);
    add_edge(3, 5, g);
    //*/
    g[0].label = "2";
    g[1].label = "5";
    g[2].label = "4";
    g[3].label = "3";
    g[4].label = "1";
    g[5].label = "2";
    /*/
    g[0].label = "pink";
    g[1].label = "green";
    g[2].label = "grey";
    g[3].label = "yellow";
    g[4].label = "red";
    g[5].label = "pink";
    **/

    g[graph_bundle].m_graph_index = 2;
    g[graph_bundle].m_density = density(g);

    return g;
}

#endif //BGL_TEST_GRAPHS_H
