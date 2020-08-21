//
// Created by loryt on 8/21/2020.
//

#ifndef BGL_TEST_GRAPHS_H
#define BGL_TEST_GRAPHS_H



using namespace boost;
using namespace std;
/* STRUCTS FOR KERNEL FEATURES */
typedef vector<vector<size_t>> distance_matrix_t;
typedef vector<vector<size_t>> hash_table_t;
/* PROPERTIES DEFINITIONS */
// simplest properties
struct m_graph_property
{
    m_graph_property() : m_graph_index(0), m_density(0.0) {}
    // used to identify a graph
    size_t m_graph_index;
    // used to store graph density
    double m_density;
};

struct m_edge_property
{
    // weight is positive integer
    size_t weight;
    m_edge_property() : weight(1) {}
    explicit m_edge_property(size_t w) : weight(w) {}
};

struct m_vertex_property
{
    string label;
    m_vertex_property() : label("default") {}
    explicit m_vertex_property(string l) : label(std::move(l)) {}
};

/* GRAPH TYPE DEFINITIONS */
// graph type identified by an index, edges have a weight and vertices have a name
typedef adjacency_list<
        setS, // store edges in set to avoid multi-edges, descriptor is integer and can be used as offset in external properties
        vecS, // store vertices
        undirectedS, // no direction
        m_vertex_property, // label vertices with words or letters or numbers
        m_edge_property, // natural weight on edges
        m_graph_property> // distinguish graphs by index
labeled_graph_type;
/* USE SPECIFIC TYPES    */
using vertex_descriptor = graph_traits<labeled_graph_type>::vertex_descriptor;
using edge_descriptor = graph_traits<labeled_graph_type>::edge_descriptor;
using vertex_iterator = graph_traits<labeled_graph_type>::vertex_iterator;
using edge_iterator = graph_traits<labeled_graph_type>::edge_iterator;

typedef vector<labeled_graph_type> graph_vector;

/* FUNCTION DEFINITIONS */
ostream& operator<<(ostream& os, const std::pair<string, string>& p)
{
    os << p.first << " -> " << p.second;
    return os;
}

template<class tuple_type, size_t... I>
void print_tuple(const tuple_type& _tup, std::index_sequence<I...>)
{
    std::cout << "(";
    (..., (std::cout << (I == 0 ? "" : ", ") << std::get<I>(_tup)));
    std::cout << ")\n";
}

template<class... T>
void print_tuple(const std::tuple<T...>& _tup)
{
    print_tuple(_tup, std::make_index_sequence<sizeof...(T)>());
}

inline double density(const labeled_graph_type& g)
{
    try
    {
        return  (2 * (double)num_edges(g)) / ((double)num_vertices(g) * ((double)num_vertices(g) - 1));
    }
    catch (const std::exception&)
    {
        return 0.0;
    }

}

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
    cout << "Creating test graph 1 (from kernelsurvey p.23f.7)" << endl;

    labeled_graph_type g(6);
    add_edge(0, 1, g);
    add_edge(0, 2, g);
    add_edge(0, 3, g);
    add_edge(1, 3, g);
    add_edge(2, 3, g);
    add_edge(2, 4, g);
    add_edge(2, 5, g);
    /*/
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
    /**/

    g[graph_bundle].m_graph_index = 1;
    g[graph_bundle].m_density = density(g);

    return g;
}
labeled_graph_type test_graph_8()
{
    cout << "Creating test graph 2 (from kernelsurvey p.23f.7)" << endl;

    labeled_graph_type g(6);
    add_edge(0, 1, g);
    add_edge(0, 2, g);
    add_edge(0, 3, g);
    add_edge(1, 3, g);
    add_edge(2, 3, g);
    add_edge(2, 4, g);
    add_edge(2, 5, g);
    /*/
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
    /**/

    g[graph_bundle].m_graph_index = 2;
    g[graph_bundle].m_density = density(g);

    return g;
}



#endif //BGL_TEST_GRAPHS_H
