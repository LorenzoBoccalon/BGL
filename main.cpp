#include <iostream>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <functional>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/bimap.hpp>
#include <boost/container_hash/hash.hpp>
#include <utility>


using namespace boost;
using namespace std;
/* STRUCTS FOR KERNEL FEATURES */
typedef vector<vector<size_t>> distance_matrix_t;

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

/* FUNCTIONS */
// visualize graph vertices, edges and density.
void print_graph_info(const labeled_graph_type& g);
// calc graph density. It's used to decide which algo to use later.
inline double density(const labeled_graph_type& g);
// method loads n graphs from file "shock.txt" into G_
void load_graphs(graph_vector& G_, size_t n = 2, const string& filename = "data/shock.txt");
// method prints all graphs info
void print_graph_vector_info(const graph_vector& G_);
// func to print all edge weights
void print_weights(const labeled_graph_type& g);
// func to print all vertex names
void print_vertex_labels(const labeled_graph_type& g);
// Return matrix with all pair shortest paths
distance_matrix_t get_all_pairs_shortest_paths(const labeled_graph_type& g);
// Print all pair shortest paths. If verbose is true then shows matrix of paths
template<typename VertexNameMap>
void print_distance_matrix(distance_matrix_t distance_matrix, VertexNameMap name_map);
// create a copy of the graph with an edge connecting each vertex with weight wq  equal to shortest path
labeled_graph_type floyd_warshall_transform(const labeled_graph_type& g, bool verbose);
// compute the shortest path kernel between two graphs
size_t shortest_path_kernel(const labeled_graph_type& S_1, const labeled_graph_type& S_2);
// compute the weisfeiler lehman kernel between two graphs
size_t weisfeiler_lehman_kernel(const labeled_graph_type& g_1, const labeled_graph_type& g_2, size_t depth);
// Create test graphs
labeled_graph_type test_graph_1();
labeled_graph_type test_graph_2();
labeled_graph_type test_graph_3();
labeled_graph_type test_graph_4();
labeled_graph_type test_graph_5();
labeled_graph_type test_graph_6();
labeled_graph_type test_graph_7();
labeled_graph_type test_graph_8();

int main()
{
    //cout << "how many graphs do you need?\t> ";
    //cin >> n;

    // vector to store graphs
    graph_vector G;

    load_graphs(G);
    print_graph_vector_info(G);

    const auto& g_1 = G[0];
    const auto& g_2 = G[1];

    auto test_1 = floyd_warshall_transform(g_1, false);
    print_graph_info(test_1);

    auto test_2 = floyd_warshall_transform(g_2, false);
    print_graph_info(test_2);

    cout << endl << endl << "\tShortest Path Kernel" << endl << endl;
    cout << "computed on (g_1, g_2)" << endl << shortest_path_kernel(test_1, test_2);

    size_t depth = 2;
    cout << endl << endl << "\tWeisfeiler Lehman Kernel" << endl << endl;
    cout << "computed on (g_1, g_2), at depth = " << depth << endl << weisfeiler_lehman_kernel(g_1, g_2, depth);

    return 0;
}


void print_vertex_labels(const labeled_graph_type& g)
{
    cout << "vertices(g) = { ";
    for (auto v : make_iterator_range(vertices(g)))
        std::cout << g[v].label << " ";
    cout << "}" << endl;
}

void print_weights(const labeled_graph_type& g)
{
    for (const auto& e : make_iterator_range(edges(g)))
        cout << g[e.m_source].label << " -> " << g[e.m_target].label << " : " << g[e].weight << endl;
}

distance_matrix_t get_all_pairs_shortest_paths(const labeled_graph_type& g)
{
    size_t dim = num_vertices(g);

    // Initialize empty matrix with #dim vectors of size #dim
    distance_matrix_t distance_matrix(dim, vector<size_t>(dim));
    auto weights = boost::get(&m_edge_property::weight, g);
    bool done = false;
    // use optimal algorithm based on graph density
    if (g[graph_bundle].m_density < 0.7)
    {
        // DEBUG
        cout << "using johnson algo..." << endl;
        done = boost::johnson_all_pairs_shortest_paths(g, distance_matrix, boost::weight_map(weights));
    }
    else
    {
        // DEBUG
        cout << "using floyd warshall algo..." << endl;
        done = boost::floyd_warshall_all_pairs_shortest_paths(g, distance_matrix, boost::weight_map(weights));
    }
    // Assert if something's wrong
    assert(done);
    return distance_matrix;
}

template<typename VertexNameMap>
void print_distance_matrix(distance_matrix_t distance_matrix, VertexNameMap name_map)
{
    auto dim = distance_matrix.size();
    cout << "Distance Matrix of g\ndim = " << dim << " x " << dim << endl;
    for (size_t i = 0; i < dim; i++)
        cout << name_map[i] << " ";
    cout << endl;
    std::for_each(distance_matrix.begin(), distance_matrix.end(), [](const auto& row)
    {
        std::for_each(row.begin(), row.end(), [](size_t w)
        {
            cout << w << " ";
        });
        cout << endl;
    });
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


labeled_graph_type floyd_warshall_transform(const labeled_graph_type& g, bool verbose)
{
    distance_matrix_t A = get_all_pairs_shortest_paths(g);
    labeled_graph_type S(num_vertices(g));
    // DEBUG
    if (verbose)
    {
        auto name_map = get(&m_vertex_property::label, g);
        print_distance_matrix(A, name_map);
        cout << "names:" << endl;
        for (size_t i = 0; i < A.size(); i++)
        {
            S[i].label = g[i].label;
            cout << S[i].label << " ";
        }
        cout << endl;
        cout << "weights:" << endl;
        for (size_t row = 0; row < A.size(); ++row) {
            for (size_t col = row + 1; col < A.size(); ++col) {
                if (row != col) {
                    edge_descriptor e; bool s;
                    tie(e, s) = add_edge(row, col, S);
                    S[e].weight = A[row][col];
                    cout << S[e].weight << " ";
                }
            }
        }
        cout << endl;
    }
    else
    {
        // copy the labels
        for (size_t i = 0; i < A.size(); i++)
            S[i].label = g[i].label;
        // create an edge for each pair of vertices
        for (size_t row = 0; row < A.size(); ++row) {
            for (size_t col = row + 1; col < A.size(); ++col) {
                if (row != col) {
                    edge_descriptor e; bool s;
                    tie(e, s) = add_edge(row, col, S);
                    S[e].weight = A[row][col];
                }
            }
        }
    }

    S[graph_bundle].m_density = density(S); // should be 1
    S[graph_bundle].m_graph_index = 1234;   // meaningless ID
    return S;
}

size_t shortest_path_kernel(const labeled_graph_type& S_1, const labeled_graph_type& S_2)
{
    size_t kernel = 0;

    std::set<std::tuple<string, string, size_t>> tuples;

    // insert paths in both directions
    for (const auto& e_1 : make_iterator_range(edges(S_1)))
    {
        tuples.insert(std::make_tuple(S_1[e_1.m_source].label, S_1[e_1.m_target].label, S_1[e_1].weight));
        tuples.insert(std::make_tuple(S_1[e_1.m_target].label, S_1[e_1.m_source].label, S_1[e_1].weight));
    }
    
    for (const auto& e_2 : make_iterator_range(edges(S_2)))
    {
        tuples.insert(std::make_tuple(S_2[e_2.m_source].label, S_2[e_2.m_target].label, S_2[e_2].weight));
        tuples.insert(std::make_tuple(S_2[e_2.m_target].label, S_2[e_2.m_source].label, S_2[e_2].weight));
    }

    cout << "tuples set:" << endl;
    for (const auto& t : tuples)
        print_tuple(t);

    vector<size_t> phi_S_1(tuples.size(), 0);
    vector<size_t> phi_S_2(tuples.size(), 0);

    // for both graphs, count how many paths are present
    for (const auto& e_1 : make_iterator_range(edges(S_1)))
    {
        auto t1 = std::make_tuple(S_1[e_1.m_source].label, S_1[e_1.m_target].label, S_1[e_1].weight);
        auto t2 = std::make_tuple(S_1[e_1.m_target].label, S_1[e_1.m_source].label, S_1[e_1].weight);

        auto f1 = tuples.find(t1);
        auto f2 = tuples.find(t2);

        if (f1 != tuples.end())
        {
            auto pos = std::distance(tuples.begin(), f1);
            phi_S_1[pos]++;
        }
        if (f2 != tuples.end())
        {
            auto pos = std::distance(tuples.begin(), f2);
            phi_S_1[pos]++;
        }
    }

    cout << "Phi(S_1):" << endl;
    for (const auto& n : phi_S_1)
        cout << n << " ";
    cout << endl;

    for (const auto& e_2 : make_iterator_range(edges(S_2)))
    {
        auto t1 = std::make_tuple(S_2[e_2.m_source].label, S_2[e_2.m_target].label, S_2[e_2].weight);
        auto t2 = std::make_tuple(S_2[e_2.m_target].label, S_2[e_2.m_source].label, S_2[e_2].weight);

        auto f1 = tuples.find(t1);
        auto f2 = tuples.find(t2);

        if (f1 != tuples.end())
        {
            auto pos = std::distance(tuples.begin(), f1);
            phi_S_2[pos]++;
        }
        if (f2 != tuples.end())
        {
            auto pos = std::distance(tuples.begin(), f2);
            phi_S_2[pos]++;
        }
    }

    cout << "Phi(S_2):" << endl;
    for (const auto& n : phi_S_2)
        cout << n << " ";
    cout << endl;

    // result is inner product of the two feature vectors Phi
    kernel = std::inner_product(phi_S_1.begin(), phi_S_1.end(), phi_S_2.begin(), 0);

    cout << "result = ";
    return kernel;
}

size_t helper_insert(bimap<size_t, string>& label_id_map, const string& the_label)
{
    // hash label to create unique identifier
    auto hash = hash_value(the_label);

    typedef bimap<size_t, string> bm_t;
    // add mapping hash-label to set
    label_id_map.insert(bm_t::value_type(hash, the_label));

    return hash;
}

size_t helper_insert(bimap<size_t, string>& label_id_map, const vector<size_t>& the_label)
{
    // hash label to create unique identifier
    auto hash = hash_range(the_label.begin(), the_label.end());

    string new_label_string = label_id_map.left.find(the_label[0])->second + ", ";
    for(int i = 1; i < the_label.size(); i++)
        new_label_string += label_id_map.left.find(the_label[i])->second + " ";

    typedef bimap<size_t, string> bm_t;
    // add mapping hash-label to set
    label_id_map.insert(bm_t::value_type(hash, new_label_string));

    return hash;
}

void wl_initialize(const labeled_graph_type& g, vector<vector<vertex_descriptor>>& ver_to_id, bimap<size_t, string>& label_id_map)
{
    typedef bimap<size_t, string> bm_t;
    // find all unique labels and add them to level 0
    for (const auto& v : make_iterator_range(vertices(g)))
    {
        size_t level = 0;
        // read current label
        auto the_label = g[v].label;
        // check if present
        auto present = label_id_map.right.find(the_label);
        // present        : r_iter
        // present.first  : string
        // present.second : size_t

        // if NOT present map it to a number and insert
        if ( present == label_id_map.right.end())
        {
            // add it to set and
            // insert id in vertex id map
            ver_to_id[level][v] = helper_insert(label_id_map, the_label);
            // DEBUG
            cout << "not present | ";
        }
        else
        {
            // insert id in vertex id map
            ver_to_id[level][v] = present->second;
            // DEBUG
            cout << "present | ";
        }
        // DEBUG
        cout << v << " -> " << label_id_map.left.find(ver_to_id[level][v])->second << endl;
    }
}

size_t weisfeiler_lehman_kernel(const labeled_graph_type& g_1, const labeled_graph_type& g_2, size_t depth)
{
    size_t kernel = 0;

    // https://www.boost.org/doc/libs/1_74_0/doc/html/hash/reference.html#boost.hash_combine

    // map integer to a label
    // bm.left = num -> string
    // bm.right = string -> num
    // X -> Y, where X is a counter and Y is a label
    typedef bimap<size_t, string> bm_t;
    bm_t label_id_map;

    typedef vector<vertex_descriptor> ver_to_id_t;
    // 2-D array: first dimension represents level of wl, second dimension represents the vertex, cell stores the label id
    vector<ver_to_id_t> ver_to_id_1(depth, ver_to_id_t(num_vertices(g_1))), ver_to_id_2(depth, ver_to_id_t(num_vertices(g_2)));
    // used for current depth
    size_t level = 0;

    wl_initialize(g_1, ver_to_id_1, label_id_map);
    wl_initialize(g_2, ver_to_id_2, label_id_map);

    for (level = 1; level < depth; level++)
    {
        // compare label of previous level 
        for (const auto& v_1 : make_iterator_range(vertices(g_1)))
        {
            // container for neighbor label id
            std::multiset<size_t> neighbor_id;
            // id of vertex label in previous level
            size_t this_node_id = ver_to_id_1[level-1][v_1];
            // DEBUG
            cout << "level " << level-1 << " | [" << this_node_id << ", " << label_id_map.left.find(this_node_id)->second << "] : ";
            for (const auto& v_adj : make_iterator_range(adjacent_vertices(v_1, g_1)))
            {
                // label of current neighbor
                string this_neighbor_label = g_1[v_adj].label;
                // id of label of current neighbor
                auto find = label_id_map.right.find(this_neighbor_label);
                // label must be present
                assert(find != label_id_map.right.end());
                size_t this_neighbor_id = find->second;
                // add the id to container
                neighbor_id.insert(this_neighbor_id);
                cout << "[" << this_neighbor_id << ", " << this_neighbor_label << "] ";
            }
            cout<<endl;
            // transform multiset to vector where in first position is present current node label id
            vector<size_t> new_label(neighbor_id.begin(), neighbor_id.end());
            new_label.insert(new_label.begin(), this_node_id);

            helper_insert(label_id_map, new_label);
        }

        for (const auto& v_2 : make_iterator_range(vertices(g_2)))
        {

            for (const auto& v_adj : make_iterator_range(adjacent_vertices(v_2, g_2)))
            {

            }
        }
    }

    // print all labels
    cout << "number of distinct labels: " << label_id_map.size() << endl;
    cout << "label mapping [ID -> label]:" << endl;
    for (const auto& p : label_id_map.left)
        cout << p.first << " -> " << p.second << endl;
    cout << "label mapping [label -> ID]:" << endl;
    for (const auto& p : label_id_map.right)
        cout << p.first << " -> " << p.second << endl;

    for (const auto& id : ver_to_id_1[0])
        cout << "g_1 | id: " << id << " -> label: " << label_id_map.left.find(id)->second << endl;
    for (const auto& id : ver_to_id_2[0])
        cout << "g_2 | id: " << id << " -> label: " << label_id_map.left.find(id)->second << endl;


    // count common labels
/*
    vector<size_t> phi_g_1(label_id_map.size(), 0);
    vector<size_t> phi_g_2(label_id_map.size(), 0);

    for (level = 0; level < depth; level++)
    {
        cout << "depth: " << level << ", features in g_1" << endl;
        for (const auto& [this_v, this_hash] : ver_to_id_1[level])
        {
            cout << "vertex descriptor: " << this_v << " | label ID: " << this_hash << " <-> " << id_to_lab[this_hash] << endl;
            phi_g_1[hash_pos_map[this_hash]]++;
        }

        cout << "depth: " << level << ", features in g_2" << endl;
        for (const auto& [this_v, this_hash] : ver_to_id_2[level])
        {
            cout << "vertex descriptor: " << this_v << " | label ID: " << this_hash << " <-> " << id_to_lab[this_hash] << endl;
            //cout << "deb " << this_hash << " " << hash_pos_map[this_hash] << endl;
            phi_g_2[hash_pos_map[this_hash]]++;
        }
    }

    cout << "Phi g_1:" << endl;
    for (auto p : phi_g_1)
        cout << p << " ";
    cout << endl;

    cout << "Phi g_1:" << endl;
    for (auto p : phi_g_2)
        cout << p << " ";
    cout << endl;

    // result is inner product of the two feature vectors Phi
    kernel = std::inner_product(phi_g_1.begin(), phi_g_1.end(), phi_g_2.begin(), 0);
*/
    cout << "result = ";
    return kernel;
}

void load_graphs(graph_vector& G_, size_t n, const string& filename)
{
    std::ifstream read(filename);

    if (filename == "shock.txt" && (n < 0 || n>149))
        n = 1;
    if (filename == "ppi.txt" && (n < 0 || n>85))
        n = 1;

    int i = 0;
    // read until you reach the end of the file or the number of graphs wanted
    for (string line; getline(read, line) && i != n; ++i) {

        // inserting the line into a stream that helps to parse the content
        stringstream ss(line);

        size_t num_vertices, source, target;
        char comma;

        ss >> num_vertices;

        //cout << "num of vertices: " << num_vertices << endl; //DEBUG
        labeled_graph_type g_(num_vertices);
        //vertex_name_map_t name_map = get(vertex_name, g_);
        //edge_weight_map_t weight_map = get(edge_weight, g_);

        while (ss >> source >> comma >> target)
        {
            edge_descriptor ed;
            bool inserted;

            // in txt file nodes start from 1
            source--; target--;

            //cout << "(" << source << "," << target << ")" << " "; //DEBUG
            tie(ed, inserted) = add_edge(source, target, g_);

            if(inserted)
            {
                vertex_descriptor vs, vt;
                vs = boost::source(ed, g_);
                vt = boost::target(ed, g_);
                g_[vs].label = "v" + std::to_string(source);
                g_[vt].label = "v" + std::to_string(target);
                // DEBUG
                //cout << "s: " << g_[vs].label << endl; // DEBUG
                //cout << "t: " << g_[vt].label << endl; // DEBUG
                //cout << "w: " << g_[ed].weight << endl; // DEBUG
            }
        }

        g_[graph_bundle].m_graph_index = i;
        g_[graph_bundle].m_density = density(g_);
        G_.push_back(g_);
    }

    read.close();
}

void print_graph_vector_info(const graph_vector& G_)
{
    for (const auto& g_ : G_)
        print_graph_info(g_);
}

void print_graph_info(const labeled_graph_type& g)
{
    cout << "\n**************************************************************************************************\n";
    cout << "Graph ID = " << g[graph_bundle].m_graph_index << endl;
    cout << "|vertices(g)| = " << num_vertices(g) << ",\t";
    cout << "|edges(g)| = " << num_edges(g) << endl;

    print_vertex_labels(g);

    cout << endl;

    print_weights(g);

    double d = g[graph_bundle].m_density;

    cout.precision(3);

    if (d > 0.7)
        cout << "\n\nGraph is dense (D = " << fixed << d << ")";
    else
        cout << "\n\nGraph is sparse (D = " << fixed << d << ")";

    cout << "\n**************************************************************************************************\n";
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