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
#include <boost/container_hash/hash.hpp>
#include <utility>
#include "test_graphs.h"


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
    //const auto& g_1 = test_graph_7();
    //const auto& g_2 = test_graph_8();

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

void initialize_hash_table(hash_table_t& hash_table, const labeled_graph_type& g)
{
    // initialize first level of hash table
    int level = 0;
    // for each node, hash its label and insert it in corresponding cell in hash table
    for (const auto& v : make_iterator_range(vertices(g)))
    {
        hash_table[level][v] = hash_value(g[v].label);
        cout << "hash: " << hash_table[level][v] << " - label: " << g[v].label << endl;
    }

}

void next_level_hash_table(hash_table_t& hash_table, const labeled_graph_type& g, size_t current_level)
{
    auto previous_level = current_level - 1;
    cout << "previous level: " << previous_level << endl;
    // starting from each node examine its hashed label and its neighbors' hashed label
    for (const auto& v : make_iterator_range(vertices(g)))
    {
        size_t this_node_hash = hash_table[previous_level][v];
        cout << "this node hash: " << this_node_hash;
        multiset<size_t> neighbor_node_hashes;
        // add each adjacent hashed label to multiset
        cout << " | neighbor hash: ";
        for (const auto& v_adj : make_iterator_range(adjacent_vertices(v, g)))
        {
            auto neighbor_hash = hash_table[previous_level][v_adj];
            neighbor_node_hashes.insert(neighbor_hash);
            cout << neighbor_hash << " ";
        }
        cout << endl;
        // transform multiset in vector
        vector<size_t> new_node_label(neighbor_node_hashes.begin(), neighbor_node_hashes.end());
        // add in first position its previous hashed label
        new_node_label.insert(new_node_label.begin(), this_node_hash);
        // hash the new label and insert it in corresponding cell in hash table
        hash_table[current_level][v] = hash_range(new_node_label.begin(), new_node_label.end());

        // DEBUG
        cout << "new label: ";
        for (const auto& h : new_node_label)
            cout << h << " ";
        cout << "with hash: " << hash_table[current_level][v] << endl;
    }
}

size_t weisfeiler_lehman_kernel(const labeled_graph_type& g_1, const labeled_graph_type& g_2, size_t depth)
{
    size_t kernel = 0;

    // hash tables containing label hash for each vertex for each level of depth in W-L algorithm
    hash_table_t hash_table_g_1(depth, vector<size_t>(num_vertices(g_1)));
    hash_table_t hash_table_g_2(depth, vector<size_t>(num_vertices(g_2)));

    // initialize hash tables
    initialize_hash_table(hash_table_g_1, g_1);
    initialize_hash_table(hash_table_g_2, g_2);

    // compute the new hashed labels for each level for each graph
    for (size_t level = 1; level < depth; level++)
    {
        cout << "level: " << level << " graph: 1" << endl;
        next_level_hash_table(hash_table_g_1, g_1, level);
        cout << "level: " << level << " graph: 2" << endl;
        next_level_hash_table(hash_table_g_2, g_2, level);
    }

    // create a set containing every unique hashed label
    set<size_t> hashed_labels;

    // scan first graph hash table
    for(const auto& row : hash_table_g_1)
        hashed_labels.insert(row.begin(), row.end());
    // scan second graph hash table
    for(const auto& row : hash_table_g_2)
        hashed_labels.insert(row.begin(), row.end());

    // print all labels
    cout << "number of distinct hashed labels: " << hashed_labels.size() << endl;
/*
    cout << "list of hashed labels:" << endl;
    for (const auto& l : hashed_labels)
        cout << '\t' << l << endl;
*/
    // DEBUG
    cout << "hash table 1: " << endl;
    for(int level = 0; level < depth; level++)
    {
        cout << "\tlevel " << level << endl << '\t';
        for(const auto& h : hash_table_g_1[level])
            cout << h << " ";
        cout << endl;
    }
    cout << "hash table 2: " << endl;
    for(int level = 0; level < depth; level++)
    {
        cout << "\tlevel " << level << endl << '\t';
        for(const auto& h : hash_table_g_2[level])
            cout << h << " ";
        cout << endl;
    }



    // count common labels
    vector<size_t> phi_g_1(hashed_labels.size(), 0);
    vector<size_t> phi_g_2(hashed_labels.size(), 0);

    // for each label
    for (auto [position, label_itr] = std::tuple{0, hashed_labels.begin()}; label_itr != hashed_labels.end(); label_itr++, position++)
    {
        // count how many times label_itr's present in the hash tables
        for(const auto& row : hash_table_g_1)
                phi_g_1[position] = count(row.begin(), row.end(), *label_itr);

        for(const auto& row : hash_table_g_2)
                phi_g_2[position] = count(row.begin(), row.end(), *label_itr);
    }

    cout << "Phi g 1:" << endl;
    for (auto p : phi_g_1)
        cout << p << " ";
    cout << endl;

    cout << "Phi g 2:" << endl;
    for (auto p : phi_g_2)
        cout << p << " ";
    cout << endl;

    // result is inner product of the two feature vectors Phi
    kernel = std::inner_product(phi_g_1.begin(), phi_g_1.end(), phi_g_2.begin(), 0);

    // check if hash is order dependent
    vector<size_t> a{hash_value("v3"),hash_value("v2")};
    vector<size_t> b{hash_value("v2"),hash_value("v3")};

    auto ha = hash_range(a.begin(), a.end());
    auto hb = hash_range(b.begin(), b.end());
    cout << boolalpha << (ha == hb) << endl;

    cout << "result = ";
    return kernel;
}

void load_graphs(graph_vector& G_, size_t n, const string& filename)
{
    std::ifstream read(filename);

    if (filename == "data/shock.txt" && (n < 0 || n>149))
        n = 1;
    if (filename == "data/ppi.txt" && (n < 0 || n>85))
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
