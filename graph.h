//
// Created by loryt on 8/23/2020.
//

#ifndef BGL_GRAPH_H
#define BGL_GRAPH_H

#include <fstream>

using namespace boost;
using namespace std;
/* STRUCTS FOR KERNEL FEATURES */
typedef vector<vector<size_t>> distance_matrix_t;
typedef vector<vector<size_t>> hash_table_t;
typedef std::set<std::tuple<string, string, size_t>> tuple_set_t;
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
void print_graph_info(const labeled_graph_type& g);
inline double density(const labeled_graph_type& g);
void load_graphs(graph_vector& G_, size_t n = 2, const string& filename = "data/shock.txt");
void print_graph_vector_info(const graph_vector& G_);

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

// calc graph density. It's used to decide which algo to use later.
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

// func to print all vertex names
void print_vertex_labels(const labeled_graph_type& g)
{
    cout << "vertices(g) = { ";
    for (auto v : make_iterator_range(vertices(g)))
        std::cout << g[v].label << " ";
    cout << "}" << endl;
}

// func to print all edge weights
void print_weights(const labeled_graph_type& g)
{
    for (const auto& e : make_iterator_range(edges(g)))
        cout << g[e.m_source].label << " -> " << g[e.m_target].label << " : " << g[e].weight << endl;
}

// method loads n graphs from file "shock.txt" into a graph vector
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

// method prints all graphs info contained in vector
void print_graph_vector_info(const graph_vector& G_)
{
    for (const auto& g_ : G_)
        print_graph_info(g_);
}

// visualize graph vertices, edges and density.
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

#endif //BGL_GRAPH_H
