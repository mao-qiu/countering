#include <bits/stdc++.h>
#define LL long long
using namespace std;
#define PII pair<int, int>

struct Con
{
    int x;
    double con;
    bool operator<(const Con &x) const
    {
        return con > x.con;
    }
};

class Edge
{
public:
    int to;
    double p;
};

struct hashFunction
{
    size_t operator()(const pair<int, int> &x) const
    {
        return (x.first << 17) ^ x.second;
    }
};

class Index
{
public:
    int x, dist2seed;
    unordered_map<int, unordered_set<int>> rh, h;
    unordered_map<int, double> con;
    unordered_map<int, int> dist;
    unordered_set<PII, hashFunction> live_edge, dead_edge;
};

class Dynamic
{
private:
    int n, m, budget, num_seeds;
    int r, num_sample;
    double eps, l, kpt;
    unordered_set<int> seed_list;
    unordered_map<int, unordered_set<int>> seeds;
    unordered_map<int, unordered_map<int, double>> e, re;
    unordered_map<int, int> seed_bool;
    unordered_map<int, double> sum_con;
    unordered_set<int> valid_vertex_id;
    bool check(int x);
    void generate_sample(int index_id);
    void reconstruct(int index_id);
    void update(int u, Index &index, int index_id);
    void change_parameter();

public:
    vector<Index> index;
    void init(int n, int m, int budget, vector<vector<Edge>> &e, vector<vector<int>> &seeds, double eps, double l);
    void insert_vertex(int v);
    void erase_vertex(int v);
    void insert_edge(int u, int v, double p);
    void erase_edge(int u, int v);
    void insert_seed(int v, int campaign_id);
    void erase_seed(int v);
    void change_budget(int b);
    void change_probability(int u, int v, double p);

    vector<Con> top_con();
};