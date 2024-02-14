#include <bits/stdc++.h>
using namespace std;
#define LL long long
#define LD long double
#define PII pair<int, int>
#define inf 0x3fffffff

struct Edge
{
    int to;
    double p;
    Edge(int a = 0, double b = 0) : to(a), p(b) {}
};
mt19937 rand_num(2023);
vector<vector<Edge>> e;

bool check(int x)
{
    uniform_int_distribution<int> dist(1, x);
    if (dist(rand_num) == 1)
        return true;
    return false;
}

void compute_seed_sample(vector<vector<int>> &e, vector<double> &seed_contri, vector<vector<int>> &seeds, int n)
{
    queue<int> Q;
    vector<int> dist;
    vector<vector<int>> e_new;
    vector<int> in_degree;
    vector<int> valid;
    vector<double> contribution;
    e_new.resize(n);
    for (int i = 0; i < n; i++)
        dist.push_back(0), in_degree.push_back(0),
            valid.push_back(0), contribution.push_back(0.0);
    for (auto seed_list : seeds)
        for (auto seed : seed_list)
            Q.push(seed),
                dist[seed] = 1;
    while (!Q.empty())
    {
        int x = Q.front();
        valid[x] = 1;
        Q.pop();
        for (auto y : e[x])
            if (dist[y] == 0)
            {
                e_new[y].push_back(x);
                in_degree[x]++;
                Q.push(y);
                dist[y] = dist[x] + 1;
            }
            else if (dist[y] == dist[x] + 1)
            {
                e_new[y].push_back(x);
                in_degree[x]++;
            }
    }
    for (int i = 0; i < n; i++)
        if (in_degree[i] == 0 && valid[i] == 1)
            Q.push(i);
    while (!Q.empty())
    {
        int x = Q.front();
        Q.pop();
        contribution[x] += 1.0;
        double z = (contribution[x]) / ((double)e_new[x].size());
        for (auto y : e_new[x])
        {
            contribution[y] += z;
            in_degree[y]--;
            if (in_degree[y] == 0)
                Q.push(y);
        }
    }

    for (auto seed_list : seeds)
        for (auto seed : seed_list)
            seed_contri[seed] += contribution[seed];
}

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        printf("Usage: ./mic data_file seed_file budget 𝜀 l\n");
        return 0;
    }
    int n, m;
    char input[100], output[100], seed_input[100];
    int budget = atoi(argv[3]);
    sprintf(input, "../data/%s", argv[1]);
    sprintf(seed_input, "%s", argv[2]);
    sprintf(output, "../result/mic-b=%d-%s", budget, argv[1]);
    ifstream in(input);
    ifstream seedin(seed_input);
    ofstream out(output);
    double eps = atof(argv[4]) / 2;
    double l = atof(argv[5]);

    in >> n >> m;
    e.resize(n);
    for (int i = 1; i <= m; i++)
    {
        int x, y;
        double z;
        in >> x >> y >> z;
        e[x].push_back(Edge(y, z));
    }

    vector<int> seed_bool;
    seed_bool.resize(n);
    for (int i = 0; i < n; i++)
        seed_bool[i] = 0;
    int num_campaign, num_seeds = 0;
    vector<int> num_seed;
    vector<vector<int>> seeds;
    seedin >> num_campaign;
    seeds.resize(num_campaign);
    for (int i = 0; i < num_campaign; i++)
    {
        int num, x;
        seedin >> num;
        num_seeds += num;
        num_seed.push_back(num);
        while (num--)
        {
            seedin >> x;
            seed_bool[x] = i + 1;
            seeds[i].push_back(x);
        }
    }

    double beginTime = clock();
    vector<double> seed_contri;
    for (int i = 0; i < n; i++)
        seed_contri.push_back(0.0);
    int theta = int((2.0 + eps) * ((double)n) * (l * log((double)n)) / ((double)budget) / eps / eps + 1);
    out << "theta: " << theta << endl;
    cout << "theta: " << theta << endl;

    int T = max(0, theta);
    double Times = T;
    while (T--)
    {
        vector<vector<int>> e_sample;
        e_sample.resize(n);
        for (int i = 0; i < n; i++)
            for (auto j : e[i])
            {
                int x = i, y = j.to;
                double z = j.p;
                if (check((int)(1.0 / z)))
                    e_sample[x].push_back(y);
            }
        compute_seed_sample(e_sample, seed_contri, seeds, n);
    }
    for (int i = 0; i < n; i++)
        seed_contri[i] /= Times;
    for (auto seed : seeds[0])
        seed_contri[seed] = 0;

    out << "counters:" << endl;
    cout << "counters:" << endl;
    for (int i = 0; i < budget; i++)
    {
        double Max = -1;
        int seed = 0;
        for (int j = 0; j < n; j++)
            if (seed_contri[j] > Max)
                Max = seed_contri[j], seed = j;
        out << i << "\t" << seed << "\t" << seed_contri[seed] << endl;
        printf("%d\t%d\t%.10f\n", i, seed, seed_contri[seed]);
        seed_contri[seed] = 0;
    }

    double endTime = clock();
    out << "time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;
    cout << "time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;
    return 0;
}