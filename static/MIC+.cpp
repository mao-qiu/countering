#include <bits/stdc++.h>
using namespace std;
#define LL long long
#define LD long double
#define PII pair<int, int>
#define inf 0x3fffffff

struct Edge
{
    int to, tag;
    double p;
    Edge(int a = 0, double b = 0, int c = 0) : to(a), p(b), tag(c) {}
};
mt19937 rand_num(2023);
vector<vector<Edge>> e, re;
vector<int> vis;

bool check(int x)
{
    uniform_int_distribution<int> dist(1, x);
    if (dist(rand_num) == 1)
        return true;
    return false;
}

void compute_seed_sample(int x, vector<double> &seed_contri, vector<int> &seed_bool, int n, int tag)
{
    if (seed_bool[x] != 0)
    {
        seed_contri[x] += 1.0;
        return;
    }
    vector<int> d;
    d.resize(n);

    vector<vector<int>> e_sample;
    e_sample.resize(n);
    queue<int> Q;
    d[x] = 1;
    Q.push(x);
    vis[x] = tag;
    vector<int> seed;
    int seed_dist = inf;
    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        if (seed_dist != inf && d[u] == seed_dist)
            continue;
        for (auto j : re[u])
        {
            int v = j.to;
            double z = j.p;
            if (check((int)(1.0 / z)))
            {
                e_sample[v].push_back(u);
                if (vis[v] != tag)
                {
                    vis[v] = tag;
                    d[v] = d[u] + 1;
                    Q.push(v);
                    if (seed_bool[v] != 0)
                    {
                        seed.push_back(v);
                        seed_dist = min(seed_dist, d[v]);
                    }
                }
            }
        }
    }
    if (seed.size() == 0)
        return;
    // cout<<"x: "<<x<<endl;
    // for (auto u:seed)
    //     cout<<u<<" ";
    // cout<<endl;
    vector<int> scan, atime;
    vector<double> contribution;
    atime.resize(n);
    contribution.resize(n);
    for (auto u : seed)
        Q.push(u), atime[u] = 1, vis[u] = -tag;
    vector<vector<int>> e_influence;
    e_influence.resize(n);
    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        scan.push_back(u);
        contribution[u] = 0.0;
        for (auto v : e_sample[u])
        {
            if (vis[v] != -tag)
            {
                vis[v] = -tag;
                Q.push(v);
                atime[v] = atime[u] + 1;
            }
            if (atime[v] == atime[u] + 1)
                e_influence[v].push_back(u);
        }
    }
    contribution[x] = 1.0;
    int len = scan.size();
    for (int i = len - 1; i >= 0; i--)
    {
        int u = scan[i];

        // if (contribution[u] > 0)
        //     cout << u << " " << contribution[u] << " "<<x<<endl;
        double p = e_influence[u].size();
        for (auto v : e_influence[u])
            contribution[v] += contribution[u] / p;
    }

    for (auto u : seed)
        seed_contri[u] += contribution[u];
}

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        printf("Usage: ./micp data_file seed_file budget 𝜀 l\n");
        return 0;
    }
    int n, m;
    char input[100], output[100], seed_input[100];
    int budget = atoi(argv[3]);
    sprintf(input, "../data/%s", argv[1]);
    sprintf(seed_input, "%s", argv[2]);
    sprintf(output, "../result/micp-b=%d-%s", budget, argv[1]);
    ifstream in(input);
    ifstream seedin(seed_input);
    ofstream out(output);
    double eps = atof(argv[4]);
    double l = atof(argv[5]);
    in >> n >> m;
    e.resize(n);
    re.resize(n);
    for (int i = 1; i <= m; i++)
    {
        int x, y;
        double z;
        in >> x >> y >> z;
        e[x].push_back(Edge(y, z, 0));
        re[y].push_back(Edge(x, z, 0));
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
    int r = (int)((eps + 2.0) * ((double)n) * l * log((double)n) / ((double)num_seeds - seeds[0].size()) / eps / eps + 1);
    out << "r: " << r << endl;
    cout << "r: " << r << endl;
    vector<double> seed_contri;
    seed_contri.resize(n);
    vis.resize(n);
    for (int i = 0; i < n; i++)
        seed_contri[i] = 0.0, vis[i] = 0;
    int RRnum = 0;
    int T = r;
    double Times = T;
    while (T--)
    {
        uniform_int_distribution<int> dist(0, n - 1);
        int x = dist(rand_num);
        RRnum++;
        compute_seed_sample(x, seed_contri, seed_bool, n, RRnum);
    }
    for (auto seed : seeds[0])
        seed_contri[seed] = 0;
    double KPT = 0;
    for (int i = 0; i < n; i++)
        KPT += seed_contri[i] / Times;
    KPT = KPT * ((double)n) * ((double)budget) / ((double)num_seeds);
    KPT /= (1 + eps);

    int theta = int((eps + 4.0) * 2.0 * ((double)n) * (l * log((double)n)) / KPT / eps / eps + 1);
    out << "KPT: " << KPT << "\ttheta: " << theta << endl;
    cout << "KPT: " << KPT << "\ttheta: " << theta << endl;

    T = max(0, theta - r);
    Times = max(Times, (double)theta);
    while (T--)
    {
        uniform_int_distribution<int> dist(0, n - 1);
        int x = dist(rand_num);
        RRnum++;
        compute_seed_sample(x, seed_contri, seed_bool, n, RRnum);
    }
    for (int i = 0; i < n; i++)
    {
        seed_contri[i] /= Times;
        seed_contri[i] *= (double)n;
    }
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
