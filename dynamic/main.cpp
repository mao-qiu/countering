#include "dynamic.h"
using namespace std;
mt19937 rand_edge(2023);

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        printf("Usage: ./main data_file seed_file budget 𝜀 l\n");
        return 0;
    }
    int n, m;
    char input[100], output[100], seed_input[100];
    int budget = atoi(argv[3]);
    sprintf(input, "../data/%s", argv[1]);
    sprintf(seed_input, "%s", argv[2]);
    sprintf(output, "../result/dynamic-b=%d-%s", budget, argv[1]);
    ifstream in(input);
    ifstream seedin(seed_input);
    ofstream out(output);
    double eps = atof(argv[4]);
    double l = atof(argv[5]);

    in >> n >> m;
    vector<vector<Edge>> e;
    vector<pair<PII, double>> edges;
    e.resize(n);
    for (int i = 1; i <= m; i++)
    {
        int x, y;
        double z;
        in >> x >> y >> z;
        Edge tmp;
        tmp.to = y;
        tmp.p = z;
        e[x].push_back(tmp);
        edges.push_back(make_pair(PII(x, y), z));
    }
    int num_campaign, num_seeds = 0;
    vector<vector<int>> seeds;
    seedin >> num_campaign;
    seeds.resize(num_campaign);
    for (int i = 0; i < num_campaign; i++)
    {
        int num, x;
        seedin >> num;
        num_seeds += num;
        while (num--)
        {
            seedin >> x;
            seeds[i].push_back(x);
        }
    }

    double beginTime = clock();
    Dynamic dynamic;
    dynamic.init(n, m, budget, e, seeds, eps, l);
    vector<Con> result = dynamic.top_con();
    for (auto x : result)
        cout << x.x << "\t" << x.con << endl,
            out << x.x << "\t" << x.con << endl;
    double endTime = clock();
    out << "build index time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;
    cout << "build index time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;

    double maxTime = 0, minTime = 10000000000;
    beginTime = clock();
    int change_edge = 3;
    vector<int> vis_edge;
    vis_edge.resize(m);
    for (int i = 0; i < change_edge; i++)
    {
        double singleBegin = clock();
        uniform_int_distribution<int> dist(0, m - 1);
        int x = dist(rand_edge);
        while (vis_edge[x] == 1)
            x = dist(rand_edge);
        vis_edge[x] = 1;
        int u = edges[x].first.first, v = edges[x].first.second;
        printf("remove edge (%d -> %d)\n", u, v);
        dynamic.erase_edge(u, v);
        double singleEnd = clock();
        maxTime = max(maxTime, singleEnd - singleBegin);
        minTime = min(minTime, singleEnd - singleBegin);
    }
    endTime = clock();
    out << "remove time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;
    cout << "remove time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;
    out << "remove average time : " << (endTime - beginTime) / CLOCKS_PER_SEC / change_edge << "s." << endl;
    cout << "remove average time : " << (endTime - beginTime) / CLOCKS_PER_SEC / change_edge << "s." << endl;
    out << "remove max time : " << maxTime / CLOCKS_PER_SEC << "s." << endl;
    cout << "remove max time : " << maxTime / CLOCKS_PER_SEC << "s." << endl;
    out << "remove min time : " << minTime / CLOCKS_PER_SEC << "s." << endl;
    cout << "remove min time : " << minTime / CLOCKS_PER_SEC << "s." << endl;

    maxTime = 0, minTime = 10000000000;
    beginTime = clock();
    for (int i = 0; i < m; i++)
        if (vis_edge[i] == 1)
        {
            double singleBegin = clock();
            int u = edges[i].first.first, v = edges[i].first.second;
            double p = edges[i].second;
            printf("insert edge (%d -> %d), p = %.5f\n", u, v, p);
            dynamic.insert_edge(u, v, p);
            double singleEnd = clock();
            maxTime = max(maxTime, singleEnd - singleBegin);
            minTime = min(minTime, singleEnd - singleBegin);
        }
    endTime = clock();
    out << "insert time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;
    cout << "insert time : " << (endTime - beginTime) / CLOCKS_PER_SEC << "s." << endl;
    out << "insert average time : " << (endTime - beginTime) / CLOCKS_PER_SEC / change_edge << "s." << endl;
    cout << "insert average time : " << (endTime - beginTime) / CLOCKS_PER_SEC / change_edge << "s." << endl;
    out << "insert max time : " << maxTime / CLOCKS_PER_SEC << "s." << endl;
    cout << "insert max time : " << maxTime / CLOCKS_PER_SEC << "s." << endl;
    out << "insert min time : " << minTime / CLOCKS_PER_SEC << "s." << endl;
    cout << "insert min time : " << minTime / CLOCKS_PER_SEC << "s." << endl;

    return 0;
}