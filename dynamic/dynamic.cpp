#include "dynamic.h"
using namespace std;
mt19937 rand_num(2023);
#define inf 0x3fffffff

bool Dynamic::check(int x)
{
    uniform_int_distribution<int> dist(1, x);
    if (dist(rand_num) == 1)
        return true;
    return false;
}

void Dynamic::generate_sample(int index_id)
{
    uniform_int_distribution<int> dist(0, n - 1);
    int x = dist(rand_num);
    Index sample;
    sample.x = x;
    if (seed_bool.find(x) != seed_bool.end())
    {
        sample.con[x] = 1;
        sum_con[x] += 1;
        if (index_id == -1)
            index.push_back(sample), num_sample++;
        else
            index[index_id] = sample;
        return;
    }
    vector<int> d;
    d.resize(n);

    vector<vector<int>> e_sample;
    vector<int> vis;
    e_sample.resize(n);
    vis.resize(n);
    queue<int> Q;
    d[x] = 1;
    Q.push(x);
    vis[x] = 1;
    int seed_dist = inf;
    vector<int> Seed;
    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        if (seed_dist != inf && d[u] == seed_dist)
            continue;
        for (auto j : re[u])
        {
            int v = j.first;
            double z = j.second;
            if (check((int)(1.0 / z)))
            {
                sample.live_edge.insert(PII(u, v));
                e_sample[v].push_back(u);
                if (vis[v] != 1)
                {
                    vis[v] = 1;
                    d[v] = d[u] + 1;
                    Q.push(v);
                    if (seed_bool.find(v) != seed_bool.end())
                    {
                        Seed.push_back(v);
                        seed_dist = min(seed_dist, d[v]);
                    }
                }
            }
            else
                sample.dead_edge.insert(PII(u, v));
        }
    }
    if (Seed.size() == 0)
    {
        if (index_id == -1)
            index.push_back(sample), num_sample++;
        else
            index[index_id] = sample;
        return;
    }
    sample.dist2seed = seed_dist;
    vector<int> scan;
    vector<double> contribution;
    contribution.resize(n);
    for (auto u : Seed)
        Q.push(u), sample.dist[u] = seed_dist + 1, vis[u] = -1;
    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        scan.push_back(u);
        contribution[u] = 0.0;
        for (auto v : e_sample[u])
        {
            if (vis[v] != -1)
            {
                vis[v] = -1;
                Q.push(v);
                sample.dist[v] = sample.dist[u] - 1;
            }
            if (sample.dist[v] == sample.dist[u] - 1)
                sample.h[v].insert(u),
                    sample.rh[u].insert(v);
        }
    }

    contribution[x] = 1.0;
    int len = scan.size();
    for (int i = len - 1; i >= 0; i--)
    {
        int u = scan[i];
        sample.con[u] = contribution[u];
        sum_con[u] += contribution[u];
        double p = sample.h[u].size();
        for (auto v : sample.h[u])
            contribution[v] += contribution[u] / p;
    }
    if (index_id == -1)
        index.push_back(sample), num_sample++;
    else
        index[index_id] = sample;
}

void Dynamic::init(int n, int m, int budget, vector<vector<Edge>> &e, vector<vector<int>> &seeds, double eps, double l)
{
    this->n = n;
    this->m = m;
    this->budget = budget;
    for (int i = 0; i < (int)seeds.size(); i++)
    {
        for (auto x : seeds[i])
            this->seeds[i].insert(x);
    }
    this->num_sample = 0;
    this->l = l;
    this->eps = eps;
    this->e.clear();
    this->re.clear();
    this->seed_bool.clear();
    this->sum_con.clear();
    for (int x = 0; x < n; x++)
    {
        this->valid_vertex_id.insert(x);
        for (auto y : e[x])
        {
            this->e[x][y.to] = y.p;
            this->re[y.to][x] = y.p;
        }
    }
    this->seed_list.clear();
    int campaign_num = seeds.size();
    for (auto campaign_id = 0; campaign_id < campaign_num; campaign_id++)
        for (auto x : seeds[campaign_id])
            seed_bool[x] = campaign_id + 1,
            this->seed_list.insert(x);
    this->num_seeds = this->seed_list.size();
    int r = (int)((eps + 2.0) * ((double)n) * l * log((double)n) / ((double)num_seeds - seeds[0].size()) / eps / eps + 1);
    this->r = r;
    cout << "r: " << r << endl;
    for (int i = 0; i < r; i++)
        generate_sample(-1);
    vector<double> con;
    con.resize(n);
    for (auto s : seed_list)
        con[s] = sum_con[s];
    for (auto s : seeds[0])
        con[s] = 0;
    double KPT = 0;
    for (auto s : seed_list)
        KPT += con[s] / ((double)r);
    KPT = KPT * ((double)n) * ((double)budget) / ((double)num_seeds);
    KPT /= (1 + eps);
    int theta = int((4.0 + eps) * 2.0 * ((double)n) * (l * log((double)n)) / KPT / eps / eps + 1);
    this->kpt = KPT;
    cout << "KPT: " << KPT << endl;
    cout << "theta: " << theta << endl;
    for (int i = r; i < max(theta, r); i++)
        generate_sample(-1);
}

vector<Con> Dynamic::top_con()
{
    vector<double> con;
    con.resize(n);
    for (auto s : seed_list)
        con[s] = sum_con[s];
    for (auto s : seeds[0])
        con[s] = 0;
    vector<Con> S;
    for (auto s : seed_list)
    {
        Con tmp;
        tmp.x = s;
        tmp.con = con[s] * ((double)n) / ((double)num_sample);
        S.push_back(tmp);
    }
    sort(S.begin(), S.end());
    vector<Con> result;
    for (int i = 0; i < budget; i++)
        result.push_back(S[i]);
    return result;
}

void Dynamic::update(int u, Index &ind, int index_id)
{
    queue<int> Q;
    Q.push(u);
    vector<int> vis;
    vis.resize(n);
    vis[u] = 1;
    int new_dist2seed = inf;
    vector<int> seed;
    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        if (seed_bool.find(u) != seed_bool.end())
        {
            new_dist2seed = min(new_dist2seed, ind.dist[u]);
            seed.push_back(u);
        }
        if (ind.dist[u] >= new_dist2seed || ind.dist[u] >= ind.dist2seed)
            continue;
        for (auto j : re[u])
        {
            int v = j.first;
            double p = j.second;
            if (ind.dead_edge.find(PII(u, v)) != ind.dead_edge.end())
                continue;
            if (ind.live_edge.find(PII(u, v)) == ind.live_edge.end())
            {
                if (!check(1.0 / p))
                {
                    ind.dead_edge.insert(PII(u, v));
                    continue;
                }
                else
                    ind.live_edge.insert(PII(u, v));
            }
            if (ind.dist[u] + 1 <= ind.dist[v] || ind.rh[v].size() == 0)
            {
                if (vis[v] != 1)
                    vis[v] = 1, Q.push(v);
                ind.dist[v] = ind.dist[u] + 1;
            }
        }
    }
    if (new_dist2seed > ind.dist2seed)
        return;
    for (auto s : seed_list)
        sum_con[s] -= index[index_id].con[s];
    reconstruct(index_id);
}

void Dynamic::reconstruct(int index_id)
{
    int x = index[index_id].x;
    Index sample;
    sample.x = x;
    if (seed_bool.find(x) != seed_bool.end())
    {
        sample.con[x] = 1;
        sum_con[x] += 1;
        index[index_id] = sample;
        return;
    }
    vector<int> d;
    d.resize(n);

    vector<vector<int>> e_sample;
    vector<int> vis;
    e_sample.resize(n);
    vis.resize(n);
    queue<int> Q;
    d[x] = 1;
    Q.push(x);
    vis[x] = 1;
    int seed_dist = inf;
    vector<int> Seed;
    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        if (seed_dist != inf && d[u] == seed_dist)
            continue;
        for (auto j : re[u])
        {
            int v = j.first;
            double z = j.second;
            if (index[index_id].dead_edge.find(PII(u, v)) != index[index_id].dead_edge.end())
            {
                sample.dead_edge.insert(PII(u, v));
                continue;
            }
            if (index[index_id].live_edge.find(PII(u, v)) == index[index_id].live_edge.end())
            {
                if (check((int)(1.0 / z)))
                    sample.live_edge.insert(PII(u, v));
                else
                    sample.dead_edge.insert(PII(u, v));
            }
            else
                sample.live_edge.insert(PII(u, v));
            if (sample.live_edge.find(PII(u, v)) != sample.live_edge.end())
            {
                e_sample[v].push_back(u);
                if (vis[v] != 1)
                {
                    vis[v] = 1;
                    d[v] = d[u] + 1;
                    Q.push(v);
                    if (seed_bool.find(v) != seed_bool.end())
                    {
                        Seed.push_back(v);
                        seed_dist = min(seed_dist, d[v]);
                    }
                }
            }
        }
    }
    if (Seed.size() == 0)
    {
        index[index_id] = sample;
        return;
    }
    sample.dist2seed = seed_dist;
    vector<int> scan;
    vector<double> contribution;
    contribution.resize(n);
    for (auto u : Seed)
        Q.push(u), sample.dist[u] = seed_dist + 1, vis[u] = -1;
    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        scan.push_back(u);
        contribution[u] = 0.0;
        for (auto v : e_sample[u])
        {
            if (vis[v] != -1)
            {
                vis[v] = -1;
                Q.push(v);
                sample.dist[v] = sample.dist[u] - 1;
            }
            if (sample.dist[v] == sample.dist[u] - 1)
                sample.h[v].insert(u),
                    sample.rh[u].insert(v);
        }
    }

    contribution[x] = 1.0;
    int len = scan.size();
    for (int i = len - 1; i >= 0; i--)
    {
        int u = scan[i];
        sample.con[u] = contribution[u];
        sum_con[u] += contribution[u];
        double p = sample.h[u].size();
        for (auto v : sample.h[u])
            contribution[v] += contribution[u] / p;
    }
    index[index_id] = sample;
}

void Dynamic::insert_edge(int u, int v, double p)
{
    if (e[u].find(v) != e[u].end())
    {
        printf("Fail. (%d,%d) exists.\n", u, v);
        return;
    }
    m++;
    e[u][v] = p;
    re[v][u] = p;
    for (int i = 0; i < num_sample; i++)
    {
        if (!check((int)(1.0 / p)))
        {
            index[i].dead_edge.insert(PII(v, u));
            continue;
        }
        index[i].live_edge.insert(PII(v, u));
        if (index[i].con.find(v) == index[i].con.end())
        {
            if (index[i].con.size() == 0)
            {
                for (auto s : seed_list)
                    sum_con[s] -= index[i].con[s];
                reconstruct(i);
            }
            continue;
        }
        if (index[i].con.find(u) == index[i].con.end())
        {
            if (index[i].dist[v] == index[i].dist2seed)
                continue;
            index[i].dist[u] = index[i].dist[v] + 1;
            update(u, index[i], i);
        }
        if (index[i].dist[u] < index[i].dist[v] + 1)
            continue;
        if (index[i].dist[u] == index[i].dist[v] + 1)
        {
            if (index[i].h[v].size() > 0)
            {
                queue<pair<int, pair<double, double>>> Q;
                double con_old, con_new;
                con_new = index[i].con[v] / ((double)index[i].h[v].size() + 1);
                con_old = index[i].con[v] / ((double)index[i].h[v].size());
                Q.push(make_pair(v, make_pair(con_old, con_new)));
                while (!Q.empty())
                {
                    int x = Q.front().first;
                    double con_old = Q.front().second.first;
                    double con_new = Q.front().second.second;
                    Q.pop();
                    for (auto y : index[i].h[x])
                    {
                        index[i].con[y] = index[i].con[y] - con_old + con_new,
                        sum_con[y] = sum_con[y] - con_old + con_new;
                        if (index[i].h[y].size() == 0)
                            continue;
                        double p = index[i].h[y].size();
                        Q.push(make_pair(y, make_pair(con_old / p, con_new / p)));
                    }
                }
            }
            index[i].h[v].insert(u);
            index[i].rh[u].insert(v);
            for (auto s : seed_list)
                sum_con[s] -= index[i].con[s];
            reconstruct(i);
        }
        else
        {
            index[i].dist[u] = index[i].dist[v] + 1;
            update(u, index[i], i);
        }
    }
    vector<double> con;
    con.resize(n);
    for (auto s : seed_list)
        con[s] = sum_con[s];
    for (auto s : seeds[0])
        con[s] = 0;
    double KPT = 0;
    for (auto s : seed_list)
        KPT += con[s] / ((double)num_sample);
    KPT = KPT * ((double)n) * ((double)budget) / ((double)num_seeds);
    KPT /= (1 + eps);
    this->kpt = KPT;
    cout << "kpt:" << KPT << endl;
    int theta = int((4.0 + eps) * 2.0 * ((double)n) * (l * log((double)n)) / KPT / eps / eps + 1);
    int more_sample = theta - num_sample;
    for (int i = 0; i < more_sample; i++)
    {
        if (num_sample < (int)index.size())
            generate_sample(num_sample), num_sample++;
        else
            generate_sample(-1);
    }
}

void Dynamic::erase_edge(int u, int v)
{
    if (e[u].find(v) == e[u].end())
    {
        printf("Fail. (%d,%d) not exists.\n", u, v);
        return;
    }
    for (int i = 0; i < num_sample; i++)
    {
        if (index[i].con.find(v) == index[i].con.end())
            continue;
        if (index[i].con.find(u) == index[i].con.end())
            continue;
        if (index[i].dead_edge.find(PII(v, u)) != index[i].dead_edge.end())
        {
            index[i].dead_edge.erase(PII(v, u));
            continue;
        }
        if (index[i].live_edge.find(PII(v, u)) == index[i].live_edge.end())
        {
            if (check((int)(1.0 / e[u][v])))
                index[i].live_edge.insert(PII(v, u));
        }
        if (index[i].live_edge.find(PII(v, u)) == index[i].live_edge.end())
            continue;
        index[i].live_edge.erase(PII(v, u));
        if (index[i].dist[v] + 1 == index[i].dist[u])
        {
            index[i].h[v].erase(u);
            index[i].rh[u].erase(v);
            int num = index[i].rh[u].size();
            if (num > 0)
            {
                if (index[i].h[v].size() > 0)
                {
                    queue<pair<int, pair<double, double>>> Q;
                    double con_old, con_new;
                    con_new = 0;
                    con_old = index[i].con[v] / ((double)index[i].h[v].size() + 1);
                    Q.push(make_pair(u, make_pair(con_old, con_new)));
                    while (!Q.empty())
                    {
                        int x = Q.front().first;
                        double con_old = Q.front().second.first;
                        double con_new = Q.front().second.second;
                        Q.pop();
                        index[i].con[x] = index[i].con[x] - con_old + con_new,
                        sum_con[x] = sum_con[x] - con_old + con_new;
                        if (index[i].h[x].size() == 0)
                            continue;
                        double p = index[i].h[x].size();
                        con_old /= p;
                        con_new /= p;
                        for (auto y : index[i].h[x])
                        {
                            Q.push(make_pair(y, make_pair(con_old, con_new)));
                        }
                    }

                    con_old = index[i].con[v] / ((double)index[i].h[v].size() + 1);
                    con_new = index[i].con[v] / ((double)index[i].h[v].size());
                    Q.push(make_pair(v, make_pair(con_old, con_new)));
                    while (!Q.empty())
                    {
                        int x = Q.front().first;
                        double con_old = Q.front().second.first;
                        double con_new = Q.front().second.second;
                        Q.pop();
                        for (auto y : index[i].h[x])
                            if (y != u)
                            {
                                index[i].con[y] = index[i].con[y] - con_old + con_new,
                                sum_con[y] = sum_con[y] - con_old + con_new;
                                if (index[i].h[y].size() == 0)
                                    continue;
                                double p = index[i].h[y].size();
                                Q.push(make_pair(y, make_pair(con_old / p, con_new / p)));
                            }
                    }
                }
                else
                {
                    queue<int> Q;
                    Q.push(v);
                    while (!Q.empty())
                    {
                        int x = Q.front();
                        Q.pop();
                        for (auto y : index[i].rh[x])
                        {
                            if (index[i].h[y].size() == 1)
                                Q.push(y);
                            index[i].h[y].erase(x);
                        }
                        index[i].rh[x].clear();
                    }
                    vector<int> vis;
                    vis.resize(n);
                    Q.push(index[i].x);
                    vis[index[i].x] = 1;
                    for (auto tmp : index[i].con)
                    {
                        int s = tmp.first;
                        sum_con[s] -= tmp.second;
                    }
                    index[i].con.clear();
                    index[i].con[index[i].x] = 1;
                    while (!Q.empty())
                    {
                        int x = Q.front();
                        Q.pop();
                        sum_con[x] += index[i].con[x];
                        double p = index[i].h[x].size();
                        for (auto y : index[i].h[x])
                        {
                            index[i].con[y] += index[i].con[x] / p;
                            if (vis[y] != 1)
                                vis[y] = 1, Q.push(y);
                        }
                    }
                }
            }
            else
            {
                for (auto s : seed_list)
                    sum_con[s] -= index[i].con[s];
                reconstruct(i);
            }
        }
    }
    m--;
    e[u].erase(v);
    re[v].erase(u);
    vector<double> con;
    con.resize(n);
    for (auto s : seed_list)
        con[s] = sum_con[s];
    for (auto s : seeds[0])
        con[s] = 0;
    double KPT = 0;
    for (auto s : seed_list)
        KPT += con[s] / ((double)num_sample);
    KPT = KPT * ((double)n) * ((double)budget) / ((double)num_seeds);
    KPT /= (1 + eps);
    this->kpt = KPT;
    cout << "kpt:" << KPT << endl;
    int theta = int((4.0 + eps) * 2.0 * ((double)n) * (l * log((double)n)) / KPT / eps / eps + 1);
    int more_sample = theta - num_sample;
    for (int i = 0; i < more_sample; i++)
    {
        if (num_sample < (int)index.size())
            generate_sample(num_sample), num_sample++;
        else
            generate_sample(-1);
    }
}

void Dynamic::change_probability(int u, int v, double p)
{
    if (e[u].find(v) == e[u].end())
    {
        printf("Fail. (%d,%d) not exists.\n", u, v);
        return;
    }
    erase_edge(u, v);
    insert_edge(u, v, p);
}

void Dynamic::change_parameter()
{
    int newr = (int)((eps + 2.0) * ((double)n) * l * log((double)n) / ((double)num_seeds - seeds[0].size()) / eps / eps + 1);
    int more_sample = newr - num_sample;
    for (int i = 0; i < more_sample; i++)
    {
        if (num_sample < (int)index.size())
            generate_sample(num_sample), num_sample++;
        else
            generate_sample(-1);
    }
    this->r = newr;
    vector<double> con;
    con.resize(n);
    for (auto s : seed_list)
        con[s] = sum_con[s];
    for (auto s : seeds[0])
        con[s] = 0;
    double KPT = 0;
    for (auto s : seed_list)
        KPT += con[s] / ((double)r);
    KPT = KPT * ((double)n) * ((double)budget) / ((double)num_seeds);
    KPT /= (1 + eps);
    this->kpt = KPT;
    int theta = int((4.0 + eps) * 2.0 * ((double)n) * (l * log((double)n)) / KPT / eps / eps + 1);
    more_sample = theta - num_sample;
    for (int i = 0; i < more_sample; i++)
    {
        if (num_sample < (int)index.size())
            generate_sample(num_sample), num_sample++;
        else
            generate_sample(-1);
    }
}

void Dynamic::insert_vertex(int v)
{
    if (sum_con.find(v) != sum_con.end())
    {
        printf("Fail. Vertex %d exists.\n", v);
        return;
    }
    n++;
    sum_con[v] = 0;
    e[v].clear();
    valid_vertex_id.insert(v);
    if (seed_bool[v] > 0)
        insert_seed(v, seed_bool[v]);
    change_parameter();
    for (int i = 0; i < num_sample; i++)
    {
        if (check(n))
        {
            index[i].x = v;
            index[i].h.clear();
            index[i].rh.clear();
            index[i].con.clear();
            index[i].dist.clear();
            index[i].live_edge.clear();
            index[i].dead_edge.clear();
            for (auto s : seed_list)
                sum_con[s] -= index[i].con[s];
            reconstruct(i);
        }
    }
}

void Dynamic::erase_vertex(int v)
{
    if (sum_con.find(v) == sum_con.end())
    {
        printf("Fail. Vertex %d not exists.\n", v);
        return;
    }
    n--;
    if (seed_bool[v] > 0)
        erase_seed(v);
    for (auto y : e[v])
    {
        int u = y.first;
        erase_edge(v, u);
    }
    for (auto y : re[v])
    {
        int u = y.first;
        erase_edge(u, v);
    }
    valid_vertex_id.erase(v);
    sum_con.erase(v);
    change_parameter();
    for (int i = 0; i < num_sample; i++)
        if (index[i].x == v)
        {
            vector<int> x;
            x.resize(1);
            sample(valid_vertex_id.begin(), valid_vertex_id.end(), x.begin(), 1, rand_num);
            index[i].x = x[0];
            index[i].h.clear();
            index[i].rh.clear();
            index[i].con.clear();
            index[i].dist.clear();
            index[i].live_edge.clear();
            index[i].dead_edge.clear();
            for (auto s : seed_list)
                sum_con[s] -= index[i].con[s];
            reconstruct(i);
        }
}

void Dynamic::erase_seed(int v)
{
    if (seed_bool.find(v) == seed_bool.end())
    {
        printf("Fail. Seed %d not exists.\n", v);
        return;
    }
    if ((int)seed_bool.size() < budget)
    {
        printf("Fail. The budget is larger than the number of seeds.\n");
        return;
    }
    sum_con.erase(v);
    seeds[seed_bool[v] - 1].erase(v);
    seed_bool.erase(v);
    seed_list.erase(v);
    num_seeds--;
    change_parameter();
    for (int i = 0; i < num_sample; i++)
    {
        if (index[i].con.find(v) == index[i].con.end())
            continue;
        if (index[i].con[v] == 0)
            continue;
        for (auto s : seed_list)
            sum_con[s] -= index[i].con[s];
        reconstruct(i);
    }
}

void Dynamic::insert_seed(int v, int campaign_id)
{
    if (seed_bool.find(v) != seed_bool.end())
    {
        printf("Fail. Seed %d exists.\n", v);
        return;
    }
    if ((int)seed_bool.size() < budget)
    {
        printf("Fail. The budget is larger than the number of seeds.\n");
        return;
    }
    sum_con[v] = 0;
    seed_bool[v] = campaign_id + 1;
    seeds[campaign_id].insert(v);
    seed_list.insert(v);
    num_seeds++;
    change_parameter();
    int new_seed = v;
    for (int i = 0; i < num_sample; i++)
    {
        int x = index[i].x;
        int flag = 0;
        if (seed_bool.find(x) != seed_bool.end())
        {
            if (x == v)
            {
                for (auto s : seed_list)
                    sum_con[s] -= index[i].con[s];
                reconstruct(i);
            }
            continue;
        }
        vector<int> d;
        d.resize(n);
        vector<int> vis;
        vis.resize(n);
        queue<int> Q;
        d[x] = 1;
        Q.push(x);
        vis[x] = 1;
        int seed_dist = inf;
        while (!Q.empty())
        {
            int u = Q.front();
            Q.pop();
            if (seed_dist != inf && d[u] == seed_dist)
                continue;
            for (auto j : re[u])
            {
                int v = j.first;
                double z = j.second;
                if (index[i].dead_edge.find(PII(u, v)) != index[i].dead_edge.end())
                    continue;
                if (index[i].live_edge.find(PII(u, v)) == index[i].live_edge.end())
                {
                    if (check((int)(1.0 / z)))
                        index[i].live_edge.insert(PII(u, v));
                    else
                        index[i].dead_edge.insert(PII(u, v));
                }
                if (index[i].live_edge.find(PII(u, v)) != index[i].live_edge.end())
                {
                    if (vis[v] != 1)
                    {
                        vis[v] = 1;
                        d[v] = d[u] + 1;
                        Q.push(v);
                        if (seed_bool.find(v) != seed_bool.end())
                        {
                            if (v == new_seed)
                                flag = 1;
                            seed_dist = min(seed_dist, d[v]);
                        }
                    }
                }
            }
        }
        if (flag == 1)
        {
            for (auto s : seed_list)
                sum_con[s] -= index[i].con[s];
            reconstruct(i);
        }
    }
}

void Dynamic::change_budget(int b)
{
    if ((int)seed_bool.size() < budget)
    {
        printf("Fail. The budget is larger than the number of seeds.\n");
        return;
    }
    this->budget = b;
    change_parameter();
}