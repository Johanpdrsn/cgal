template<class T>
struct CentroidDecomposition {

    std::map<int, std::set<int>> tree;
    std::map<int, int> size;
    int id;

    CentroidDecomposition(const std::vector<T> &data) : tree(create_adjacency_list(data)) {
        int n = tree.size();
        build(tree.begin()->first, -1);
    }

    int dfs(int u, int p) {
        size[u] = 1;
        for (auto v: tree[u])
            if (v != p) {
                size[u] += dfs(v, u);
            }
        return size[u];
    }

    int centroid(int u, int p, int n) {
        for (auto v: tree[u])
            if (v != p) {
                if (size[v] > n / 2) return centroid(v, u, n);
            }
        return u;
    }

    void build(int u, int p) {
        int n = dfs(u, p);
        int c = centroid(u, p, n);
        if (p == -1)
            p = c;
        id = p;
    }

    std::map<int, std::set<int>> create_adjacency_list(const std::vector<T> &data) {

        std::set<int> ids;
        std::for_each(data.begin(), data.end(),
                      [&ids](const T &A) { ids.insert(A->info().id); });

        std::map<int, std::set<int>> t;

        for (const auto &f: data) {
            std::set<int> s;
            for (int i = 0; i < 3; i++) {
                auto neigh = f->neighbor(i);
                if (neigh != nullptr && neigh->info().in_domain() && ids.count(neigh->info().id) > 0) {
                    s.insert(neigh->info().id);
                }
            }
            t[f->info().id] = s;
        }
        return t;
    }
};
