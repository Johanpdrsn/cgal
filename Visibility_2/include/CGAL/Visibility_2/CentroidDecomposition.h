template<class T>
struct CentroidDecomposition {

    std::unordered_map<int, std::unordered_set<int>> tree;
    std::unordered_map<int, int> size;
    int id;
    std::unordered_set<int> left, right;

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

        left.insert(p);
        for (auto v: tree[id]) {
            tree[v].erase(id);
        }

        if (tree[p].size() == 0) {
            return;
        } else if (tree[p].size() == 1) {
            right.insert(*(tree[p].begin()));
        } else {
            left.insert(*tree[p].begin());
            right.insert(*next(tree[p].begin()));

            tree.erase(p);

            find_graph(left);
            find_graph(right);

            left.size() < right.size() ? left.insert(p) : right.insert(p);
        }
    }

    void find_graph(std::unordered_set<int> &s) {
        int temp;
        do {
            temp = s.size();
            for (auto a: s) {
                s.insert(tree[a].begin(), tree[a].end());
            }
        } while (temp != s.size());

    }

    std::unordered_map<int, std::unordered_set<int>> create_adjacency_list(const std::vector<T> &data) {

        std::unordered_set<int> ids;
        std::for_each(data.begin(), data.end(),
                      [&ids](const T &A) { ids.insert(A->info().id); });

        std::unordered_map<int, std::unordered_set<int>> t;

        for (const auto &f: data) {
            std::unordered_set<int> s;
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
