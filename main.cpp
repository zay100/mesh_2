#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <map>

using namespace std;

#ifndef SAFE_MODE
#define SAFE_MODE
#endif

#ifdef SAFE_MODE
#define SAFE_ASSERT(X,...) if(!(X)) exit(printf("Assert failed: "))
#else
#define SAFE_ASSERT(X,...)
#endif

struct Node {
    double x;
    double y;
    Node(double _x, double _y){
        x = _x;
        y = _y;
    }
};

struct Edge {
    unsigned int first, second;
    Edge(unsigned int _first, unsigned int _second){
        first = _first;
        second = _second;
    }
    Edge() {
        first = 0;
        second = 0;
    }
};

struct Boundary {
    Edge edge;
    unsigned int direction;
    Boundary(unsigned int _first, unsigned int _second, unsigned int _direction){
        edge.first = _first;
        edge.second = _second;
        direction = _direction;
    }
};

typedef int tInt; // index

enum tIdxType {
    IDX_N   = 0,  // nodes
    IDX_E   = 1,  // elements
    IDX_S   = 2  // segments
};

template<tIdxType T>
class tIdx {
private:
    tInt i;
    inline tIdx(const tInt j) : i(j) {}

public:
    inline tIdx() : i(-1) {}

    inline void Set(int j) { i = j; }

    inline tInt idx() const { return i; }

    inline tIdx<T>& operator=(const tInt j) { i = j; return *this; }
    inline tIdx<T>& operator=(const size_t j) { i = (tInt)j; return *this; }
    inline tIdx<T>& operator=(const tIdx<T>& j){ i = j.i; return *this; }

    inline tIdx<T>& operator++() { ++i; return *this; }
    inline tIdx<T>& operator--() { --i; return *this; }

    inline tIdx<T> operator++(int) { return tIdx<T>(i++); }//????? почему int в скобке, а не tidx
    inline tIdx<T> operator--(int) { return tIdx<T>(i--); }

    inline tIdx<T> operator+(tIdx<T> j) { return tIdx<T>(i + j.i); }
    inline tIdx<T> operator-(tIdx<T> j) { return tIdx<T>(i - j.i); }

    inline tIdx<T> operator+(tInt j) { return tIdx<T>(i + j); }
    inline tIdx<T> operator-(tInt j) { return tIdx<T>(i + j); }
    inline tIdx<T> operator*(tInt j) { return tIdx<T>(i * j); }
    inline tIdx<T> operator/(tInt j) { return tIdx<T>(i / j); }
    inline tIdx<T> operator%(tInt j) { return tIdx<T>(i % j); }

    inline bool operator<(tInt j) const { return i < j; }
    inline bool operator>(tInt j) const { return i > j; }
    inline bool operator<=(tInt j) const { return i <= j; }
    inline bool operator>=(tInt j) const { return i >= j; }
    inline bool operator==(tInt j) const { return i == j; }
    inline bool operator!=(tInt j) const { return i != j; }

    inline bool operator<(tIdx<T> j) const { return i < j.i; }
    inline bool operator>(tIdx<T> j) const { return i > j.i; }
    inline bool operator<=(tIdx<T> j) const { return i <= j.i; }
    inline bool operator>=(tIdx<T> j) const { return i >= j.i; }
    inline bool operator==(tIdx<T> j) const { return i == j.i; }
    inline bool operator!=(tIdx<T> j) const { return i != j.i; }

    inline bool operator!() { return !i; }

    template<tIdxType TT>
    tIdx(const tIdx<TT> &J);//не понятно
};

template<typename T>
class tBlock {
private:
    T *p;
    const char *name; // pointer to parent object's name
    int N;

    inline tBlock() : N(0), p(nullptr), name(nullptr) {}
    inline tBlock<T>& operator=(const tBlock<T>& B) {
        if (this != &B) { N = B.N; p = B.v; name = B.name; }
        return *this;
    }

public:
    inline tBlock(const tBlock<T> &B) : N(B.N), p(B.p), name(B.name) {}
    inline tBlock(T* ptr, int n, const char* label = nullptr) : p(ptr), N(n), name(label) {}

    inline ~tBlock() {}

    inline void Reset() { N = 0; p = nullptr; name = nullptr; }

    inline operator T* () { return p; }//не понятен синтаксис, почему нет возвращаемого значения типа
    inline operator const T* () const { return p; }

    inline int Size() const { return N; }
    inline T* Body() const { return p; }

    inline bool Allocated() const { return (N > 0) && (p != nullptr); }

    inline T& operator[](int i) {
        SAFE_ASSERT((p && i >= 0 && i < N), "tBlock: wrong index [%d] size=%d (%s)", i, N, name ? name : "");
        return p[i];//посмотреть как работает
    }
    inline const T& operator[](int i) const {
        SAFE_ASSERT((p && i >= 0 && i < N), "tBlock: wrong index [%d] size=%d (%s)", i, N, name ? name : "");
        return p[i];
    }
};

template <tIdxType V, typename T> // T - type of kept value, V - idxType of kept value
struct Topo {
    vector<tInt> IA;
    vector<T> JA;
    string name;

    void clear() {
        IA.clear();
        JA.clear();
    }

    int GetBlockSize(tIdx<V> i) {
        SAFE_ASSERT((i >= 0 && i <= IA.size()), "Topo: wrong index [%d] for block in (%s)", i, N, name ? name : "");
        if (i == IA.size() - 1) {
            return JA.size() - IA[i.idx()];
        }
        return IA[i.idx() + 1] - IA[i.idx()];
    }

    inline tBlock<T> operator[](const tIdx<V> id) {
        SAFE_ASSERT((id >= 0 && id < IA.size()), "Topo: wrong index [%d] blocks_count=%d (%s)", id, IA.size(), name);
        return tBlock<T>(&JA[IA[id.idx()]], GetBlockSize(id));
    }

    T get(tIdx<V> i, int j) {
        return (*this)[i][j];
    }

    void add(vector<T> &a) {
        IA.push_back(JA.size());
        for (int i = 0; i < a.size(); ++i) {
            JA.push_back(a[i]);
        }
    }

    int elem_number() {
        return IA.size();
    }
};

struct Mesh {
public:
    vector<Node> nodes;
    Topo< IDX_E, tIdx<IDX_N> > T; // T == ENadj
    Topo< IDX_N, tIdx<IDX_E> > RT; // RT == NEadj
    Topo< IDX_N, tIdx<IDX_N> > NNadj;
    Topo< IDX_N, tIdx<IDX_S> > NSadj;
    Topo< IDX_S, tIdx<IDX_E> > SEadj;
    Topo< IDX_E, tIdx<IDX_S> > ESadj;
    vector<Boundary> boundary;
    vector<Edge> edges;
    Mesh *child, *father;
    map<unsigned int, unsigned int> nodeConnectionWithChild;
    vector<unsigned int> nodeConnectionWithFather;

    void print_all_debug() {
        cout << "NODES" << endl;
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            cout << nodes[i].x << " " << nodes[i].y << endl;
        }
        cout << "ELEMENTS->NODES" << endl;
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                cout << T[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "NODES->ELEMENTS" << endl;
        tIdx<IDX_N> ne_beg, ne_end;
        ne_beg = 0;
        ne_end = RT.elem_number();
        for (tIdx<IDX_N> i = ne_beg; i < ne_end; ++i) {
            for (int j = 0; j < RT.GetBlockSize(i); ++j) {
                cout << RT[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "Node to Node adjacency" << endl;
        tIdx<IDX_N> nn_beg, nn_end;
        nn_beg = 0;
        nn_end = NNadj.elem_number();
        for (tIdx<IDX_N> i = nn_beg; i < nn_end; ++i) {
            for (int j = 0; j < NNadj.GetBlockSize(i); ++j) {
                cout << NNadj[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "Node to Segment adjacency" << endl;
        tIdx<IDX_N> ns_beg, ns_end;
        ns_beg = 0;
        ns_end = NSadj.elem_number();
        for (tIdx<IDX_N> i = ns_beg; i < ns_end; ++i) {
            for (int j = 0; j < NSadj.GetBlockSize(i); ++j) {
                cout << NSadj[i][j].idx() << " ";
            }
            cout << endl;
        }
        cout << "EDGES" << endl;
        for (unsigned int i = 0; i < edges.size(); ++i) {
            cout << edges[i].first << " " << edges[i].second << endl;
        }
        cout << "BOUNDARY" << endl;
        for (unsigned int i = 0; i < boundary.size(); ++i) {
            cout << boundary[i].edge.first << " " << boundary[i].edge.second << " " << boundary[i].direction << endl;
        }
    }

    void gen_inverse_top () {
        RT.clear();
        vector<tInt> e(nodes.size());
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                ++e[T[i][j].idx()];
            }
        }
        for (unsigned int i = 0; i < e.size(); ++i) {
            tInt curSize = RT.JA.size();
            RT.IA.push_back(curSize);
            tIdx<IDX_E> elem;
            elem = e[i];
            RT.JA.push_back(elem);
            RT.JA.resize(curSize + e[i]);
        }

        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                tIdx<IDX_N> elem = T.get(i, j);
                --RT.JA[RT.IA[elem.idx()]];
                RT.JA[ RT.IA[elem.idx()] + RT.JA[RT.IA[elem.idx()]].idx() ] = i;
            }
        }
    }

    void gen_NNadj() {
        NNadj.clear();
        tIdx<IDX_N> ne_beg, ne_end;
        ne_beg = 0;
        ne_end = RT.elem_number();
        for (tIdx<IDX_N> i = ne_beg; i < ne_end; ++i) {
            unordered_set<int> un;
            vector< tIdx<IDX_N> > v;
            for (int j = 0; j < RT.GetBlockSize(i); ++j) {
                tIdx<IDX_E> elem = RT.get(i, j);
                int elem_size = T.GetBlockSize(elem);
                int in_elem;
                for (int k = 0; k < elem_size; ++k) {
                    if (i == T.get(elem, k)) {
                        in_elem = k;
                        break;
                    }
                }
                for (int k = 0; k < elem_size; ++k) {
                    tIdx<IDX_N> nn = T.get(elem, k);
                    if ((k != in_elem) && (elem_size == 3 || ((k % 3 == 0) ^ (in_elem % 3 == 0)))) {
                        if (un.find(nn.idx()) == un.end()) {
                            un.insert(nn.idx());
                            v.push_back(nn);
                        }
                    }
                }
            }
            NNadj.add(v);
        }
    }

    void gen_NSadj() {
        NSadj.clear();
        vector<tInt> e(nodes.size());
        tIdx<IDX_N> ns_beg, ns_end;
        ns_beg = 0;
        ns_end = NSadj.elem_number();
        for (unsigned int i = 0; i < edges.size(); ++i) {
            ++e[edges[i].first];
            ++e[edges[i].second];
        }
        for (unsigned int i = 0; i < e.size(); ++i) {
            tInt curSize = NSadj.JA.size();
            NSadj.IA.push_back(curSize);
            tIdx<IDX_S> elem;
            elem = e[i];
            NSadj.JA.push_back(elem);
            NSadj.JA.resize(curSize + e[i]);
        }

        for (unsigned int i = 0; i < edges.size(); ++i) {
            unsigned int elem = edges[i].first;
            --NSadj.JA[NSadj.IA[elem]];
            NSadj.JA[ NSadj.IA[elem] + NSadj.JA[NSadj.IA[elem]].idx() ] = i;
            elem = edges[i].second;
            --NSadj.JA[NSadj.IA[elem]];
            NSadj.JA[ NSadj.IA[elem] + NSadj.JA[NSadj.IA[elem]].idx() ] = i;
        }
    }

    void gen_edges() {
        edges.clear();
        tIdx<IDX_N> nn_beg, nn_end;
        nn_beg = 0;
        nn_end = NNadj.elem_number();
        for (tIdx<IDX_N> i = nn_beg; i < nn_end; ++i) {
            for (int j = 0; j < NNadj.GetBlockSize(i); ++j) {
                if (i < NNadj.get(i, j)) {
                    edges.push_back(Edge(i.idx(), NNadj.get(i, j).idx()));
                }
            }
        }
    }

    void gen_ESadj() {
        ESadj.clear();
        vector< tIdx<IDX_S> > ans;
        unordered_set<int> v;
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            v.clear();
            ans.clear();
            for (int j = 0; j < T.GetBlockSize(i); ++j) { // Ð´Ð¾Ð±Ð°Ð²Ð¸Ð¼ Ð²ÑÐµ ÑƒÐ·Ð»Ñ‹ ÑÐ»ÐµÐ¼ÐµÐ½Ñ‚Ð° Ð² ÑÐµÑ‚
                v.insert(T.get(i, j).idx());
            }
            for (int j = 0; j < T.GetBlockSize(i); ++j) {
                tIdx<IDX_N> cur = T.get(i, j); // Ð¿Ð¾Ð»ÑƒÑ‡Ð¸Ð»Ð¸ Ð½Ð¾Ð¼ÐµÑ€ Ñ‚ÐµÐºÑƒÑ‰ÐµÐ³Ð¾ ÑƒÐ·Ð»Ð° ÑÐ»ÐµÐ¼ÐµÐ½Ñ‚Ð°
                for (int k = 0; k < NSadj.GetBlockSize(cur); ++k) {
                    tIdx<IDX_S> c_edge = NSadj.get(cur, k); // Ð½Ð¾Ð¼ÐµÑ€ Ñ‚ÐµÐºÑƒÑ‰ÐµÐ³Ð¾ Ñ€ÐµÐ±Ñ€Ð°
                    if (v.count(edges[c_edge.idx()].first) && v.count(edges[c_edge.idx()].second)) {
                        ans.push_back(c_edge);
                    }
                }
                v.erase(cur.idx());
            }
            ESadj.add(ans);
        }
    }

    void gen_SEadj() {
        SEadj.clear();
        SEadj.IA.resize(edges.size());
        for (unsigned int i = 0; i < edges.size(); ++i) {
            SEadj.IA[i] = 2 * i;
        }
        vector< tIdx<IDX_E> > ja(2 * edges.size());
        tIdx<IDX_E> es_beg, es_end;
        es_beg = 0;
        es_end = ESadj.elem_number();
        for (tIdx<IDX_E> i = es_beg; i < es_end; ++i) { // i-Ñ‹Ð¹ ÑÐ»ÐµÐ¼ÐµÐ½Ñ‚
            for (int j = 0; j < ESadj.GetBlockSize(i); ++j) { // Ð¿Ð¾ ÐµÐ³Ð¾ Ñ€ÐµÐ±Ñ€Ð°Ð¼ Ð¸Ð´ÐµÐ¼
                ja[SEadj.IA[ESadj.get(i, j).idx()]++] = i;
            }
        }
        for (unsigned int i = 0; i < edges.size(); ++i) {
            unsigned int n_size = SEadj.JA.size();
            for (tInt j = 2 * i; j < SEadj.IA[i]; ++j) {
                SEadj.JA.push_back(ja[j]);
            }
            SEadj.IA[i] = n_size;
        }
    }

    void print_coord() {
        ofstream out("coordinate.msh");
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            out << nodes[i].x << " " << nodes[i].y << endl;
        }
        out.close();
    }

    void print_back_topo() {
        ofstream out("backtopo.msh");
        tIdx<IDX_N> ne_beg, ne_end;
        ne_beg = 0;
        ne_end = RT.elem_number();
        for (tIdx<IDX_N> i = ne_beg; i < ne_end; ++i) {
            out << RT.GetBlockSize(i);
            for (int j = 0; j < RT.GetBlockSize(i); ++j) {
                out << " " << RT.get(i, j).idx();
            }
            out << endl;
        }
        out.close();
    }

    void print_mesh() {
        ofstream out("mesh.txt");
        out << "Nn " << nodes.size() << endl << "Nt " << T.elem_number() << endl << "NFaceBC " << boundary.size() << endl << "NumCoords 2" << endl;
        out.close();
    }

    void print_nbf() {
        ofstream out("bctopo.msh");
        for (unsigned int i = 0; i < boundary.size(); ++i) {
            out << 2 << " " << boundary[i].edge.first << " " << boundary[i].edge.second << " " << boundary[i].direction << endl;
        }
        out.close();
    }

    void print_topo() {
        ofstream out("topo.msh");
        tIdx<IDX_E> en_beg, en_end;
        en_beg = 0;
        en_end = T.elem_number();
        for (tIdx<IDX_E> i = en_beg; i < en_end; ++i) {
            out << T.GetBlockSize(i);
            for (int j = 0; j < 3; ++j) {
                out << " " << T.get(i, j).idx();
                if (T.GetBlockSize(i) == 4 && j == 1) {
                    out << " " << T.get(i, 3).idx();
                }
            }
            out << endl;
        }
        out.close();
    }

    void print_ESadj() {
        ofstream out("ESadj.msh");
        tIdx<IDX_E> es_beg, es_end;
        es_beg = 0;
        es_end = ESadj.elem_number();
        for (tIdx<IDX_E> i = es_beg; i < es_end; ++i) {
            out << ESadj.GetBlockSize(i);
            for (int j = 0; j < ESadj.GetBlockSize(i); ++j) {
                out << " " << ESadj.get(i, j).idx();
            }
            out << endl;
        }
        out.close();
    }

    void print_SEadj() {
        ofstream out("SEadj.msh");
        tIdx<IDX_S> se_beg, se_end;
        se_beg = 0;
        se_end = SEadj.elem_number();
        for (tIdx<IDX_S> i = se_beg; i < se_end; ++i) { // Ñƒ Ñ€ÐµÐ±Ñ€Ð° Ð»Ð¸Ð±Ð¾ 2 ÑÐ»ÐµÐ¼ÐµÐ½Ñ‚Ð° Ð»Ð¸Ð±Ð¾ 1
            out << SEadj.GetBlockSize(i);
            for (int j = 0; j < SEadj.GetBlockSize(i); ++j) {
                out << " " << SEadj.get(i, j).idx();
            }
            out << endl;
        }
        out.close();
    }

    void print_all() {
        print_coord();
        print_back_topo();
        print_mesh();
        print_nbf();
        print_topo();
        print_ESadj();
        print_SEadj();
    }

    void gen_all() {
        gen_inverse_top();
        gen_NNadj();
        gen_edges();
        gen_NSadj();
        gen_ESadj();
        gen_SEadj();
    }
};

int meshGen(Mesh& m, string filePath) {
    unsigned int Nx, Ny, K, M;
    double Lx, Ly;
    ifstream fin;
    try {
        fin.open(filePath);
    }
    catch (...) {
        cout << "No such file" << endl;
        return 1;
    }
    fin >> Nx >> Ny >> Lx >> Ly >> K >> M;
    m.child = new Mesh;
    m.child->father = &m;

    //Nodes
    vector< tIdx<IDX_N> > allChildsNodes;
    m.nodes.clear();
    double xSize = Lx / (Nx - 1), ySize = Ly / (Ny - 1);
    for (unsigned int i = 0; i < Ny; ++i) {
        for (unsigned int j = 0; j < Nx; ++j) {
            m.nodes.push_back(Node(j * xSize, i * ySize));
            if (i == 0 || i == Ny - 1|| j == 0 || j == Nx - 1) {
                tIdx<IDX_N> childsNodes;
                childsNodes = m.child->nodes.size();
                int fathersNodes = m.nodes.size() - 1;
                m.child->nodes.push_back(Node(j * xSize, i * ySize));
                allChildsNodes.push_back(childsNodes);
                m.nodeConnectionWithChild[fathersNodes] = childsNodes.idx();
                m.child->nodeConnectionWithFather.push_back(fathersNodes);
            }
        }
    }
    //Elements
    m.T.clear();
    m.child->T.clear();
    m.child->T.add(allChildsNodes);
    unsigned int now = 1;
    bool split = true, direction = true;
    for (unsigned int i = 0; i < Ny - 1; ++i) {
        for (unsigned int j = 0; j < Nx - 1; ++j) {
            vector< tIdx<IDX_N> > tm;
            tIdx<IDX_N> leftUp;
            leftUp = i * Nx + j;
            tIdx<IDX_N> leftDown;
            leftDown = leftUp + Nx;
            if (split) {
                if (direction) {
                    tm = {leftUp, leftDown, leftDown + 1};
                    m.T.add(tm);
                    tm = {leftUp, leftUp + 1, leftDown + 1};
                    m.T.add(tm);
                }
                else {
                    tm = {leftUp, leftUp + 1, leftDown};
                    m.T.add(tm);
                    tm = {leftUp + 1, leftDown, leftDown + 1};
                    m.T.add(tm);
                }
            }
            else {
                tm = {leftUp, leftUp + 1, leftDown, leftDown + 1};
                m.T.add(tm);
            }
            if (split && now % K == 0) {
                direction ^= true;
                now = 0;
                split = false;
            }
            else if (!split && now % M == 0) {
                now = 0;
                split = true;
            }
            ++now;
        }
    }
    //Boundary
    m.boundary.clear();
    for (unsigned int i = 0; i < Ny - 1; ++i) {
        Boundary b1 = Boundary(i * Nx, (i + 1) * Nx, 1), b2 = Boundary((i + 1) * Nx - 1, (i + 2) * Nx - 1, 2);
        m.boundary.push_back(b1);
        m.boundary.push_back(b2);
        m.child->boundary.push_back(b1);
        m.child->boundary.push_back(b2);
    }
    for (unsigned int i = 0; i < Nx - 1; ++i) {
        Boundary b1 = Boundary(i, i + 1, 3), b2 = Boundary(i + (Ny - 1) * Nx, i + 1 + (Ny - 1) * Nx, 4);
        m.boundary.push_back(b1);
        m.boundary.push_back(b2);
        m.child->boundary.push_back(b1);
        m.child->boundary.push_back(b2);
    }
    return 0;
}

int main(int argc, char* argv[]) {
    string s;
    if (argc == 1) {
        s = "in.txt";
        /*
        cout << "Enter path to file, that contains 6 numbers:" << endl;
        cout << "Number of nodes by X" << endl << "Number of nodes by Y" << endl;
        cout << "Width" << endl << "Height" << endl;
        cout << "Number of divided cells" << endl << "Number of undivided cells" << endl;
        return 0;*/
    } else {
        s = argv[1];
    }
    Mesh* mesh = new Mesh;
    if (meshGen(*mesh, s)) {
        return 1;
    }
    mesh->gen_all();
    mesh->print_all();
    mesh->print_all_debug();
    return 0;
}
