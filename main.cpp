#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include <list>
#include <ctime>

using namespace std;

class Faces {
public:
    vector<vector<int>> interior;
    vector <int> external;
    int n = 0;
    // Проверка полученных граней
    Faces(vector<vector<int>>* interior, vector <int>* external) {
        if(interior != nullptr && external != nullptr) {
            this->interior = *interior;
            this->external = *external;
            n = (*interior).size() + (*external).size();
        } else {
            n = 0;
        }
    }
    // Преобразование данных для вывода
    string toString() {
        string result = "Faces size = " + to_string(n) + "\nExternal face:\n";
        for (int i : external) {
            result += to_string(i) + " ";
        }
        result += "\nInterior faces:\n";
        for (auto & i : interior) {
            for (int j : i) {
                result += to_string(j) + " ";
            }
            result += "\n";
        }

        return result;
    }
};

class Graph
{
public:
    vector<vector<int>> matrix;
    int n = 0;
    // explicit - для предотвращения неявного преобразования типов при инициализации
    explicit Graph (const vector<vector<int>>& m) {
        matrix = m;
        n = (int)m.size();
    }

    explicit Graph (int n) {
        this->n = n;
        matrix.assign(n, vector<int>(n));
    }

    string toString()  {
        string res;
        for(int i = 0;i < n;i++) {
            for(int j = i;j < n;j++) {
                if(matrix[i][j] == 1) {
                    res += to_string(i) + " -- " + to_string(j) + ";" + "\n";
                }
            }
        }
        return res;
    }

    static vector<vector <int>> RemoveAll(vector<vector<int>> *a, const vector<int>& b) { //+
        for (int i = 0; i < a->size(); i++) {
            if ((*a)[i] == b) {
                a->erase(a->begin() + i);
            }
        }
        return *a;
    }

    bool CheckGraphForm() {
        if (matrix.size() > n) {
            return true;
        }
        for (int x = 0; x < n; x++) {
            for (int y = x+1; y < n; y++) {
                if (matrix[x][y] == 0 && matrix[y][x] !=0) {
                    return true;
                }
                if (matrix[x][y]!=0 && matrix[y][x] == 0) {
                    return true;
                }
            }
        }
        vector<int>* c = GetCycle();
        if (c== nullptr) {
            return true;
        }
        return false;
    }

    void NewEdge (int k, int m) {
        matrix [k][m] = matrix [m][k] = 1;
    }

    bool ContEdge (int k, int m) {
        return (matrix[k][m] == 1);
    }

    bool FindCycle (vector <int>* result, vector<int> used, int parent, int v) {
        used.resize(n, 0);
        used[v] = 1;
        for (int i = 0; i < n; i++) {
            if (i == parent) continue;
            if (matrix[v][i] == 0) continue;
            if (used[i] == 0) {
                result->push_back(v);
                if (FindCycle(result, used, v, i)) {
                    // Цикл найден
                    return true;
                }
                else {
                    result->pop_back();
                }
            }
            if (used[i] == 1) {
                result->push_back(v);
                auto cycle = new vector <int>;
                for (int j = 0; j < result->size(); j++) {
                    if ((*result)[j] == i) {
                        for (int k = j; k < result->size(); k++) {
                            (*cycle).push_back((*result)[k]);
                        }
                        result->clear();
                        for (int & l : *cycle) {
                            result->push_back(l);
                        }
                        return true;
                    }
                }
                return true;
            }
        }
        used[v] = 2;
        return false;
    }
    // Возврат найденного цикла
    vector <int>* GetCycle() {
        auto cycle = new vector<int>();
        bool isCycle = FindCycle(cycle, *new vector<int>, -1, 0);
        if(!isCycle) {
            return nullptr;
        }
        else {
            auto* result = new vector<int>();
            for (int v:*cycle) result->push_back(v);
            return result;
        }
    }
    //Поиск сегментов
    void dfsSegments(vector <int> used, vector<bool>* laidVertexes, Graph* result, int v) {
        used[v] = 1;
        for (int i = 0; i < n; i++) {
            if (matrix[v][i]==1){
                result->NewEdge(v,i);
                if (used[i] == 0 && !(*laidVertexes)[i]) dfsSegments(used, laidVertexes, result, i);
            }
        }
    }
    // Возврат сегментов
    vector <Graph>* getSegments (vector<bool>* laidVertexes, vector<vector<bool>>* edges) {
        auto segments = new vector<Graph> ();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (matrix[i][j] == 1 && !(*edges)[i][j] && (*laidVertexes)[i] && (*laidVertexes)[j]) {
                    auto *t = new Graph(n);
                    t->NewEdge(i, j);
                    (*segments).push_back(*t);
                }
            }
        }
        auto* used = new vector<int>(n);
        for (int i = 0; i < n; i++) {
            if ((*used)[i] == 0 && !(*laidVertexes)[i]) {
                auto* res = new Graph(n);
                dfsSegments(*used, laidVertexes, res, i);
                segments->push_back(*res);
            }
        }
        return segments;
    }
    // Поиск цепи
    void dfsChain(vector<int>* used, vector<bool>* laidVertexes, vector<int>* chain, int v) {
        (*used)[v] = 1;
        chain->push_back(v);
        for (int i = 0; i < n; i++) {
            if (matrix[v][i] == 1 && (*used)[i] == 0) {
                if (!(*laidVertexes)[i]) {
                    dfsChain(used, laidVertexes, chain, i);
                } else {
                    chain->push_back(i);
                }
                return;
            }
        }
    }
    // Получение цепей графа
    vector<int>* getChain(vector<bool>* laidVertexes) {
        auto result = new vector<int>();
        for (int i = 0; i < n; i++) {
            if ((*laidVertexes)[i]) {
                bool inGraph = false;
                for (int j = 0; j < n; j++) {
                    if (ContEdge(i, j))
                        inGraph = true;
                }
                if (inGraph) {
                    dfsChain(new vector<int>(n), laidVertexes, result, i);
                    break;
                }
            }
        }
        return result;
    }
    // Укладка цепей графа
    static void layChain(vector<vector<bool>>* result, vector <int> chain, bool cyclic) {
        for (int i = 0; i < chain.size() - 1; i++) {
            (*result)[chain[i]][chain[i + 1]] = true;
            (*result)[chain[i + 1]][chain[i]] = true;
        }
        if (cyclic) {
            (*result)[chain[0]][chain[chain.size() - 1]] = true;
            (*result)[chain[chain.size() - 1]][chain[0]] = true;
        }
    }
    // Проверка наличия у грани нужных сегментов
    bool isFaceContainsSegment(vector <int> face, Graph segment, vector<bool>* laidVertexes) const {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (segment.ContEdge(i, j)) {
                    if (((*laidVertexes)[i] && (find(face.begin(), face.end(), i) == face.end()))
                        || ((*laidVertexes)[j] && (find(face.begin(), face.end(), j) == face.end()))) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    // Подсчет гранй относящихся к сегменту
    vector <int> calcNumOfFacesContainedSegments (vector<vector<int>> *intFaces, const vector<int>& extFace,
                                                  vector <Graph> segments, vector<bool>* laidVertexes,vector<vector<int>>* destFaces) {
        auto count = new vector<int>(segments.size());
        for (int i = 0; i < segments.size(); i++) {
            for (const vector<int>& face: *intFaces) {
                if (isFaceContainsSegment(face, segments[i], laidVertexes)) {
                    (*destFaces)[i] = face;
                    (*count)[i]++;
                }
            }
            if (isFaceContainsSegment(extFace, segments[i], laidVertexes)) {
                (*destFaces)[i] = extFace;
                (*count)[i]++;
            }
        }
        return *count;
    }

    Faces* getPlanarLaying() {
        if (n == 1) {
            auto faces = new vector<vector<int>>();
            auto outerFace = new vector<int>();
            (*outerFace).push_back(0);

            for (int i = 0; i < outerFace->size(); i++) {
                faces->push_back(outerFace[i]);
            }
            for (int i = 0; i < outerFace->size(); i++) {
                faces->push_back(outerFace[i]);
            }
            return new Faces(faces, outerFace);
        }
        vector<int>* c = GetCycle();
        if(c->empty()) {
            return nullptr;
        }

        auto intFaces = new vector<vector<int>>;
        vector<int>* extFace = c;
        intFaces->push_back(*c);
        intFaces->push_back(*extFace);

        auto laidVertexes = new vector<bool> (n);
        auto laidEdges = new vector<vector<bool>> (n, vector<bool> (n));
        for (int i : *c) {
            (*laidVertexes)[i] = true;
        }
        layChain(laidEdges, *c, true);

        while (true) {
            vector<Graph>* segments = getSegments(laidVertexes, laidEdges);
            if (segments->empty()) {
                break;
            }
            auto destFaces = new vector<vector<int>>(segments->size());
            auto count = calcNumOfFacesContainedSegments(intFaces, *extFace, *segments, laidVertexes, destFaces);
            int mi = 0;
            for (int i = 0; i < segments->size(); i++) {
                if (count[i] < count[mi])
                    mi = i;
            }
            if (count[mi] == 0) {
                return nullptr;
            } else { //++++++
                //Укладка выбранного сегмента
                //Выделяем цепь между двумя контактными вершинами

                auto chain = ((*segments)[mi]).getChain(laidVertexes);
                for (int i : *chain) {
                    (*laidVertexes)[i] = true;
                }
                layChain(laidEdges, *chain, false);
                vector <int> face = (*destFaces)[mi];

                auto face1 = new vector<int>();
                auto face2 = new vector<int>();
                int contactFirst = 0, contactSecond = 0;
                for (int i = 0; i < face.size(); i++) {
                    if (face[i] == (*chain)[0]) {
                        contactFirst = i;
                    }
                    if (face[i] == (*chain)[chain->size() - 1]) {
                        contactSecond = i;
                    }
                }
                auto reverseChain = *chain;
                reverse(reverseChain.begin(), reverseChain.end());
                int faceSize = face.size();
                if (face != *extFace) {

                    if (contactFirst < contactSecond) {
                        copy((*chain).begin(), (*chain).end(), back_inserter(*face1));
                        for (int i = (contactSecond + 1) % faceSize; i != contactFirst; i = (i + 1) % faceSize) {
                            face1->push_back(face[i]);
                        }
                        copy(reverseChain.begin(), reverseChain.end(), back_inserter(*face2));
                        for (int i = (contactFirst + 1) % faceSize; i != contactSecond; i = (i + 1) % faceSize) {
                            face2->push_back(face[i]);
                        }
                    } else {
                        copy(reverseChain.begin(), reverseChain.end(), back_inserter(*face1));
                        for (int i = (contactFirst + 1) % faceSize; i != contactSecond; i = (i + 1) % faceSize) {
                            face1->push_back(face[i]);
                        }
                        copy((*chain).begin(), (*chain).end(), back_inserter(*face2));
                        for (int i = (contactSecond + 1) % faceSize; i != contactFirst; i = (i + 1) % faceSize) {
                            face2->push_back(face[i]);
                        }
                    }

                    RemoveAll(intFaces,face);
                    intFaces->push_back(*face1);
                    intFaces->push_back(*face2);


                } else {
                    auto newOuterFace = new vector<int>();
                    if (contactFirst < contactSecond) {
                        copy((*chain).begin(), (*chain).end(), back_inserter(*newOuterFace));
                        for (int i = (contactSecond + 1) % faceSize; i != contactFirst; i = (i + 1) % faceSize) {
                            newOuterFace->push_back(face[i]);
                        }
                        copy((*chain).begin(), (*chain).end(), back_inserter(*face2));
                        for (int i = (contactSecond - 1 + faceSize) % faceSize; i != contactFirst; i = (i - 1 + faceSize) % faceSize) {
                            face2->push_back(face[i]);
                        }
                    } else {
                        copy(reverseChain.begin(), reverseChain.end(), back_inserter(*newOuterFace));
                        for (int i = (contactFirst + 1) % faceSize; i != contactSecond; i = (i + 1) % faceSize) {
                            newOuterFace->push_back(face[i]);
                        }
                        copy(reverseChain.begin(), reverseChain.end(), back_inserter(*face2));
                        for (int i = (contactFirst - 1 + faceSize) % faceSize; i != contactSecond; i = (i - 1 + faceSize) % faceSize) {
                            face2->push_back(face[i]);
                        }
                    }
                    RemoveAll(intFaces, *extFace);
                    intFaces->push_back(*newOuterFace);
                    intFaces->push_back(*face2);
                    extFace = newOuterFace;
                }
                delete face1;
                delete face2;
            }

        }
        return new Faces(intFaces, extFace);
    }
};

int main() {
    ifstream fin("input.txt");
    ofstream fout("visualize.gv");
    unsigned int start_time =  clock();
    int n = 0;
    fin >> n;
    vector<vector<int>> m(n,vector <int> (n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fin >> m[i][j];
        }
    }
    auto* gr = new Graph(m);
    if (gr->CheckGraphForm()) {
        cout << "Incorrect graph format" << endl;
        return 0;
    }
    fout << "graph G {" << endl;

    cout << "Graph:" << endl;
    cout << gr->toString() << endl;

    fout << gr->toString();
    fout << "}";
    Faces* planar = (*gr).getPlanarLaying();
    if(planar != nullptr) {
        cout << "Graph is planarity. Edges:" << endl;
        cout << planar->toString();
    }
    else {
        cout << "Graph is not planarity";
    }
    unsigned int end_time = clock();
    cout << endl << (end_time - start_time)/1000.0 << " sec." << endl;
    fin.close();
    fout.close();
    return 0;
}
