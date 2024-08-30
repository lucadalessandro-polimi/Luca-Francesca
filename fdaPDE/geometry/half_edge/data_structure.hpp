#include <Eigen/Dense>
#include <iostream>

template <typename T>
class HalfEdge;

template <typename T>
class Vertex {
private:
    T p;
    HalfEdge<T>* e;
public:
    Vertex(const T& point) : p(point), e(nullptr) {}
    void setPoint(const T& point) { p = point; }
    void setHalfEdge(HalfEdge<T>* edge) { e = edge; }
    T getPoint() const { return p; }
    HalfEdge<T>* getHalfEdge() const { return e; }
};

template <typename T>
class HalfEdge {
private:
    Vertex<T>* v;
    HalfEdge<T>* twin;
    HalfEdge<T>* next;
    HalfEdge<T>* prev;
public:

    //HalfEdge()= default;

    HalfEdge(Vertex<T>* vert = nullptr)
        : v(vert), twin(nullptr), next(nullptr), prev(nullptr) {}

    // Copy constructor
    HalfEdge(const HalfEdge& other) : v(other.v), twin(other.twin), next(other.next), prev(other.prev) {}

    // Copy assignment operator
    HalfEdge& operator=(const HalfEdge& other) {
        if (this != &other) {
            v = other.v;
            twin = other.twin;
            next = other.next;
            prev = other.prev;
        }
        return *this;
    }

    // Move constructor
    HalfEdge(HalfEdge&& other) noexcept 
        : v(other.v), twin(other.twin), next(other.next), prev(other.prev) {
        other.v = nullptr;
        other.twin = nullptr;
        other.next = nullptr;
        other.prev = nullptr;
    }

    // Move assignment operator
    HalfEdge& operator=(HalfEdge&& other) noexcept {
        if (this != &other) {
            v = other.v;
            twin = other.twin;
            next = other.next;
            prev = other.prev;
            other.v = nullptr;
            other.twin = nullptr;
            other.next = nullptr;
            other.prev = nullptr;
        }
        return *this;
    }

    void setVertex(Vertex<T>* vert) { v = vert; }
    void setTwin(HalfEdge<T>* tw) { twin = tw; }
    void setNext(HalfEdge<T>* ne) { next = ne; }
    void setPrev(HalfEdge<T>* pr) { prev = pr; }

    Vertex<T>* getVertex() const { return v; }
    HalfEdge<T>* getTwin() const { return twin; }
    HalfEdge<T>* getNext() const { return next; }
    HalfEdge<T>* getPrev() const { return prev; }

    // Prevent copying to avoid accidental shallow copies
    //HalfEdge(const HalfEdge&) = delete;
    //HalfEdge& operator=(const HalfEdge&) = delete;
};

template <typename T>
class Edge {
private:
    HalfEdge<T>* e1;
    HalfEdge<T>* e2;
public:

    Edge()= default;

    Edge(Vertex<T>* v1, Vertex<T>* v2) {
        e1 = new HalfEdge<T>(v1);
        e2 = new HalfEdge<T>(v2);
        e1->setTwin(e2);
        e2->setTwin(e1);
    }

    // Move constructor
    Edge(Edge&& other) noexcept 
        : e1(other.e1), e2(other.e2) {
        other.e1 = nullptr;
        other.e2 = nullptr;
    }

    // Move assignment operator
    Edge& operator=(Edge&& other) noexcept {
        if (this != &other) {
            delete e1;
            delete e2;
            e1 = other.e1;
            e2 = other.e2;
            other.e1 = nullptr;
            other.e2 = nullptr;
        }
        return *this;
    }

    ~Edge() {
        delete e1;
        delete e2;
    }

    HalfEdge<T>* getEdge1() const { return e1; }
    HalfEdge<T>* getEdge2() const { return e2; }

    // Prevent copying to avoid accidental shallow copies
    Edge(const Edge&) = delete;
    Edge& operator=(const Edge&) = delete;
};

// Funzione per controllare i collegamenti del HalfEdge
template <typename T>
void checkHalfEdgeConnections(const HalfEdge<T>* edge) {
    std::cout << "Twin: " << (edge->getTwin() ? "Correct" : "Incorrect") << std::endl;
    std::cout << "Next: " << (edge->getNext() ? "Correct" : "Not present") << std::endl;
    std::cout << "Prev: " << (edge->getPrev() ? "Correct" : "Not present") << std::endl;
}









