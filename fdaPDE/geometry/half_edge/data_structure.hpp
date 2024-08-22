#include <iostream>
#include "../../../eigen/Eigen/Dense"   

class Vertex;
class HalfEdge;
class Edge;

class Vertex{
private:
          Eigen::Vector3d p;
          HalfEdge* e;
public:
          Vertex(const Eigen::Vector3d& point): p(point), e(nullptr){};
          void setPoint(Eigen::Vector3d point){p=point;};
          void setEdge(HalfEdge* ed){e=ed;};

          //destructor
};

class HalfEdge{
private:
          Vertex v;
          HalfEdge* twin;
          HalfEdge* next;
          HalfEdge* prec;
public:
          HalfEdge(const Eigen::Vector3d& ver): v(ver), twin(nullptr), next(nullptr), prec(nullptr){}; 
          HalfEdge();        

          void setVertex(Vertex ver){v=ver;};
          void setTwin(HalfEdge* tw){twin=tw;};
          void setNext(HalfEdge* ne){next=ne;};
          void setPrec(HalfEdge* pr){prec=pr;};

          //destructor
          
};


class Edge{
private:
          HalfEdge *e1, *e2;
public:
          Edge(const Eigen::Vector3d& ver1, const Eigen::Vector3d& ver2): //BOH CON TUTTI I PUNTATORI NON CAPISCO
          {
                    e1->setVertex(ver1); 
                    e2->setVertex(ver2);
          };

//CAPIRE SE PUO ANDARE BENE
//SERVE UNA PARENTELA?
          Edge() {
                    e1 = new HalfEdge();
                    e2 = new HalfEdge();
                    e1->setTwin(e2);
                    e2->setTwin(e1);
          }

          void setEdge1(HalfEdge* ed1){e1=ed1;};
          void setEdge2(HalfEdge* ed2){e2=ed2;};
          ~Edge() {
                    delete e1;
                    delete e2;
          }

};


//FARE TEST SU PARENTELE TWIN, PREC NEXT PER CAPIRE SE VA BENE