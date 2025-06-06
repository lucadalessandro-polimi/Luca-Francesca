#include <iostream>
#define TRIREAL double
#include <cstring>


// Per usare funzioni C in codice C++
extern "C" {
    #include "triangle.h"
}

int main() {
    struct triangulateio in, out;

    // Zero inizializzazione
    memset(&in, 0, sizeof(in));
    memset(&out, 0, sizeof(out));

    // === INPUT ===
    in.numberofpoints = 3;
    in.pointlist = (TRIREAL *) malloc(2 * 3 * sizeof(TRIREAL));
    in.pointlist[0] = 0.0; in.pointlist[1] = 0.0;
    in.pointlist[2] = 1.0; in.pointlist[3] = 0.0;
    in.pointlist[4] = 0.0; in.pointlist[5] = 1.0;

    in.numberofsegments = 3;
    in.segmentlist = (int *) malloc(2 * 3 * sizeof(int));
    in.segmentlist[0] = 0; in.segmentlist[1] = 1;
    in.segmentlist[2] = 1; in.segmentlist[3] = 2;
    in.segmentlist[4] = 2; in.segmentlist[5] = 0;

    // === CHIAMATA A TRIANGULATE ===
    triangulate(const_cast<char*>("pzQ"), &in, &out, nullptr);

    // === OUTPUT: triangoli ===
    std::cout << "Numero triangoli: " << out.numberoftriangles << std::endl;
    for (int i = 0; i < out.numberoftriangles; ++i) {
        std::cout << "Triangolo " << i << ": ";
        for (int j = 0; j < out.numberofcorners; ++j) {
            std::cout << out.trianglelist[i * out.numberofcorners + j] << " ";
        }
        std::cout << std::endl;
    }

    // === Libera memoria ===
    free(in.pointlist);
    free(in.segmentlist);
    free(out.pointlist);
    free(out.trianglelist);

    return 0;
}