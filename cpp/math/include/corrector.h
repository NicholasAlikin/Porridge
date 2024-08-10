/*
Corrector algorithms:
    o default: Newton-Raphson algorithm
    o normal-flow: Gauss-Newton algorithm
                   using Moore-Penrose inverse
                   (pseudoinverse) matrix
    o Psevdo-arc-lenght algorithms:
        - sphere scheme,
        - orthogonal scheme.
*/

#include "vector.h"

template <typename T>
class Corrector {
using namespace math;
public:
    size_t iter = 0;
    size_t global_iter = 0;
    int exitflag = 0;

    Corrector() {}
    virtual ~Corrector() {}

    // based on corrector algorithm
    vector_t<double,1> process_iteration(vector_t<double,1>, vector_t<double,2>);
    void process_callback();

    


};