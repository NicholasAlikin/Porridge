#include "vector.h"

int main () {
    using namespace math;

    vector_t<int,2> vv = repmat2<int>(3,4,5);
    std::cout << vv;

}