#pragma once

#include "model.hpp"
#include "elements.hpp"
#include "nodes.hpp"

namespace fem {


struct Point;
struct Line;


struct Point {
    static size_t EmptyID;

};


struct StateQ {
	math::vector<double> displacement;
};

struct StateQV: StateQ {
	math::vector<double> velocity;
};

struct StateQVA: StateQV {
	math::vector<double> acceleration;
};

} // namespace fem