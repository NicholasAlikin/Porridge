#pragma once

#include <iostream>
#include <string>

namespace npath {
    
struct EXITFLAG {
    static const int OK = 0;
    static const int NORM_VAR_AND_FUN = 1;
    static const int MAX_ITER = 2;
    static const int PARAM_END = 3;
    static const int PREV_POINT = 4;

};

struct EXEFLAG {
    static constexpr char HELP[] {"-h"};
    static constexpr char STEP[] {"-ds"};
    static constexpr char PRINT[] {"-print"};
    static constexpr char PARAM_START[] {"-p0"};
    static constexpr char PARAM_END[] {"-p1"};
    
};

void print_exe_help();

int parse_program_options(int argc, char* argv[]
                            ,double& ds
                            ,double& param_start
                            ,double& param_end);

struct PathFollowing {
    static constexpr double track_backward = 0.9;
};


/*ifdefs:

o CORRECTOR_PRINTITER - \numerical_path_following\source\corrector.hpp


*/


} // namespace npath