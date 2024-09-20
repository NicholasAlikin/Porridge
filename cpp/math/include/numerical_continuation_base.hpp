#pragma once

namespace math {
    
struct EXITFLAG {
    static const int OK = 0;
    static const int NORM_VAR_AND_FUN = 1;
    static const int MAX_ITER = 2;
    static const int PARAM_END = 3;

};

struct EXEFLAG {
    static constexpr char HELP[] {"-h"};
    static constexpr char STEP[] {"-ds"};
    static constexpr char PRINT[] {"-print"};
    static constexpr char PARAM_START[] {"-param0"};
    static constexpr char PARAM_END[] {"-param1"};
    
};

void print_exe_help() {
    std::cout
        << "Usage: executive_file [options]\n"
        << "Options:"
        << "\n\t-h\t\tDisplay this information"
        << "\n\t-ds\t\tSet continuation path step"
        << "\n\t-param0\t\tSet begin parameter value"
        << "\n\t-param1\t\tSet end parameter value"
        << std::endl;
}

bool parse_program_options(int argc, char* argv[]
                            ,double& ds
                            ,double& param_start
                            ,double& param_end) {
    if (argc == 1)
        throw std::logic_error("Incorrect number of options!");
    if (std::string(argv[1]) == EXEFLAG::HELP) {
        print_exe_help();
        return false; 
    }
    for (int i = 1; i < argc-1; ++i) {
        if (std::string(argv[i]) == EXEFLAG::STEP) {
            ds = std::stod(argv[++i]);
        } else if (std::string(argv[i]) == EXEFLAG::PARAM_START) {
            param_start = std::stod(argv[++i]);
        } else if (std::string(argv[i]) == EXEFLAG::PARAM_END) {
            param_end = std::stod(argv[++i]);
        } else {
            continue;
        }
    }
    return true;
}



} // namespace math