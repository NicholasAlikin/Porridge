#include "numerical_continuation_base.hpp"

namespace npath {

void print_exe_help() {
    std::cout
        << "Usage: executive_file [options]\n"
        << "Options:"
        << "\n\t-h\t\tDisplay this information"
        << "\n\t-ds\t\tSet continuation path step"
        << "\n\t-p0\t\tSet begin parameter value"
        << "\n\t-p1\t\tSet end parameter value"
        << std::endl;
}

int parse_program_options(int argc, char* argv[]
                            ,double& ds
                            ,double& param_start
                            ,double& param_end)
{
    if (argc == 1)
        throw std::logic_error("Incorrect number of options!");
    if (std::string(argv[1]) == EXEFLAG::HELP) {
        print_exe_help();
        return 0; 
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
    return 1;
}

} // namespace npath