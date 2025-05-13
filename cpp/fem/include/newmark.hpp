#pragma once

#include "linalg.hpp"
#include "fem.hpp"
#include "analyzes.hpp"

#include <array>

namespace fem {

/* Numerical time integration using Newmark method

    1. For corotational finite elements
*/

class Newmark {
public:
    enum flags_time_step { OK, RESIDAL, MAX_TIME_STEP_ITERS };
    
    struct IntegrationResults {
        friend class Newmark;

        math::vector_t<double,2> displacements;         /* nodal displacements for each time */
        math::vector_t<double,2> velocities;      /* nodal velicities for each time */
        math::vector_t<double,2> accelerations;    /* nodal accelerations for each time */
        math::vector_t<double,1> times;             /* integrate time moments */

    private:
    public:
    };

    struct StepControl {
        friend class Newmark;

        double step_max             // max time step
             , step_min             // min time step
             , step_inc;            // incriase step: step *= step_inc
                                    // decriase step: step /= step_inc

        size_t min_subiters         // max time step corrector iters to consider time step `successful` 
             , successful_steps_max // number of `successful` time steps in a row
                                    // after which next `successful` step will increase time step size
             , max_subiters;        // max time step corrector iters
    // private:
        size_t successful_steps;                        // number of successful time steps in a row
        size_t count_point_rejected_time_steps;         // number of rejected time steps for current time moment
        size_t count_point_step_reduction;              // number of time step decriasing for current time moment
        size_t max_point_rejected_time_steps;           // max allowed number of time step rejection.
                                                        // If count_point_rejected_time_steps > max_point_rejected_time_steps
                                                        // time step back wiil be done.
    public:
        StepControl(double step_max, double step_min, double step_inc
                        , size_t min_subiters, size_t successful_steps_max, size_t max_subiters
                        , size_t max_point_rejected_time_steps = 3);
    };

private:

    double dt;      // time step
    double gamma;   // integration parameter
    double beta;    // integration parameter
    double dt_min;  
    double dt_max;                  // max time step
    double dt_inc;                  // time step increment*
    size_t max_substeps;            // max subiters in time step
    size_t min_substeps;            // min subiters in time step to increase step size
    size_t max_step_try_at_point;   // max number of tries to step
    
    math::vector<double> temp_theta;
    math::vector_t<double,2> temp_rotTensor;
    math::vector_t<double,2> temp_Rsumi;
    math::vector_t<double,3> Rsum;         /* nodal total rotation tensors for current time */

    math::vector_t<double,2> displacements;     
    math::vector_t<double,2> velocities;  
    math::vector_t<double,2> accelerations;
    math::vector_t<double,4> rotation_tensors;
    math::vector_t<double,1> times;

    math::vector<double> inc_displacement;
    math::vector<double> displacement; 
    math::vector<double> velocity;
    math::vector<double> acceleration;
    math::vector<double> inc_acceleration;

    math::vector<double> load_external;
    math::vector<double> load_inertia;
    math::vector<double> load_damping;
    math::vector<double> load_inner;
    math::vector<double> load_general;
    math::vector<double> load_reaction;
    

    math::vector_t<double,2> matrix_mass;
    math::vector_t<double,2> matrix_damp;
    math::vector_t<double,2> matrix_gyro;
    math::vector_t<double,2> matrix_stif;
    math::vector_t<double,2> matrix_general; // general * dq = fgeneral

    double damp_mass, damp_stif; // Rayleight damping coefficients

    double gamma__beta_dt
          ,beta_dt_1
          ,beta_dt2_1
          ,one__gemma_beta
          ,dt_gamma__one_2beta_1
          ,one__twobeta_1
          ,one__gamma_twobeta_1
          ,dt__one__gamma_twobeta_1;

    /* Temporal variables to avoid not necessary memory allocations */
    math::vector_t<double,2> temp_L;
    math::vector<double> temp_D;

    /* Method to calculate Newmark integrate constants*/
    double calculate_gamma(double alpha);
    double calculate_beta(double alpha);

    /* Time step iteration process max residals */
    double max_residal_displacement = 1e-5
          ,max_residal_force        = 1e-4
          ,max_residal_work         = 1e-5  // force * displacement
          ,residal_displacement             // current residal
          ,residal_force                    // current residal
          ,residal_work;                    // current residal
    int subiter_flag;     // flag of time step iteration process

    size_t iter, iter_total;       // number of time steps
    size_t subiter;  // number of iteration during time step
    
    StepControl step_ctrl;

    bool incremental_large_rotation;
    /* Store integration results 
    TODO: store in q only linear displacements and not angles
    */


public:
    Newmark() = default;

    Newmark(double alpha
            ,const Model& model, const Assemble& assemble
            ,double damp_mass, double damp_stif, const StepControl& step_ctrl
            ,bool incremental_large_rotation = false);


    void integrate(IntegrationResults& results
                  ,Model& model
                  ,const Assemble& assemble
                  ,const math::vector<double>& displacement0
                  ,const math::vector<double>& velocity0
                  ,const math::vector_t<double,3>& Rsum0
                  ,double start_time
                  ,double end_time
                  ,double first_step
                  ,double max_time_error = 1e-8);
private:
    void integrate_precomputing(const math::vector<double>& displacement0
                              , const math::vector<double>& velocity0
                              , const math::vector_t<double,3>& Rsum0
                              , double start_time
                              , double end_time
                              , double first_step);
    void add_solve(double cur_time);
    void print_time_step(double cur_time, double end_time);

    void time_step(Model& model
                  ,const Assemble& assemble
                  ,double time);
    bool is_correct_time_step();
    void accept_time_step(double cur_time, double end_time);
    double reject_time_step();
    bool increase_time_step();
    bool decrease_time_step();
    void time_step_back();

    void time_step_iteration(const Model& model
                            ,const Assemble& assemble);
    void calc_flag_time_step();
    void print_time_step_iteration();

    void calc_matrix_general();
    void calc_load_general();

    void calc_state();
    void predict_state();
    void inc_state(const Model& model
                         ,const Assemble& assemble);
    void update_step(double new_step);

    void calc_general_load_n_matrix(const Model& model
                                   ,const Assemble& assemble);
    void calc_load_damp();
    void calc_init_acceleration();
    
};



} // namespace fem