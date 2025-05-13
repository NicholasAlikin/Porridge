#include "newmark.hpp"

#ifndef PRINT_TIME_STEP_UPDATE_INFO
#define PRINT_TIME_STEP_UPDATE_INFO
#endif

namespace fem {

Newmark::StepControl::StepControl(double step_max
                                        , double step_min
                                        , double step_inc
                                        , size_t min_subiters
                                        , size_t max_subiters
                                        , size_t successful_steps_max
                                        , size_t max_point_rejected_time_steps)
        : step_max(step_max), step_min(step_min), step_inc(step_inc)
        , min_subiters(min_subiters)
        , successful_steps_max(successful_steps_max), max_subiters(max_subiters)
        , max_point_rejected_time_steps(max_point_rejected_time_steps)
{

}

Newmark::Newmark(double alpha
                ,const Model& model, const Assemble& assemble
                ,double damp_mass, double damp_stif
                ,const StepControl& step_ctrl
                ,bool incremental_large_rotation)
        :gamma(calculate_gamma(alpha))
        ,beta(calculate_beta(beta))
        
        ,temp_theta(BaseNode::DIM)
        ,temp_rotTensor(math::zeros<double>(BaseNode::DIM,BaseNode::DIM))
        ,temp_Rsumi(math::zeros<double>(BaseNode::DIM,BaseNode::DIM))
        ,Rsum(AnalysisTraits::setup_Rsum(model))
        
        ,inc_displacement(assemble.ndofs)
        ,displacement(assemble.ndofs)
        ,velocity(assemble.ndofs)
        ,acceleration(assemble.ndofs)
        
        ,load_external(assemble.ndofs)
        ,load_inertia(assemble.ndofs)
        ,load_damping(assemble.ndofs)
        ,load_inner(assemble.ndofs)
        ,load_general(assemble.ndofs)
        ,load_reaction(assemble.ndofs)

        ,matrix_mass(math::zeros<double>(assemble.ndofs,assemble.ndofs))
        ,matrix_damp(math::zeros<double>(assemble.ndofs,assemble.ndofs))
        ,matrix_gyro(math::zeros<double>(assemble.ndofs,assemble.ndofs))
        ,matrix_stif(math::zeros<double>(assemble.ndofs,assemble.ndofs))
        ,matrix_general(math::zeros<double>(assemble.ndofs,assemble.ndofs))
        
        ,damp_mass(damp_mass), damp_stif(damp_stif)

        ,one__gemma_beta(1.0 - gamma/beta)
        ,one__twobeta_1(1.0 - 1.0/(2.0*beta))
        ,one__gamma_twobeta_1(1.0 - gamma/(2.0*beta))

        ,temp_L(math::zeros<double>(assemble.ndofs,assemble.ndofs))
        ,temp_D(assemble.ndofs)

        ,step_ctrl(step_ctrl)
        ,incremental_large_rotation(incremental_large_rotation)
{
}

double Newmark::calculate_gamma(double alpha) {
    return  0.5 + alpha;
}

double Newmark::calculate_beta(double alpha) {
    return 0.25 * (1.0+alpha)*(1.0+alpha);
}


void Newmark::integrate(Newmark::IntegrationResults& results
                        ,Model& model
                        ,const Assemble& assemble
                        ,const math::vector<double>& displacement0
                        ,const math::vector<double>& velocity0
                        ,const math::vector_t<double,3>& Rsum0
                        ,double cur_time // start time
                        ,double end_time
                        ,double first_step
                        ,double max_time_error)
{
    integrate_precomputing(displacement0,velocity0,Rsum0,cur_time,end_time,first_step);
    
    /* Calculate inital acceleration */
    ModelTraits::store_load_vector(model,assemble,load_external,cur_time);
    calc_general_load_n_matrix(model,assemble);
    calc_init_acceleration();
    print_time_step_iteration();
    
    /* Save inital state */
    accept_time_step(cur_time, end_time);
    
    /* Main loop over time span*/
    while (cur_time < end_time) {
        ++iter_total;

        /* Update current time */
        cur_time += dt;
        if (cur_time > end_time) {
            dt = end_time - (cur_time-dt);
            std::cout << "#Last time step, dt = " << dt << "\n";
            if (dt < max_time_error) {
                break;
            }
            update_step(dt);
            cur_time = end_time;
        }

        /* Predict state for current time */    
        predict_state();
        
        /* Correct predicted state for current time */
        time_step(model,assemble,cur_time);

        /* If state was not corrected step is rejected */
        if (!is_correct_time_step()) {
            cur_time = reject_time_step();
            continue;
        }
        
        /* Accept corrected state */
        accept_time_step(cur_time, end_time);
        ++iter;
    }

    /* Place results */
    results.displacements = std::move(displacements);
    results.velocities = std::move(velocities);
    results.accelerations = std::move(accelerations);
    results.times = std::move(times);
}

void Newmark::integrate_precomputing(const math::vector<double>& displacement0
                                    ,const math::vector<double>& velocity0
                                    ,const math::vector_t<double,3>& Rsum0
                                    ,double start_time
                                    ,double end_time
                                    ,double first_step)
{
    iter = 0;
    std::copy(displacement0.begin(),displacement0.end(),displacement.begin());
    std::copy(velocity0.begin(),velocity0.end(),velocity.begin());
    if (incremental_large_rotation)
        std::copy(Rsum0.begin(),Rsum0.end(),Rsum.begin());
    std::fill(acceleration.begin(),acceleration.end(),0.0);

    update_step(first_step);
}



void Newmark::time_step(Model& model
                       ,const Assemble& assemble
                       ,double time) {
    subiter = 0;
    // external load containts constant while one time step
    ModelTraits::store_load_vector(model,assemble,load_external,time);

    do {
        time_step_iteration(model,assemble);
        calc_flag_time_step();
        inc_state(model,assemble);

        ++subiter;
        print_time_step_iteration();

    } while (subiter_flag == flags_time_step::OK);
}

bool Newmark::is_correct_time_step() {
    /* Time step is `correct` only
    if correction was stopped because
    residal constraints was satisfied */
    if (subiter_flag == flags_time_step::RESIDAL) {
        return true;
    }
    return false;
}

void Newmark::accept_time_step(double cur_time
                             , double end_time) {
    /* Save state for current time */
    add_solve(cur_time);
    print_time_step(cur_time,end_time);
    step_ctrl.count_point_rejected_time_steps = 0;
    step_ctrl.count_point_step_reduction = 0;
    increase_time_step();
}

double Newmark::reject_time_step() {
    ++step_ctrl.count_point_rejected_time_steps;
    if (step_ctrl.count_point_rejected_time_steps >= step_ctrl.max_point_rejected_time_steps) {
        time_step_back();
        step_ctrl.count_point_step_reduction = 0;
        step_ctrl.count_point_rejected_time_steps = 0;
    
    } else {
        if (decrease_time_step()) {
            ++step_ctrl.count_point_step_reduction;
        }
    }
    return times.last();
}

bool Newmark::increase_time_step() {
    if (dt >= step_ctrl.step_max || subiter > step_ctrl.min_subiters) {
        step_ctrl.successful_steps = 0;
        return false;
    }

    if (step_ctrl.successful_steps < step_ctrl.successful_steps_max) {
        ++step_ctrl.successful_steps;
        return false;
    }

    update_step(dt * step_ctrl.step_inc);
    step_ctrl.successful_steps = 0;

#ifdef PRINT_TIME_STEP_UPDATE_INFO
    std::cout << "# Time step is increased! New step size: "
              << dt << '\n';
#endif
    return true;
}

bool Newmark::decrease_time_step() {
    step_ctrl.successful_steps = 0;

    if (dt <= step_ctrl.step_min) return false;

    update_step(dt / step_ctrl.step_inc);
#ifdef PRINT_TIME_STEP_UPDATE_INFO
    std::cout << "# Time step is decreased! New step size: "
              << dt << '\n';
#endif
    return true;
}

void Newmark::time_step_back() {
    /* Delete last accepted time step */
    displacements.erase(displacements.end()-1);
    velocities.erase(velocities.end()-1);
    accelerations.erase(accelerations.end()-1);
    times.erase(times.end()-1);

    /* Update state to new current time */
    displacement = displacements.last();
    velocity = velocities.last();
    acceleration = accelerations.last();
    
    if (incremental_large_rotation) {
        rotation_tensors.erase(rotation_tensors.end()-1);
        Rsum = rotation_tensors.last();
    }
    
    /* Update time step size */
    if (step_ctrl.count_point_step_reduction > 1) {
        update_step(dt * math::pow(step_ctrl.step_inc,step_ctrl.count_point_step_reduction-1));
    }
    

    std::cout << "# Time step back! Current time: "
              << times.last() << '\n';
}

/* Calculate general load and matrix (combination of tangent matrices)
    and solve linear equation */
void Newmark::time_step_iteration(const Model& model
                                 ,const Assemble& assemble)
{
    calc_general_load_n_matrix(model,assemble);
    // math::solve_ldlt(matrix_general,load_general,inc_displacement,temp_L,temp_D,load_general.size());
    math::solve_lu(matrix_general,load_general,inc_displacement,load_general.size());
}

void Newmark::calc_flag_time_step() {
    residal_displacement = math::norm(inc_displacement);
    residal_force = math::norm(load_general);///math::norm(load_reaction);
    residal_work = math::dot(load_general,inc_displacement);
    if ((residal_force < max_residal_force)
        && (residal_work < max_residal_work)
        && (residal_displacement < max_residal_displacement))
    {
        subiter_flag = flags_time_step::RESIDAL;
    } else if (subiter >= step_ctrl.max_subiters) {
        subiter_flag = flags_time_step::MAX_TIME_STEP_ITERS;
    } else {
        subiter_flag = flags_time_step::OK;
    }
}

/* Calculate general load and general matrix */
void Newmark::calc_general_load_n_matrix(const Model& model
                                        ,const Assemble& assemble)
{
    /* Calculate tangent mass matrix and inertia load */
    math::fill(matrix_mass.begin(),matrix_mass.end(),0.0);
    math::fill(load_inertia.begin(),load_inertia.end(),0.0);
    /* Function to calculate element matrix and load vector */
#if 0
    ModelTraits::assemble(model,assemble,matrix_mass,load_inertia
                                ,displacement,velocity,acceleration,Rsum
                                ,&BaseElement::tangentMass_inertiaLoad);
#else
    math::fill(matrix_gyro.begin(),matrix_gyro.end(),0.0);
    if (incremental_large_rotation) {
        ModelTraits::assemble(model,assemble,matrix_mass,matrix_gyro,load_inertia
                                    ,displacement,velocity,acceleration,Rsum
                                    ,&BaseElement::tangentMassGyro_inertiaLoad);
    } else {
        ModelTraits::assemble(model,assemble,matrix_mass,matrix_gyro,load_inertia
                                    ,displacement,velocity,acceleration
                                    ,&BaseElement::tangentMassGyro_inertiaLoad);
    }
#endif
    /* Calculate tangent stiffness matrix and inner load */
    math::fill(matrix_stif.begin(),matrix_stif.end(),0.0);
    math::fill(load_inner.begin(),load_inner.end(),0.0);
    if (incremental_large_rotation) {
        ModelTraits::assemble(model,assemble,matrix_stif,load_inner
                                    ,displacement,Rsum
                                    ,&BaseElement::tangentStiffness_innerLoad);
    } else {
        ModelTraits::assemble(model,assemble,matrix_stif,load_inner
                                    ,displacement
                                    ,&BaseElement::tangentStiffness_innerLoad);
    }
    /* Rayleigh damping is used,
    so tangent damp matrix and damp load is not calculated
    by assemble elements matrices and loads*/
    
    /* Calculate general matrix as linear combination of mass and stiffness matrices.
    Damp matrix is calculated here using Rayleigh damping */
    calc_matrix_general();

    /* Calculate damping load as dot product of damp matrix and velocity */
    calc_load_damp();

    /* Calculate general load as linear combination of external, inertia, inner and damping loads */
    calc_load_general();
}

void Newmark::add_solve(double cur_time)
{
    displacements.push_back(displacement);
    velocities.push_back(velocity);
    accelerations.push_back(acceleration);
    
    if (incremental_large_rotation)
        rotation_tensors.push_back(Rsum);

    times.push_back(cur_time);
}

void Newmark::calc_matrix_general() {
    /* general = mass/(beta*dt^2) + damp*gamma/(beta*dt) + stif */
    
    auto mass_row = matrix_mass.begin()
        ,damp_row = matrix_damp.begin()
        ,gyro_row = matrix_gyro.begin()
        ,stif_row = matrix_stif.begin();
    decltype(matrix_mass[0].begin()) mass_col
                            , damp_col
                            , gyro_col
                            , stif_col;
    auto general_row = matrix_general.begin()
        ,general_row_end = matrix_general.end();
    decltype(matrix_general[0].begin()) general_col
                                       ,general_col_end;

    while (general_row < general_row_end) {
        general_col = general_row->begin();
        general_col_end = general_row->end();
        mass_col = mass_row->begin();
        damp_col = damp_row->begin();
        gyro_col = gyro_row->begin();   
        stif_col = stif_row->begin();

        while (general_col < general_col_end) {
            *damp_col = (*mass_col) * damp_mass
                        +(*stif_col) * damp_stif;
            
            /* If used damp matrix */
            *general_col = (*mass_col)*beta_dt2_1
                         + (*damp_col + *gyro_col)*gamma__beta_dt
                         + (*stif_col);

            /* If used Rayleigh damping */
            // *general_col = (*mass_col)*(beta_dt2_1 + gamma__beta_dt*damp_mass)
                        //  + (*stif_col)*(1.0        + gamma__beta_dt*damp_stif);

            ++general_col;
            ++mass_col;
            ++damp_col;
            ++gyro_col;
            ++stif_col;
        }

        ++general_row;
        ++mass_row;
        ++damp_row;
        ++gyro_row;
        ++stif_row;
    }
    
}


void Newmark::calc_load_general() {
    auto general        = load_general.begin()
        ,general_end    = load_general.end()
        ,freact         = load_reaction.begin() 
        ,fext           = load_external.begin()
        ,finner         = load_inner.begin()
        ,finert         = load_inertia.begin()
        ,fdamp          = load_damping.begin();

        while (general < general_end) {
            *freact = *finner + *finert + *fdamp;
            *general = *fext - *freact;

            ++general;
            ++freact;
            ++fext;
            ++finner;
            ++finert;
            ++fdamp;
        }
}

/* Calculate dqidt and d2qidt2 using qi*/
void Newmark::calc_state() {
    auto displ      = displacement.begin()
        ,displ_end  = displacement.end()
        ,displ_prev = displacements.last().begin()
        ,vel        = velocity.begin()
        ,vel_prev   = velocities.last().begin()
        ,accel      = acceleration.begin()
        ,accel_prev = accelerations.last().begin();

    double sub_displs;
    while (displ < displ_end) {
        sub_displs = *displ - *displ_prev;

        *vel   = gamma__beta_dt*sub_displs + one__gemma_beta*(*vel_prev) + dt__one__gamma_twobeta_1*(*accel_prev);
        *accel = beta_dt2_1    *sub_displs - beta_dt_1      *(*vel_prev) + one__twobeta_1          *(*accel_prev);

        ++displ;
        ++displ_prev;
        ++vel;
        ++vel_prev;
        ++accel;
        ++accel_prev;
    }
}

/* Calculate state prediction dqidt and d2qidt2 assuming q(i) = q(i-1)*/
void Newmark::predict_state() {
#if 1
    auto vel        = velocity.begin()
        ,vel_end    = velocity.end()
        ,vel_prev   = velocities.last().begin()
        ,displ      = displacement.begin()
        // ,accel      = acceleration.begin()
        ,accel_prev = accelerations.last().begin();

    while (vel < vel_end) {

        // *vel   =   one__gemma_beta*(*vel_prev) + dt__one__gamma_twobeta_1*(*accel_prev);
        // *accel =  -beta_dt_1      *(*vel_prev) + one__twobeta_1          *(*accel_prev);
        *displ += dt* (*vel_prev) + dt*dt/2.0*(*accel_prev);
        *vel   += dt* (*accel_prev);
        
        ++vel;
        ++vel_prev;
        ++displ;
        // ++accel;
        ++accel_prev;
    }
    // std::cout << "predicted"
    //     <<  " |q| = " << math::norm(displacement)
    //     <<  ", |v| = " << math::norm(velocity)
    //     <<  ", |a| = " << math::norm(acceleration) << std::endl;
#endif
}

void Newmark::inc_state(const Model& model
                       ,const Assemble& assemble)
{
    auto inc_displ      = inc_displacement.begin()
        ,inc_displ_end  = inc_displacement.end()
        ,displ          = displacement.begin()
        ,vel            = velocity.begin()
        ,accel          = acceleration.begin();
    double dq;
    while (inc_displ < inc_displ_end) {
        dq = *inc_displ;

        *displ += dq;
        *vel   += gamma__beta_dt*dq;
        *accel += beta_dt2_1    *dq;

        ++inc_displ;
        ++displ;
        ++vel;
        ++accel;
    }

    if (incremental_large_rotation) {
        /* Update Rsum after increment velocity and acceleration,
            because they don`t have to be less than 2pi */
        AnalysisTraits::update_Rsum( model
                                    ,assemble
                                    ,Rsum
                                    ,displacement
                                    ,temp_theta
                                    ,temp_rotTensor
                                    ,temp_Rsumi);
    }
}

void Newmark::update_step(double new_step) {
    dt = new_step;

    gamma__beta_dt = gamma/(beta*dt);
    beta_dt_1 = 1./(beta*dt);
    beta_dt2_1 = beta_dt_1/dt;
    dt_gamma__one_2beta_1 = dt*gamma*(one__twobeta_1);
    dt__one__gamma_twobeta_1 = dt*one__gamma_twobeta_1;

}

/* Calculate damping load using Rayleigh damping */
void Newmark::calc_load_damp() {
    /* fdamp = matrix_damp * velocity */
    math::fill(load_damping.begin(), load_damping.end(), 0.0);
    math::dot(matrix_damp,velocity,load_damping);
}

/* Calculate inital acceleration using inital displacement and velocity*/
void Newmark::calc_init_acceleration() {
    math::solve_ldlt(matrix_mass,load_general,acceleration,temp_L,temp_D,load_general.size());
}

void Newmark::print_time_step(double cur_time, double end_time) {
    std::cout << "#Iter " << iter
            << ", time " << cur_time << "/" << end_time
            << ", time step " << dt
            << ", |fext| = " << math::norm(load_external)
            <<  "\n";
#if 0
    std::cout << cur_time << ' ' << displacement << '\n'; 
#endif
}

void Newmark::print_time_step_iteration() {
    std::cout <<"#\tsub iter " << subiter
            <<  ", |dq| = " << residal_displacement
            <<  ", |f| = " << residal_force
            <<  ", |f*dq| = " << residal_work
#if 0
            <<  ", |q| = " << math::norm(displacement)
            <<  ", |v| = " << math::norm(velocity)
            <<  ", |a| = " << math::norm(acceleration)
#endif
#if 0
            <<  ", |finert| = " << math::norm(load_inertia)
            <<  ", |finner| = " << math::norm(load_inner)
            <<  ", |fdamp| = " << math::norm(load_damping)
            <<  ", |fext| = " << math::norm(load_external)
            <<  ", |M-MT| = " << math::norm(matrix_mass-math::transpose(matrix_mass))
            <<  ", |K-KT| = " << math::norm(matrix_stif-math::transpose(matrix_stif))
            <<  ", |G-GT| = " << math::norm(matrix_gyro-math::transpose(matrix_gyro))
#endif
            <<  "\n";
}

} // namespace fem