#include "fftw_plan_manager.h"
unsigned int fftw_plan_manager::plan_flag = FFTW_ESTIMATE;
std::map<uint64_t, fftw_plan> fftw_plan_manager::r2c_plan_list;
std::map<uint64_t, fftw_plan> fftw_plan_manager::c2r_plan_list;

void fftw_plan_manager::set_flag(unsigned int flag)
{
    plan_flag = flag;
}
void fftw_plan_manager::exec_r2c(uint64_t n, double *input, fftw_complex *output)
{
    fftw_plan plan;
    auto iter = r2c_plan_list.find(n);
    if (iter == r2c_plan_list.end())
    {
        plan = fftw_plan_dft_r2c_1d(n, input, output, plan_flag);
        r2c_plan_list.emplace(n, plan);
    }
    else
    {
        plan = iter->second;
    }
    fftw_execute_dft_r2c(plan, input, output);
}
void fftw_plan_manager::exec_c2r(uint64_t n, fftw_complex *input, double *output)
{
    fftw_plan plan;
    auto iter = c2r_plan_list.find(n);
    if (iter == c2r_plan_list.end())
    {
        plan = fftw_plan_dft_c2r_1d(n, input, output, plan_flag);
        c2r_plan_list.emplace(n, plan);
    }
    else
    {
        plan = iter->second;
    }
    fftw_execute_dft_c2r(plan, input, output);
}
void fftw_plan_manager::release_all_plan()
{
    for (auto i : r2c_plan_list)
    {
        fftw_destroy_plan(i.second);
    }
    for (auto i : c2r_plan_list)
    {
        fftw_destroy_plan(i.second);
    }
    r2c_plan_list.clear();
    c2r_plan_list.clear();
}