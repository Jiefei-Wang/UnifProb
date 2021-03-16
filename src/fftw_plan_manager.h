#include<map>
#include<fftw3.h>
class fftw_plan_manager
{
	static std::map<uint64_t, fftw_plan> r2c_plan_list;
	static std::map<uint64_t, fftw_plan> c2r_plan_list;
	static unsigned int plan_flag;
public:
	static void set_flag(unsigned int flag);
	static void exec_r2c(uint64_t n, double *input, fftw_complex *output);
	static void exec_c2r(uint64_t n, fftw_complex *input, double *output);
	static void release_all_plan();
};