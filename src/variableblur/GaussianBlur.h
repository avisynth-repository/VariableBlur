
#include "BinomialBlur.h"

typedef struct fftwf_plan_s *fftwf_plan;
typedef float fftw_complex[2];
typedef fftwf_plan (*fftwf_plan_dft_r2c_2d_proc) (int n0, int n1, 
	float *in, fftw_complex *out, unsigned flags);
typedef fftwf_plan (*fftwf_plan_dft_c2r_2d_proc) (int n0, int n1, 
	fftw_complex *in, float *out, unsigned flags);
typedef void (*fftwf_execute_plan_proc) (fftwf_plan);
typedef void (*fftwf_destroy_plan_proc) (fftwf_plan);
typedef int (*fftwf_init_threads_proc)(void);
typedef void (*fftwf_plan_nthreads_proc)(int nthreads);
#define FFTW_DESTROY_INPUT (1U << 0)
#define FFTW_MEASURE (0U)

const double gamma_table[8][5] =
{
	// lv, av, pv, thresh1, thresh2
	12.92,0.055,1.0/2.4,0.0031308,0.04045,	// sRGB
	4.5,0.099,1.0/(1.0/0.45),0.018,0.081,	// BT.709, SMPTE 170M
	4.0,0.1115,1.0/(1.0/0.45),0.0228,0.0912,// SMPTE 240M
	1.0,0.0,1.0/2.2,0.0,0.0,				// BT.470-2 System M   (2.2 gamma, no linear segment)
	1.0,0.0,1.0/2.8,0.0,0.0,				// BT.470-2 System B,G (2.8 gamma, no linear segment)
	1.0,0.0,0.45,0.0,0.0,					// general purpose, (2.22222 gamma, no linear segment)
	1.0,0.0,1/1.8,0.0,0.0,					// general purpose, (1.8 gamma, no linear segment)
	1.0,0.0,1.0,0.0,0.0						// linear, no gamma compensation
};

class GausBlur {
public:
	GausBlur(int width, int height, double sd, int border, bool integrate, double unsharpstrength,
		int _gfunc, int _pcr, int _nthreads, IScriptEnvironment* env);
	~GausBlur();
	void Blur(const unsigned char *srcp, const int spitch, unsigned char *dstp, 
		const int dpitch, IScriptEnvironment* env);
protected:
	void CalcGaus(const bool integrate, const double strength, IScriptEnvironment *env);
	double GausBlur::RombergIntegral(double a, double b, double sd);
	HINSTANCE hLib;
	fftwf_plan_dft_r2c_2d_proc fftwf_plan_dft_r2c_2d;
	fftwf_plan_dft_c2r_2d_proc fftwf_plan_dft_c2r_2d;
	fftwf_execute_plan_proc fftwf_execute;
	fftwf_destroy_plan_proc fftwf_destroy;
	fftwf_init_threads_proc fftwf_init_threads;
	fftwf_plan_nthreads_proc fftwf_plan_nthreads;
	int height;
	int width;
	double sd;
	int border;
	int radius;
	float *in;
	float *realout;
	float norm;
	float g_table[256];
	int gfunc, pcr;
	float *out;
	int outwidth;
	int srcwidth;
	int srcheight;
	float* gaus;
	fftwf_plan plan, plani;
	int nthreads;
};

class GaussianBlur : public GenericVideoFilter {
public:
	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
	GaussianBlur(PClip clip, double sdy, double sdc, int border, bool integrate, double unsharpstrength,
		int _Y, int _U, int _V, int gfunc, int gfuncc, int pcr, int pcrc, int _nthreads, IScriptEnvironment* env);
	~GaussianBlur();
protected:
	GausBlur *gbY;
	GausBlur *gbC;
	int Y, U, V, nthreads;
	unsigned char *tempA, *tempB;
};