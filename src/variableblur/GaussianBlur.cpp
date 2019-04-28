
#include "GaussianBlur.h"
#include <omp.h>

const float PI=3.14159265358979323846;

int mod16(int n)
{
	if (n&15)
		return n+16-(n&15);
	return n;
}

GaussianBlur::GaussianBlur(PClip clip, double sdy, double sdc, int border, bool integrate, double unsharpstrength,
	int _Y, int _U, int _V, int gfunc, int gfuncc, int pcr, int pcrc, int _nthreads, IScriptEnvironment* env):
	GenericVideoFilter(clip), Y(_Y), U(_U), V(_V), nthreads(_nthreads)
{
	gbY = gbC = NULL;
	tempA = tempB = NULL;
	if (nthreads < 1)
		env->ThrowError("gaussianblur:  nthreads must be >= 1!");
	if (!(env->GetCPUFlags()&CPUF_SSE2))
		env->ThrowError("gaussianblur:  sse2 support required!");
	if (!vi.IsYV12() && !vi.IsYUY2() && !vi.IsRGB24() && !vi.IsRGB32())
		env->ThrowError("gaussianblur:  input must be yv12, yuy2, rgb24, or rgb32!");
	if (border < -255 || border > 4)
		env->ThrowError("gaussianblur:  border must be in the range [-255,4]!");
	if (gfunc < -1 || gfunc > 7)
		env->ThrowError("gaussianblur:  gfunc must be in the range [-1,7]!");
	if (gfuncc < -1 || gfuncc > 7)
		env->ThrowError("gaussianblur:  gfuncc must be in the range [-1,7]!");
	if (pcr < 0 || pcr > 2)
		env->ThrowError("gaussianblur:  pcr must be in the range [0,2]!");
	if (pcrc < 0 || pcrc > 2)
		env->ThrowError("gaussianblur:  pcrc must be in the range [0,2]!");
	if ((Y < -255 || Y > 3) || (U < -255 || U > 3) || (V < -255 || V > 3))
		env->ThrowError("gaussianblur:  Y, U, and V must be in the range [-255,3]!");
	if (Y==3)
		gbY = new GausBlur(vi.width,vi.height,sdy,border,integrate,unsharpstrength,
			gfunc,pcr,nthreads,env);
	if ((vi.IsYV12() || vi.IsYUY2()) && (U==3 || V==3))
		gbC = new GausBlur(vi.width/2,vi.IsYV12()?(vi.height/2):vi.height,sdc,border,
			integrate,unsharpstrength,gfuncc,pcrc,nthreads,env);
	if (vi.IsRGB() || vi.IsYUY2())
	{
		tempA = (unsigned char*)_aligned_malloc(mod16(vi.width)*vi.height*sizeof(unsigned char),16);
		tempB = (unsigned char*)_aligned_malloc(mod16(vi.width)*vi.height*sizeof(unsigned char),16);
	}
}

GaussianBlur::~GaussianBlur()
{
	delete gbY;
	delete gbC;
	_aligned_free(tempA);
	_aligned_free(tempB);
}

int extractYUY2(const unsigned char *s, const int pitch, const int width,
	const int height, const int b, unsigned char *d, const int nthreads)
{
	if (b == 0)
	{
		#pragma omp parallel for num_threads(nthreads)
		for (int y=0; y<height; ++y)
			for (int x=0; x<width; ++x)
				d[y*width+x] = s[y*pitch+x*2];
		return width;
	}
	else
	{
		const int widthd2 = width>>1;
		const int boff = 2*b-1;
		#pragma omp parallel for num_threads(nthreads)
		for (int y=0; y<height; ++y)
			for (int x=0; x<widthd2; ++x)
				d[y*widthd2+x] = s[y*pitch+x*4+boff];
		return widthd2;
	}
}

void mergeYUY2(unsigned char *d, const int pitch, const int width,
	const int height, const int b, const unsigned char *s, const int nthreads)
{
	if (b == 0)
	{
		#pragma omp parallel for num_threads(nthreads)
		for (int y=0; y<height; ++y)
			for (int x=0; x<width; ++x)
				d[y*pitch+x*2] = s[y*width+x];
	}
	else
	{
		const int widthd2 = width>>1;
		const int boff = 2*b-1;
		#pragma omp parallel for num_threads(nthreads)
		for (int y=0; y<height; ++y)
			for (int x=0; x<widthd2; ++x)
				d[y*pitch+x*4+boff] = s[y*widthd2+x];
	}
}

void extractRGB(const unsigned char *s, const int pitch, const int width,
	const int height, const int b, unsigned char *d, const int np,
	const int nthreads)
{
	#pragma omp parallel for num_threads(nthreads)
	for (int y=0; y<height; ++y)
		for (int x=0; x<width; ++x)
			d[y*width+x] = s[y*pitch+x*np+b];
}

void mergeRGB(unsigned char *d, const int pitch, const int width,
	const int height, const int b, const unsigned char *s, const int np,
	const int nthreads)
{
	#pragma omp parallel for num_threads(nthreads)
	for (int y=0; y<height; ++y)
		for (int x=0; x<width; ++x)
			d[y*pitch+x*np+b] = s[y*width+x];
}

PVideoFrame GaussianBlur::GetFrame(int n, IScriptEnvironment* env){

	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst = env->NewVideoFrame(vi);
	
	if (vi.IsYV12())
	{
		const int plane[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
		for (int b=0; b<3; ++b)
		{
			const unsigned char *srcp = src->GetReadPtr(plane[b]);
			const int spitch = src->GetPitch(plane[b]);
			unsigned char *dstp = dst->GetWritePtr(plane[b]);
			const int dpitch = dst->GetPitch(plane[b]);
			const int height = dst->GetHeight(plane[b]);
			const int width = dst->GetRowSize(plane[b]);
			const int val = b == 0 ? Y : b == 1 ? U : V;
			if (val == 3)
			{
				if (b == 0) gbY->Blur(srcp,spitch,dstp,dpitch,env);
				else gbC->Blur(srcp,spitch,dstp,dpitch,env);
			}
			else if (val == 2)
				env->BitBlt(dstp,dpitch,srcp,spitch,width,height);
			else if (val<1 && val>-256)
				memset(dstp,-val,dpitch*height);
		}
	}
	else
	{
		const unsigned char *srcp = src->GetReadPtr();
		const int spitch = src->GetPitch();
		unsigned char *dstp = dst->GetWritePtr();
		const int dpitch = dst->GetPitch();
		if (vi.IsYUY2())
		{
			for (int b=0; b<3; ++b)
			{
				const int val = b == 0 ? Y : b == 1 ? U : V;
				if (val == 3)
				{
					const int w = extractYUY2(srcp,spitch,vi.width,vi.height,b,tempA,nthreads);
					gbY->Blur(tempA,w,tempB,w,env);
					mergeYUY2(dstp,dpitch,vi.width,vi.height,b,tempB,nthreads);
				}
				else if (val == 2)
				{
					extractYUY2(srcp,spitch,vi.width,vi.height,b,tempA,nthreads);
					mergeYUY2(dstp,dpitch,vi.width,vi.height,b,tempA,nthreads);
				}
				else if (val<1 && val>-256)
				{
					memset(tempA,-val,vi.height*vi.width);
					mergeYUY2(dstp,dpitch,vi.width,vi.height,b,tempA,nthreads);
				}
			}
		}
		else
		{
			if (Y == 3)
			{
				const int np = vi.IsRGB32() ? 4 : 3;
				for (int b=0; b<np; ++b)
				{
					extractRGB(srcp,spitch,vi.width,vi.height,b,tempA,np,nthreads);
					gbY->Blur(tempA,vi.width,tempB,vi.width,env);
					mergeRGB(dstp,dpitch,vi.width,vi.height,b,tempB,np,nthreads);
				}
			}
			else if (Y == 2)
				env->BitBlt(dstp,dpitch,srcp,spitch,dst->GetRowSize(),vi.height);
			else if (Y<1 && Y>-256)
				memset(dstp,-Y,dpitch*vi.height);
		}
	}

	return dst;
}

double igfunc(const int pi, const double *gfvals, const int pcr)
{
	double p;
	if (pcr == 0) p = pi/255.0;
	else if (pcr == 1) p = (pi-16.0)/219.0;
	else p = (pi-16)/224.0;
	p = min(max(p,0.0),1.0);
	if (p < gfvals[4])
		return p/gfvals[0];
	return pow((p+gfvals[1])/(1.0+gfvals[1]),1.0/gfvals[2]);
}

int fgfunc(const double pi, const double *gfvals, const int pcr)
{
	const double p = min(max(pi,0.0),1.0);
	double t;
	if (p < gfvals[3])
		t = gfvals[0]*p;
	else
		t = (1.0+gfvals[1])*pow(p,gfvals[2])-gfvals[1];
	if (pcr == 0)
		return (int)(t*255.0+0.5);
	if (pcr == 1)
		return (int)(t*219.0+16.5);
	return (int)(t*224.0+16.5);
}

GausBlur::GausBlur(int _width, int _height, double sd, int _border, bool integrate,
	double unsharpstrength, int _gfunc, int _pcr, int _nthreads, IScriptEnvironment* env):
	width(_width), height(_height), sd(sd), border(_border), gfunc(_gfunc), pcr(_pcr),
	nthreads(_nthreads)
{
	in = out = realout = gaus = NULL;

	radius = (int)(sd*4.0+0.5);
	if (radius >= height)
		radius = height-1;
	if (radius >= width)
		radius = width-1;

	srcheight = height;
	srcwidth = width;
	if(border>2)
	{
		width+=radius*2;
		width=((width+31)>>5)*32;
		height+=radius*2;
		height=((height+31)>>5)*32;
	}
	else if(border<1)
	{
		width+=radius;
		width=((width+31)>>5)*32;
		height+=radius;
		height=((height+31)>>5)*32;
	}
	else if(border==1)
		radius=min(width,height)/2;
	outwidth = (width>>1)+1;
	norm = 1.0/(width*height);

	in = (float*)_aligned_malloc(height*width*sizeof(float),16);
	realout = (float*)_aligned_malloc(height*width*sizeof(float),16);
	out = (float*)_aligned_malloc(height*outwidth*2*sizeof(float),16);
	gaus = (float*)_aligned_malloc(height*outwidth*2*sizeof(float),16);

	hLib = LoadLibrary("libfftw3f-3.dll");
	if (hLib)
	{
		fftwf_plan_dft_r2c_2d = (fftwf_plan_dft_r2c_2d_proc)GetProcAddress(hLib,"fftwf_plan_dft_r2c_2d");
		fftwf_plan_dft_c2r_2d = (fftwf_plan_dft_c2r_2d_proc)GetProcAddress(hLib,"fftwf_plan_dft_c2r_2d");
		fftwf_execute = (fftwf_execute_plan_proc)GetProcAddress(hLib,"fftwf_execute");
		fftwf_destroy = (fftwf_destroy_plan_proc)GetProcAddress(hLib,"fftwf_destroy_plan");
		fftwf_init_threads = (fftwf_init_threads_proc)GetProcAddress(hLib,"fftwf_init_threads");
		fftwf_plan_nthreads = (fftwf_plan_nthreads_proc)GetProcAddress(hLib,"fftwf_plan_with_nthreads");
	}
	if (!hLib || !fftwf_plan_dft_r2c_2d || !fftwf_plan_dft_c2r_2d || !fftwf_execute || 
		!fftwf_destroy || !fftwf_init_threads || !fftwf_plan_nthreads)
		env->ThrowError("gaussianblur:  unable to load libfftw3f-3.dll!");

	fftwf_init_threads();
	fftwf_plan_nthreads(nthreads);
	plan = fftwf_plan_dft_r2c_2d(height,width,in,(fftw_complex*)out,FFTW_MEASURE|FFTW_DESTROY_INPUT); // direct fft 
	plani = fftwf_plan_dft_c2r_2d(height,width,(fftw_complex*)out,realout,FFTW_MEASURE|FFTW_DESTROY_INPUT); // inverse fft

	CalcGaus(integrate,unsharpstrength,env);

	if(border<1)
	{
		for (int i=0; i<width*height; i++)
			in[i] = abs(border);
	}

	if (gfunc >= 0)
	{
		for (int i=0; i<256; ++i)
			g_table[i] = igfunc(i,gamma_table[gfunc],pcr);
	}
}

GausBlur::~GausBlur()
{
	_aligned_free(in);
	_aligned_free(out);
	_aligned_free(realout);
	_aligned_free(gaus);
	fftwf_destroy(plan);
	fftwf_destroy(plani);
	FreeLibrary(hLib);
}

void GausBlur::CalcGaus(const bool integrate, const double strength, IScriptEnvironment *env)
{
	memset(in,0,height*width*sizeof(float));
	
	double totalsum = 0.0;

	for (int y=0; y<radius; y++)
	{
		if (integrate)
		{
			const double yint = RombergIntegral(y-0.5,y+0.5,sd);
			for (int x=0; x<radius; x++)
			{
				in[y*width+x] = RombergIntegral(x-0.5,x+0.5,sd)*yint;
				totalsum += in[y*width+x];
			}
		}
		else
		{
			const int ysqrd = y*y;
			for (int x=0; x<radius; x++)
			{
				in[y*width+x] = (1.0/(sd*sqrt(2.0*PI)))*exp(-(x*x+ysqrd)/(2.0*sd*sd));
				totalsum += in[y*width+x];
			}
		}
		for (int x=width-1; x>width-radius; x--)
		{
			in[y*width+x] = in[y*width+width-x];
			totalsum += in[y*width+x];
		}
	}
	for (int y=height-1; y>height-radius; y--)
	{
		for (int x=0; x<radius; x++)
		{
			in[y*width+x] = in[(height-y)*width+x];
			totalsum += in[y*width+x];
		}
		for (int x=width-1; x>width-radius; x--)
		{
			in[y*width+x] = in[(height-y)*width+width-x];
			totalsum += in[y*width+x];
		}
	}

	double factor = 1.0/totalsum;
	if (strength != 0.0)
		factor *= -strength;

	for (int y=0; y<radius; y++)
	{
		for (int x=0; x<radius; x++)
			in[y*width+x] *= factor;	
		for (int x=width-1; x>width-radius; x--)
			in[y*width+x] *= factor;	
	}
	for (int y=height-1; y>height-radius; y--)
	{
		for (int x=0; x<radius; x++)
			in[y*width+x] *= factor;
		for (int x=width-1; x>width-radius; x--)
			in[y*width+x] *= factor;	
	}

	if(strength != 0.0)
		in[0] += 1.0+strength;

	fftwf_execute(plan);

	env->BitBlt((BYTE*)gaus,outwidth*2*sizeof(float),(BYTE*)out,outwidth*2*sizeof(float),
		outwidth*2*sizeof(float),height);
}


 /***********************************************************
  * Integral of a function gaus(X) by Romberg's method		*
  * ------------------------------------------------------- *
  * INPUTS:													*
  *          a      begin value of x variable				*
  *          b      end value of x variable					*
  *			sd      standard deviation						*
  *															*
  *															*
  * RETURNED VALUE  the integral of gausian from a to b		*
  *															*
  * REFERENCE: "Mathematiques en Turbo-Pascal (Part 1) by	*
  *             Marc Ducamp and Alain Reverchon, Eyrolles,	*
  *             Paris, 1987"								*
  *                             C++ version by J-P Moreau   *
  *            http://perso.wanadoo.fr/jean-pierre.moreau/  *
  *								    modified by t petersen  *
  **********************************************************/

double rfunc(const double x, const double sd)
{
	return (1.0/(sd*sqrt(2.0*PI)))*exp(-(x*x)/(2.0*sd*sd));
}

double GausBlur::RombergIntegral(double a, double b, double sd) 
{
	const double PREC=0.00001;
	const int ITERMIN=1;	
	const int ITERMAX=100;

	double t[ITERMAX][ITERMAX];

	int n = 0;
	double obtprec;
	double ta = (rfunc(a,sd) + rfunc(b,sd))*0.5;
	double pas = b-a;
	t[0][0] = ta*pas;
	do
	{
		n = n + 1;
		pas = pas*0.5;
		double s = ta;
		for (int i=1; i<(1<<n); i++)
			s += rfunc(a+pas*i,sd);
		t[0][n] = s*pas;
		double r = 1.0;
		for (int i=1; i<n+1; i++) 
		{
			r = r*4.0;
			int j = n-i;
			t[i][j]=(r*t[i-1][j+1]-t[i-1][j])/(r-1.0);
		}
		obtprec = fabs(t[n][0]-t[n-1][0]);
		if (n>=ITERMAX)
			return t[n][0];
	} 
	while (obtprec>PREC || n<ITERMIN); 
	return t[n][0];
}


 __declspec(align(16)) const float sse_0p5[4] = { 0.5f, 0.5f, 0.5f, 0.5f };
 __declspec(align(16)) const float sse_255[4] = { 255.0f, 255.0f, 255.0f, 255.0f };
 __declspec(align(16)) const float cmask[4] = { -1.0f, 1.0f, -1.0f, 1.0f };

void uc2float(const unsigned char *s, float *d, int width)
{
	const int remain = width&7;
	width -= remain;
	if (width <= 0)
		goto extra;
	__asm
	{
		mov rax,s
		mov rdx,d
		mov esi,width
		xor rcx,rcx
		pxor xmm2,xmm2
loop8u:
		movq xmm0,QWORD PTR[rax+rcx]
		punpcklbw xmm0,xmm2
		movdqa xmm1,xmm0
		punpcklbw xmm0,xmm2
		punpckhbw xmm1,xmm2
		cvtdq2ps xmm0,xmm0
		cvtdq2ps xmm1,xmm1
		movups [rdx+rcx*4],xmm0
		movups [rdx+rcx*4+16],xmm1
		add rcx,8
		cmp rcx,rsi
		jl loop8u
	}
extra:
	for (int x=width; x<width+remain; ++x)
		d[x] = s[x];
}

void float2uc(const float *s, unsigned char *d, int width, const float *norm)
{
	const int remain = width&7;
	width -= remain;
	if (width <= 0)
		goto extra;
	__asm
	{
		mov rcx,norm
		mov rax,s
		movss xmm2,[rcx]
		mov rdx,d
		mov esi,width
		xor rcx,rcx
		shufps xmm2,xmm2,0
		movaps xmm3,sse_0p5
		movaps xmm4,sse_255
		xorps xmm5,xmm5
loop8u:
		movups xmm0,[rax+rcx*4]
		movups xmm1,[rax+rcx*4+16]
		mulps xmm0,xmm2
		mulps xmm1,xmm2
		addps xmm0,xmm3
		addps xmm1,xmm3
		minps xmm0,xmm4
		minps xmm1,xmm4
		maxps xmm0,xmm5
		maxps xmm1,xmm5
		cvttps2dq xmm0,xmm0
		cvttps2dq xmm1,xmm1
		packssdw xmm0,xmm1
		packuswb xmm0,xmm0
		movq QWORD PTR[rdx+rcx],xmm0
		add rcx,8
		cmp rcx,rsi
		jl loop8u
	}
extra:
	for (int x=width; x<width+remain; ++x)
		d[x] = max(min((int)(s[x]*norm[0]+0.5f),255),0);
}

void lineMult(float *s, const float *p, int width)
{
	const int remain = width&7;
	width -= remain;
	if (width <= 0)
		goto extra;
	__asm
	{
		mov rax,s
		mov rdx,p
		mov esi,width
		xor rcx,rcx
loop8u:
		movups xmm0,[rax+rcx*4]
		movups xmm2,[rdx+rcx*4]
		movups xmm4,[rax+rcx*4+16]
		movups xmm6,[rdx+rcx*4+16]
		movaps xmm1,xmm0
		movaps xmm3,xmm2
		movaps xmm5,xmm4
		movaps xmm7,xmm6
		shufps xmm0,xmm0,(2<<6)+(2<<4)+(0<<2)+(0<<0)
		shufps xmm1,xmm1,(3<<6)+(3<<4)+(1<<2)+(1<<0)
		shufps xmm3,xmm3,(2<<6)+(3<<4)+(0<<2)+(1<<0)
		shufps xmm4,xmm4,(2<<6)+(2<<4)+(0<<2)+(0<<0)
		shufps xmm5,xmm5,(3<<6)+(3<<4)+(1<<2)+(1<<0)
		shufps xmm7,xmm7,(2<<6)+(3<<4)+(0<<2)+(1<<0)
		mulps xmm0,xmm2
		mulps xmm1,xmm3
		mulps xmm4,xmm6
		mulps xmm5,xmm7
		mulps xmm1,cmask
		mulps xmm5,cmask
		addps xmm0,xmm1
		addps xmm4,xmm5
		movups [rax+rcx*4],xmm0
		movups [rax+rcx*4+16],xmm4
		add rcx,8
		cmp rcx,rsi
		jl loop8u
	}
extra:
	for (int x=width; x<width+remain; x+=2)
	{
		s[x+0] = s[x+0]*p[x+0]-s[x+1]*p[x+1];
		s[x+1] = s[x+0]*p[x+1]+s[x+1]*p[x+0];
	}
}

void GausBlur::Blur(const unsigned char *srcp, const int spitch, unsigned char *dstp, 
		const int dpitch, IScriptEnvironment* env)
{
	if (gfunc >= 0)
	{
		const int doff = (border==1||border==2)?0:radius*width+radius;
		#pragma omp parallel for num_threads(nthreads)
		for(int y=0;y<srcheight;y++)
			for(int x=0;x<srcwidth;x++)
				in[y*width+x+doff]=g_table[srcp[y*spitch+x]];
	}
	else
	{
		const int doff = (border==1||border==2)?0:radius*width+radius;
		#pragma omp parallel for num_threads(nthreads)
		for(int y=0;y<srcheight;y++)
			uc2float(srcp+y*spitch,in+y*width+doff,srcwidth);
	}

	if(border>=3)
	{
		const int inc = border > 3 ? 1 : 0;
		//LeftBorder and right border
		#pragma omp parallel for num_threads(nthreads)
		for(int y=0;y<srcheight;y++)
		{
			const unsigned char *srcpT = srcp+y*spitch;
			float *inT = in+(y+radius)*width;
			if (gfunc >= 0)
			{
				for(int x=0;x<radius;x++)
					inT[radius-1-x]=g_table[srcpT[(x+1)*inc]];
				for(int x=0;x<radius;x++)
					inT[radius+srcwidth+x]=g_table[srcpT[srcwidth-1-(x+1)*inc]];
			}
			else
			{
				for(int x=0;x<radius;x++)
					inT[radius-1-x]=srcpT[(x+1)*inc];
				for(int x=0;x<radius;x++)
					inT[radius+srcwidth+x]=srcpT[srcwidth-1-(x+1)*inc];
			}
		}
		//TopBorder;
		int srcoffset=width*(radius+inc);
		int dstoffset=width*(radius-1);
		for(int y=0;y<radius;y++,dstoffset-=width,srcoffset+=inc*width)
			memcpy(in+dstoffset,in+srcoffset,width*sizeof(float));
		//buttomborder;
		srcoffset=width*(radius+srcheight-1-inc);
		dstoffset=width*(radius+srcheight);
		for(int y=0;y<radius;y++,dstoffset+=width,srcoffset-=inc*width)
			memcpy(in+dstoffset,in+srcoffset,width*sizeof(float));
	}

	fftwf_execute(plan);

	#pragma omp parallel for num_threads(nthreads)	
	for(int y=0;y<height;y++)
		lineMult(out+y*outwidth*2,gaus+y*outwidth*2,outwidth*2);

	fftwf_execute(plani);

	if(border==2)
	{
		env->BitBlt(dstp,dpitch,srcp,spitch,srcwidth,srcheight);
		if (gfunc >= 0)
		{
			#pragma omp parallel for num_threads(nthreads)
			for(int y=radius;y<srcheight-radius;y++)
				for(int x=radius;x<srcwidth-radius;x++)
					dstp[y*dpitch+x]=fgfunc(realout[y*width+x]*norm,gamma_table[gfunc],pcr);
		}
		else
		{
			#pragma omp parallel for num_threads(nthreads)
			for(int y=radius;y<srcheight-radius;y++)
				float2uc(realout+y*width+radius,dstp+y*dpitch+radius,srcwidth-2*radius,&norm);
		}
	}
	else
	{
		const int soff = border==1?0:radius*width+radius;
		if (gfunc >= 0)
		{
			#pragma omp parallel for num_threads(nthreads)
			for(int y=0;y<srcheight;y++)
				for(int x=0;x<srcwidth;x++)
					dstp[y*dpitch+x]=fgfunc(realout[y*width+x+soff]*norm,gamma_table[gfunc],pcr);
		}
		else
		{
			#pragma omp parallel for num_threads(nthreads)
			for(int y=0;y<srcheight;y++)
				float2uc(realout+y*width+soff,dstp+y*dpitch,srcwidth,&norm);
		}
	}	
}

