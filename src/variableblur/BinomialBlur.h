
#include <windows.h>
#include <math.h>
#include "avisynth.h"

class BinomialBlur :public GenericVideoFilter {

public:
	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
	BinomialBlur(PClip clip,int rady,int radc,int y,int u,int v,bool gaus,bool usemmx,bool notusedefault,bool n5,IScriptEnvironment* env);
	~BinomialBlur();
protected:
	unsigned int radY,radC;
	int Y,U,V;
	bool useMMX;
	bool defmmx;
	bool gaussian;
	bool useN5;
	unsigned short* buffer;
	int bufferSize;
	void CalcFrame(int plane,PVideoFrame &src,PVideoFrame &dst,int rad);

};

void CalcHorizontal(const unsigned short* sp,unsigned char* dp,int width,int height,int spitch,int dpitch,int rad);
void CalcVertical(const unsigned char* sp,unsigned short* dp,int width,int height,int spitch,int dpitch,int rad);


void CalcGaussianN3MMX(int width,int height,unsigned short* srcp,unsigned short* SC0,unsigned short* SC1);
void CalcGaussianN5MMX(int width,int height,unsigned short* srcp,unsigned short* SC01,unsigned short* SC23);
void CalcGaussianN5(const unsigned char* sp,int spitch,int width,int height,unsigned char* dstp,int dst_pitch,int* SC0,int* SC1,int* SC2,int* SC3);
void CalcGaussianN3(const unsigned char* sp,int spitch,int width,int height,unsigned char* dstp,int dst_pitch,int* SC0,int* SC1);

int PackFrame(int width,int height,const unsigned char* srcp,unsigned short* &inputf,int pitch);
void unpackframe(int width,int height,unsigned short* &inputf,unsigned char* dstp,int spitch,int dpitch);