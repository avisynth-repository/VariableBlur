

#include "BinomialBlur.h"

BinomialBlur::BinomialBlur(PClip clip,int rady,int radc,int y,int u,int v,bool gaus,bool usemmx,bool notusedefault,bool n5,IScriptEnvironment* env):
GenericVideoFilter(clip),radY(rady),radC(radc),Y(y),U(u),V(v),useMMX(usemmx),buffer(0),bufferSize(0),gaussian(gaus),defmmx(notusedefault),useN5(n5)
{}

BinomialBlur::~BinomialBlur(){
	delete[] buffer;
}

PVideoFrame BinomialBlur::GetFrame(int n, IScriptEnvironment* env){

	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst = env->NewVideoFrame(vi);
		
	int temp=src->GetHeight(PLANAR_Y)*src->GetRowSize(PLANAR_Y)>>(Y==3?0:1);
	if(bufferSize<temp){
		delete[] buffer;
		bufferSize=temp;
		buffer=new unsigned short[temp];
	}

	if(Y==3)
		CalcFrame(PLANAR_Y,src,dst,radY);
	else if(Y==2)
		env->BitBlt(dst->GetWritePtr(PLANAR_Y), dst->GetPitch(PLANAR_Y), src->GetReadPtr(PLANAR_Y), src->GetPitch(PLANAR_Y),src->GetRowSize(PLANAR_Y), src->GetHeight(PLANAR_Y));
	else if(Y<1&&Y>-256)
		memset (dst->GetWritePtr(PLANAR_Y), -Y, dst->GetPitch(PLANAR_Y)*dst->GetHeight(PLANAR_Y));

	if(U==3)
		CalcFrame(PLANAR_U,src,dst,radC);
	else if(U==2)
		env->BitBlt(dst->GetWritePtr(PLANAR_U), dst->GetPitch(PLANAR_U), src->GetReadPtr(PLANAR_U), src->GetPitch(PLANAR_U),src->GetRowSize(PLANAR_U), src->GetHeight(PLANAR_U));
	else if(U<1&&U>-256)
		memset (dst->GetWritePtr(PLANAR_U), -U, dst->GetPitch(PLANAR_U)*dst->GetHeight(PLANAR_U));

	if(V==3)
		CalcFrame(PLANAR_V,src,dst,radC);
	else if(V==2)
		env->BitBlt(dst->GetWritePtr(PLANAR_V), dst->GetPitch(PLANAR_V), src->GetReadPtr(PLANAR_V), src->GetPitch(PLANAR_V),src->GetRowSize(PLANAR_V), src->GetHeight(PLANAR_V));
	else if(V<1&&V>-256)
		memset (dst->GetWritePtr(PLANAR_V), -V, dst->GetPitch(PLANAR_V)*dst->GetHeight(PLANAR_V));
	return dst;
}




void BinomialBlur::CalcFrame(int plane,PVideoFrame &src,PVideoFrame &dst,int rad)
{
	int src_pitch=src->GetPitch(plane);
	int dst_pitch=dst->GetPitch(plane);


	int width=src->GetRowSize(plane);
	int height=src->GetHeight(plane);

	unsigned char* dstp=dst->GetWritePtr(plane);
	const unsigned char* srcp=src->GetReadPtr(plane);

	if(gaussian){

		if(useMMX&&(defmmx||rad>2))
		{
		unsigned short *SC0=new _declspec(align(8)) unsigned short[(width+8)*2];
		unsigned short *SC1=new _declspec(align(8)) unsigned short[(width+8)*2];


		unsigned short* packedframe=0;		
		int pfwidth=PackFrame(width,height,srcp,packedframe,src_pitch);
		if(useN5){
		for(int N=0;N<(rad>>1);N++)
			CalcGaussianN5MMX(pfwidth<<1,height,packedframe,SC0,SC1);
		if(rad&1)
			CalcGaussianN3MMX(pfwidth<<1,height,packedframe,SC0,SC1);
		delete[] SC0;
		delete[] SC1;
		}
else
		for(int N=0;N<rad;N++)
			CalcGaussianN3MMX(pfwidth<<1,height,packedframe,SC0,SC1);

		unpackframe(pfwidth,height,packedframe,dstp,pfwidth,dst_pitch);

		}
		else
		{

			if(useN5)
			{				
				_declspec(align(8)) int *SC0=new int[width+2];
				_declspec(align(8)) int *SC1=new int[width+2];
				_declspec(align(8)) int *SC2=new int[width+2];
				_declspec(align(8)) int *SC3=new int[width+2];
				const unsigned char* sp=srcp;
				int spitch=src_pitch;

				for(int N=0;N<(rad>>1);N++){
					CalcGaussianN5(sp,spitch,width,height,dstp,dst_pitch,SC0,SC1,SC2,SC3);
					sp=dstp;
					spitch=dst_pitch;
				}
				if(rad&1)
					CalcGaussianN3(sp,spitch,width,height,dstp,dst_pitch,SC0,SC1);
				delete[] SC0;
				delete[] SC1;
				delete[] SC2;
				delete[] SC3;
			}
			else{

				_declspec(align(8)) int *SC0=new int[width+1];
				_declspec(align(8)) int *SC1=new int[width+1];
				const unsigned char* sp=srcp;
				int spitch=src_pitch;
				for(int N=0;N<rad;N++){
					CalcGaussianN3(sp,spitch,width,height,dstp,dst_pitch,SC0,SC1);

					sp=dstp;
					spitch=dst_pitch;
				}
				delete[] SC0;
				delete[] SC1;
			}
		}
	}
	else{
		if(!useN5){
			//calc vertical;
			for(int y=0,soffset=0,doffset=0;y<height;y++,soffset+=src_pitch,doffset+=width){
				//left border
				for(int x=0,Val=0;x<rad;x++,Val=0){
					Val=srcp[soffset]*(rad-x);
					for(int i=-x;i<=rad;i++)
						Val+=srcp[soffset+x+i];
					buffer[doffset+x]=Val;
				}
				//right border
				for(int x=width-rad,Val=0;x<width;x++,Val=0){
					Val=srcp[soffset+width-1]*(x-width+rad+1);
					for(int i=-rad;i<=width-x;i++)
						Val+=srcp[soffset+x+i];
					buffer[doffset+x]=Val;
				}

				for(int x=rad,Val=0;x<(width-rad);x++,Val=0){
					for(int i=-rad;i<=rad;i++)
						Val+=srcp[soffset+x+i];
					buffer[doffset+x]=Val;
				}
			}


			//calc horizontal;
			int boxsize=(rad*2+1)*(rad*2+1);
			//Calc Upper border;
			for(int y=0,soffset=0,doffset=0;y<rad;y++,soffset+=width,doffset+=dst_pitch)
				for(int x=0,Val=0;x<width;x++,Val=0){
					Val=buffer[x]*(rad-y);
					for(int i=-y*width;i<=rad*width;i+=width)
						Val+=buffer[soffset+x+i];
					dstp[doffset+x]=Val/boxsize;
				}
				//Calc lower border;
				for(int y=height-rad,soffset=(height-rad)*width,doffset=(height-rad)*dst_pitch;y<height;y++,soffset+=width,doffset+=dst_pitch)
					for(int x=0,Val=0;x<width;x++,Val=0){
						Val=buffer[(height-1)*width+x]*(y-height+rad+1);
						for(int i=-rad*width;i<(height-y)*width;i+=width)
							Val+=buffer[soffset+x+i];
						dstp[doffset+x]=Val/boxsize;
					}	


					for(int y=rad,soffset=rad*width,doffset=rad*dst_pitch;y<(height-rad);y++,soffset+=width,doffset+=dst_pitch)
						for(int x=0,Val=0;x<width;x++,Val=0){
							for(int i=-rad*width;i<=rad*width;i+=width)
								Val+=buffer[soffset+x+i];
							dstp[doffset+x]=Val/boxsize;
						}
		
	}
		else
	{
		CalcVertical(srcp,buffer,width,height,src_pitch,width,rad);
		CalcHorizontal(buffer,dstp,width,height,width,dst_pitch,rad);
	}
	}

}



void CalcVertical(const unsigned char* sp,unsigned short* dp,int width,int height,int spitch,int dpitch,int rad)
{
	int div=rad*2+1;
	int val,front;
	for(int y=0,doffset=0,soffset=0;y<height;y++,soffset+=spitch,doffset+=dpitch){
		front=sp[soffset];
		val=front*(rad+1);
		for(int x=0;x<rad;x++)
			val+=sp[soffset+x];
		for(int x=rad;x<div;x++){
			val+=sp[soffset+x]-front;
			dp[doffset-rad+x]=val;}
		for(int x=div;x<width;x++){
			val+=sp[soffset+x]-sp[soffset+x-div];
			dp[doffset-rad+x]=val;}
		for(int x=width,back=sp[soffset+width-1];x<width+rad;x++){
			val+=back-sp[soffset+x-div];
			dp[doffset-rad+x]=val;}
	}
}

void CalcHorizontal(const unsigned short* sp,unsigned char* dp,int width,int height,int spitch,int dpitch,int rad)
{
	int div=rad*2+1;
	int divisor=div*div;
	int divoff=div*spitch;
	int val,front;
	for(int x=0,doffset=0,soffset=0;x<width;++x,soffset=x,doffset=x){
		front=sp[x];
		val=front*(rad+1);
		for(int y=0;y<rad;y++,soffset+=spitch)
			val+=sp[soffset];
		for(int y=rad;y<div;y++,soffset+=spitch,doffset+=dpitch){
			val+=sp[soffset]-front;
			dp[doffset]=val/divisor;}
		for(int y=div;y<height;y++,soffset+=spitch,doffset+=dpitch){
			val+=sp[soffset]-sp[soffset-divoff];
			dp[doffset]=val/divisor;}
		for(int y=height,back=sp[soffset-spitch];y<height+rad;y++,soffset+=spitch,doffset+=dpitch){
			val+=back-sp[soffset-divoff];
			dp[doffset]=val/divisor;}
	}
}
int PackFrame(int width,int height,const unsigned char* srcp,unsigned short* &inputf,int pitch)
{
	
	//int xsize = ((4 - (width & 3) & 3) + width);//MOD4 
	int xsize = (width+3)&0xFFFFFFFC;
	//bool notmod4=xsize!=width;
	int misalign=xsize-width;
	int insizex =xsize;
	inputf=new _declspec(align(8)) unsigned short[insizex*height];

	int i,j,k;
	for (i = 0; i < height; i++)
		{
			for (k = 0; k < 4; k++)
				for (j = 0; j < xsize / 4; j++)
				{
					inputf[i * insizex + j * 4 + k] = srcp[i * pitch + j + k * xsize / 4];
				}
			if(misalign)
				for (j = xsize / 4-misalign; j < xsize / 4; j++)
					inputf[i * insizex + j * 4 + 3] = srcp[i*pitch+width];
		}
	return insizex;
}

void unpackframe(int width,int height,unsigned short* &inputf,unsigned char* dstp,int spitch,int dpitch)
{
		int i,j,k;
		for (i = 0; i < height; i++)
		{
			for (k = 0; k < 4; k++)
				for (j = 0; j < width / 4; j++)
				{
					dstp[i * dpitch + j + k * width / 4] = inputf[i * spitch + j * 4 + k];
				}

		}
		delete[] inputf;
}

void CalcGaussianN3(const unsigned char* sp,int spitch,int width,int height,unsigned char* dstp,int dst_pitch,int* SC0,int* SC1)
{

	int SR0=sp[0];
	int SR1=SR0<<1;
	int temp1,temp2;
	for(int x=1;x<width;x++){
		temp1=sp[x];
		temp2=SR0+temp1;
		SR0=temp1;
		temp1=SR1+temp2;
		SR1=temp2;

		SC1[x]=temp1<<1;
		SC0[x]=temp1;

	}
	temp2=SR0<<1;
	temp1=SR1+temp2;
	temp2=temp1<<1;
	SC0[width]=temp1;
	SC1[width]=temp2;
	int doffset=0;
	for(int y=1,soffset=spitch;y<height;y++,soffset+=spitch,doffset+=dst_pitch){
		SR0=sp[soffset];
		SR1=SR0<<1;
		for(int x=1;x<width;x++){
			temp1=sp[x+soffset];
			temp2=SR0+temp1;
			SR0=temp1;
			temp1=SR1+temp2;
			SR1=temp2;

			temp2=SC0[x]+temp1;
			SC0[x]=temp1;
			dstp[doffset+x-1]=(8+SC1[x]+temp2)>>4;
			SC1[x]=temp2;
		}
		temp2=SR0<<1;
		temp1=SR1+temp2;
		temp2=SC0[width]+temp1;
		SC0[width]=temp1;
		dstp[doffset+width-1]=(8+SC1[width]+temp2)>>4;
		SC1[width]=temp2;
	}
	for(int x=1;x<=width;x++){
		temp2=SC0[x]<<1;
		dstp[doffset+x-1]=(8+SC1[x]+temp2)>>4;

	}
}

void CalcGaussianN5(const unsigned char* sp,int spitch,int width,int height,unsigned char* dstp,int dst_pitch,int* SC0,int* SC1,int* SC2,int* SC3)
{
	//init SC0,1,2,3
	int temp1,temp2;
	int SR0=sp[0];
	int SR1=SR0<<1;
	int SR2=SR1<<1;
	int SR3=SR2<<1;
	//fist line
	for(int x=0;x<width;x++){
		temp1=sp[x];
		temp2=SR0+temp1;
		SR0=temp1;
		temp1=SR1+temp2;
		SR1=temp2;
		temp2=SR2+temp1;
		SR2=temp1;
		temp1=SR3+temp2;
		SR3=temp2;

		temp2=temp1<<1;
		SC0[x]=temp1;
		temp1=temp2<<1;
		SC1[x]=temp2;
		temp2=temp1<<1;
		SC2[x]=temp1;
		SC3[x]=temp2;



	}
	//second last pixel
	temp2=SR0<<1;
	temp1=SR1+temp2;
	SR1=temp2;
	temp2=SR2+temp1;
	SR2=temp1;
	temp1=SR3+temp2;
	SR3=temp2;

	SC0[width]=temp1;
	temp2=temp1<<1;
	SC1[width]=temp2;
	temp1=temp2<<1;
	SC2[width]=temp1;
	temp2=temp1<<1;
	SC3[width]=temp2;

	//last pixel
	temp2=SR0<<1;
	temp1=SR1+temp2;
	temp2=SR2+temp1;
	temp1=SR3+temp2;

	SC0[width+1]=temp1;
	temp2=temp1<<1;
	SC1[width+1]=temp2;
	temp1=temp2<<1;
	SC2[width+1]=temp1;
	temp2=temp1<<1;
	SC3[width+1]=temp2;
	//second line
	{
		int soffset=spitch;
		SR0=sp[soffset];
		SR1=SR0<<1;
		SR2=SR1<<1;
		SR3=SR2<<1;
		for(int x=0;x<2;x++){
			temp1=sp[x+soffset];
			temp2=SR0+temp1;
			SR0=temp1;
			temp1=SR1+temp2;
			SR1=temp2;
			temp2=SR2+temp1;
			SR2=temp1;
			temp1=SR3+temp2;
			SR3=temp2;

			temp2=SC0[x]+temp1;
			SC0[x]=temp1;
			temp1=SC1[x]+temp2;
			SC1[x]=temp2;
			temp2=SC2[x]+temp1;
			SC2[x]=temp1;

			SC3[x]=temp2;

		}

		for(int x=2;x<width;x++){
			temp1=sp[x+soffset];
			temp2=SR0+temp1;
			SR0=temp1;
			temp1=SR1+temp2;
			SR1=temp2;
			temp2=SR2+temp1;
			SR2=temp1;
			temp1=SR3+temp2;
			SR3=temp2;

			temp2=SC0[x]+temp1;
			SC0[x]=temp1;
			temp1=SC1[x]+temp2;
			SC1[x]=temp2;
			temp2=SC2[x]+temp1;
			SC2[x]=temp1;
			//dstp[doffset+x-2]=(128+SC3[x]+temp2)>>8;
			SC3[x]=temp2;
		}

		//second last pixel
		temp2=SR0<<1;
		temp1=SR1+temp2;
		SR1=temp2;
		temp2=SR2+temp1;
		SR2=temp1;
		temp1=SR3+temp2;
		SR3=temp2;

		temp2=SC0[width]+temp1;
		SC0[width]=temp1;
		temp1=SC1[width]+temp2;
		SC1[width]=temp2;
		temp2=SC2[width]+temp1;
		SC2[width]=temp1;
		//dstp[doffset+width-2]=(128+SC3[x]+temp2)>>8;
		SC3[width]=temp2;

		//last pixel
		temp2=SR0<<1;
		temp1=SR1+temp2;
		temp2=SR2+temp1;
		temp1=SR3+temp2;

		temp2=SC0[width+1]+temp1;
		SC0[width+1]=temp1;
		temp1=SC1[width+1]+temp2;
		SC1[width+1]=temp2;
		temp2=SC2[width+1]+temp1;
		SC2[width+1]=temp1;
		SC3[width+1]=temp2;
	}
	//main loop
	int doffset=0;
	for(int y=2,soffset=spitch<<1;y<height;y++,soffset+=spitch,doffset+=dst_pitch){
		SR0=sp[soffset];
		SR1=SR0<<1;
		SR2=SR1<<1;
		SR3=SR2<<1;
		for(int x=0;x<2;x++){
			temp1=sp[x+soffset];
			temp2=SR0+temp1;
			SR0=temp1;
			temp1=SR1+temp2;
			SR1=temp2;
			temp2=SR2+temp1;
			SR2=temp1;
			temp1=SR3+temp2;
			SR3=temp2;

			temp2=SC0[x]+temp1;
			SC0[x]=temp1;
			temp1=SC1[x]+temp2;
			SC1[x]=temp2;
			temp2=SC2[x]+temp1;
			SC2[x]=temp1;

			SC3[x]=temp2;

		}

		for(int x=2;x<width;x++){
			temp1=sp[x+soffset];
			temp2=SR0+temp1;
			SR0=temp1;
			temp1=SR1+temp2;
			SR1=temp2;
			temp2=SR2+temp1;
			SR2=temp1;
			temp1=SR3+temp2;
			SR3=temp2;

			temp2=SC0[x]+temp1;
			SC0[x]=temp1;
			temp1=SC1[x]+temp2;
			SC1[x]=temp2;
			temp2=SC2[x]+temp1;
			SC2[x]=temp1;
			dstp[doffset+x-2]=(128+SC3[x]+temp2)>>8;
			SC3[x]=temp2;
		}

		//second last pixel
		temp2=SR0<<1;
		temp1=SR1+temp2;
		SR1=temp2;
		temp2=SR2+temp1;
		SR2=temp1;
		temp1=SR3+temp2;
		SR3=temp2;

		temp2=SC0[width]+temp1;
		SC0[width]=temp1;
		temp1=SC1[width]+temp2;
		SC1[width]=temp2;
		temp2=SC2[width]+temp1;
		SC2[width]=temp1;
		dstp[doffset+width-2]=(128+SC3[width]+temp2)>>8;
		SC3[width]=temp2;

		//last pixel
		temp2=SR0<<1;
		temp1=SR1+temp2;
		temp2=SR2+temp1;
		temp1=SR3+temp2;

		temp2=SC0[width+1]+temp1;
		SC0[width+1]=temp1;
		temp1=SC1[width+1]+temp2;
		SC1[width+1]=temp2;
		temp2=SC2[width+1]+temp1;
		SC2[width+1]=temp1;
		dstp[doffset+width-1]=(128+SC3[width+1]+temp2)>>8;
		SC3[width+1]=temp2;

	}
	//last line
	for(int x=2;x<=width+1;x++){
		temp2=SC0[x]<<1;
		temp1=SC1[x]+temp2;
		SC1[x]=temp2;
		temp2=SC2[x]+temp1;
		SC2[x]=temp1;
		SC3[x]=temp2;
		dstp[doffset+x-2]=(128+SC3[x]+temp2)>>8;

	}
	doffset+=dst_pitch;
	for(int x=2;x<=width+1;x++){
		temp2=SC0[x]<<1;
		temp1=SC1[x]+temp2;

		temp2=SC2[x]+temp1;
		dstp[doffset+x-2]=(128+SC3[x]+temp2)>>8;

	}
}










void CalcGaussianN3MMX(int width,int height,unsigned short* srcp,unsigned short* SC0,unsigned short* SC1)
{

const __int64 const8=0x0008000800080008;
//const __int64 mask=  0xFFFFFFFFFFFF0000;
const __int64 mask=  0x000000000000FFFF;

int endline=width*height;
//memset(SC1,0,width);
	__asm{
		
		movq mm7,const8		
		movq mm6,mask
		

		mov ebx,width
		mov rcx,srcp
		mov rax,rcx
		add endline,eax
		xor rdx,rdx	//rdx=0
		mov rdi,SC0
		mov rsi,SC1
		//rcx=soffset
		//rdi SC0;
		//rsi SC1;
		//rdx=x
		//rax=doffset
		
		
		//mm0 temp1
		//mm1 temp2
		//mm2 SR0
		//mm3 SR1
#define temp1 mm0
#define temp2 mm1
#define SR0 mm2
#define SR1 mm3
#define SC0(index) [rdi+(index)]
#define SC1(index) [rsi+(index)]
		//setup SR0 and SR1
		movq SR0,[rcx+rbx-2*4]	//Move last qword into SR0
		
		movq SR1,[rcx+rbx-2*8]	//Move the second last qword into SR1
		
		paddw SR1,SR0			//add SR0 to SR1 SR1:OK
		
		//setup SC0 and SC1
//first Qword
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq mm5,temp1		  //save first pixel
		psrlq temp1,16		//shift first pixel out(to calculate the last)
		
		movq mm4,SR0		//move SR0 in to extract last pixel
		psrlq mm4,48		
		psllq mm4,48		
		paddw temp1,mm4		//add last pixel

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0
		
		movq SR0,temp1		//SR0=temp1

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2

		movq SR1,temp2		//SR1=temp2

		movq SC0(rdx),temp1	//SC0[x]=temp1
		paddw temp1,temp1
		movq SC1(rdx),temp1		//SC1[x]=0

		
		
		pand mm5,mm6		//select first pixel
		psllq SR0,16		//align
		paddw SR0,mm5		//SR0:ok
		
		psllq SR1,16		//align
		paddw SR1,mm5
		paddw SR1,mm5		//SR1:ok

		add rdx,8			//increase counter
initloop:
		movq temp1,[rcx+rdx]  //temp1=srcp[x]

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0
		
		movq SR0,temp1		//SR0=temp1

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2

		movq SR1,temp2		//SR1=temp2

		movq SC0(rdx),temp1	//SC0[x]=temp1
		paddw temp1,temp1	//temp1=temp1*2
		movq SC1(rdx),temp1		//SC1[x]=temp1

		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae initloop		//continue loop if less than width

		
		//rdx=index for SC0/SC1;

yloop:
		
		add rcx,rbx			//soffset=next line
		cmp ecx,endline
		jae lastline
		//setup SR0 and SR1
		movq SR0,[rcx+rbx-2*4]	//Move last 4 short into SR0
		
		movq SR1,[rcx+rbx-2*8]	//Move the second last 4 shorts into SR1
	
		paddw SR1,SR0			//add SR0 to SR1 SR1:OK
		
		xor rdx,rdx			//rdx =x=0
		
//first Qword	
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq mm5,temp1		  //save first pixel
		psrlq temp1,16		  //shift first pixel out
		movq mm4,SR0		//move SR0 in to extract last pixel
		psrlq mm4,48		
		psllq mm4,48		
		paddw temp1,mm4		//add last pixel

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0
		
		movq SR0,temp1		//SR0=temp1

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2

		movq SR1,temp2		//SR1=temp2

		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]

		movq SC0(rdx),temp1	//SC0[x]=temp1

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2  //temp1=SC1[x]+temp2
		paddw temp1,mm7		//temp1=SC1[x]+temp2+8  (+8 because off rounding)
		psrlw temp1,4		//temp1=(SC1[x]+temp2+8)/16
		movq [rcx+rdx-8],temp1		//dst=temp1
		movq SC1(rdx),temp2		//SC1[x]=temp2
		
		
		pand mm5,mm6		//select first pixel
		psllq SR0,16		//align
		paddw SR0,mm5		//SR0:ok
		
		psllq SR1,16		//align
		paddw SR1,mm5		//SR1=SR0*2
		paddw SR1,mm5		//SR1:ok

		add rdx,8			//increase counter
		

xloop:
		movq temp1,[rcx+rdx]  //temp1=srcp[x]

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0
		
		movq SR0,temp1		//SR0=temp1

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2

		movq SR1,temp2		//SR1=temp2

		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]

		movq SC0(rdx),temp1	//SC0[x]=temp1

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2  //temp1=SC1[x]+temp2
		paddw temp1,mm7		//temp1=SC1[x]+temp2+8  (+8 because off rounding)
		psrlw temp1,4		//temp1=(SC1[x]+temp2+8)/16
		movq [rax+rdx-8],temp1		//dst=temp1
		movq SC1(rdx),temp2		//SC1[x]=temp2

		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae xloop		//continue loop if less than width
		
		add rax,rbx		//increase doffset
		jmp yloop		//next line 

lastline:
		xor rdx,rdx			//rdx =x=0
				
//first qword
		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2	//temp2=SC0[x]*2

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2  //temp1=SC1[x]+temp2
		paddw temp1,mm7		//temp1=SC1[x]+temp2+8  (+8 because off rounding)
		psrlw temp1,4		//temp1=(SC1[x]+temp2+8)/16
		movq [rcx+rdx-8],temp1		//dst=temp1
		
		add rdx,8			//increase counter

xloope:
		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2

		movq SC0(rdx),temp1	//SC0[x]=temp1

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2  //temp1=SC1[x]+temp2
		paddw temp1,mm7		//temp1=SC1[x]+temp2+8  (+8 because off rounding)
		psrlw temp1,4		//temp1=(SC1[x]+temp2+8)/16
		movq [rax+rdx-8],temp1		//dst=temp1
		
		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae xloope		//continue loop if less than width



		emms
		

	}
	#undef temp1
	#undef temp2
	#undef SR0
	#undef SR1

	#undef SC0

	#undef SC1
}





void CalcGaussianN5MMX(int width,int height,unsigned short* srcp,unsigned short* SC01,unsigned short* SC23)
{

const __int64 const8=0x0080008000800080;

const __int64 mask=  0x000000000000FFFF;

int endline=width*height;

	__asm{
		
		movq mm7,const8		
		movq mm6,mask
		

		mov ebx,width
		mov rcx,srcp
		mov rax,rcx
		add endline,eax
		xor rdx,rdx	//rdx=0
		mov rdi,SC01
		mov rsi,SC23
		//rcx=soffset
		//rdi SC0;
		//rsi SC1;
		//rdx=x
		//rax=doffset
		
		
		//mm0 temp1
		//mm1 temp2
		//mm2 SR0
		//mm3 SR1
#define temp1 mm0
#define temp2 mm1
#define SR0 mm2
#define SR1 mm3
#define SR2 mm4
#define SR3 mm5
#define pix1 mm7
#define pix2 mm6
#define CONST128 mm7
#define SC0(index) [rdi+(index)*2]
#define SC1(index) [rdi+(index)*2+8]
#define SC2(index) [rsi+(index)*2]
#define SC3(index) [rsi+(index)*2+8]
		//setup SR0, SR1,SR2 and SR3
		movq SR0,[rcx+rbx-1*8]	//Move last qword into SR0 SR0:ok = 1
		
		movq SR1,[rcx+rbx-2*8]	//Move the second last qword into SR1= 1 0
		movq SR2,SR1			//SR2 = 1 0
		paddw SR1,SR0			//add SR0 to SR1 SR1:OK		1+ 1 0=1 1
		
		movq SR3,[rcx+rbx-3*8]	//SR3= 1 0 0
		paddw SR2,SR3	//add third last qword to the second last SR2= 1 0 + 1 0 0= 1 1 0
		
		paddw SR3,[rcx+rbx-4*8] //add 4. last qword to SR3= 1 0 0+ 1 0 0 0= 1 1 0 0
		paddw SR3,SR2				//add SR2+SR3= 1 1 0+1 1 0 0 =  1 2 1 0 
		
		paddw SR2,SR1				//add SR1+SR2  SR2:ok= 1 1 + 1 1 0 = 1 2 1
		paddw SR3,SR2				//add SR2+SR3 SR3:ok= 1 2 1+ 1 2 1 0 = 1 3 3 1
		
		//setup SC0 and SC1
//first Qword
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq pix1,temp1		  //save first pixel
		psrlq temp1,16		//shift first pixel out(to calculate the last)
		
		movq mm6,[rcx+rbx-1*8]		//move SR0 in to extract last pixel
		psrlq mm6,48		
		psllq mm6,48		
		paddw temp1,mm6		//add last pixel		1

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2
		
		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9

		movq temp2,temp1
		paddw temp2,temp1		
		movq SC1(rdx),temp2		//SC1[x]=SC0+temp1	10+13

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,temp2	//temp1=SC1(temp2)
		paddw temp1,temp2	//temp1=SC1(temp2)+temp2
		movq SC2(rdx),temp1	//SC2=temp1				12+15

		movq temp2,temp1	//temp2=SC2(temp1)
		paddw temp2,temp1   //temp2=SC2(temp1)+temp1
		movq SC3(rdx),temp2 //SC3=temp2				14+17

		add rdx,8			//increase counter
		
//second qword
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq temp2,temp1     //save second pixel temporarely in temp2
		psrlq temp1,16		//shift first pixel out(to calculate the last)
		
		movq mm6,[rcx+rbx-1*8]		//move last qword in to extract last pixel
		psrlq mm6,48		
		psllq mm6,48		
		paddw temp1,mm6		//add last pixel		1
		movq pix2,temp2		//Move second pixel into mm6

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2
		
		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9

		movq temp2,temp1
		paddw temp2,temp1		
		movq SC1(rdx),temp2		//SC1[x]=SC0+temp1	10+13

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,temp2	//temp1=SC1(temp2)
		paddw temp1,temp2	//temp1=SC1(temp2)+temp2
		movq SC2(rdx),temp1		//SC2=temp1				12+15

		movq temp2,temp1	//temp2=SC2(temp1)
		paddw temp2,temp1   //temp2=SC2(temp1)+temp1
		movq SC3(rdx),temp2 //SC3=temp2				14+17


		pand pix1,mask		//select first pixel     1 0
		pand pix2,mask		//select second pixel     1
		psllq SR0,16		//align
		paddw SR0,mm6		//SR0:ok				1
		
		psllq SR1,16		//align
		paddw mm6,mm7			//mm6 = 1 1

		paddw SR1,mm6		//SR1:ok
		

		psllq SR2,16		//align
		paddw mm7,mm7			//mm7 = 1 1 0
		paddw mm6,mm7			//mm6 = 1 2 1
		paddw SR2,mm6		//SR2:ok

		psllq SR3,16		//align
		paddw mm7,mm7			//mm7 = 1 2 1 0
		paddw mm6,mm7			//mm6 = 1 3 3 1
		paddw SR3,mm6		//SR3:ok

		add rdx,8			//increase counter
initloop:
		movq temp1,[rcx+rdx]  //temp1=srcp[x]		1

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2
		
		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1			
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9

		movq temp2,temp1
		paddw temp2,temp1		
		movq SC1(rdx),temp2		//SC1[x]=SC0+temp1	10+13

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,temp2	//temp1=SC1(temp2)
		paddw temp1,temp2	//temp1=SC1(temp2)+temp2
		movq SC2(rdx),temp1		//SC2=temp1				12+15

		movq temp2,temp1	//temp2=SC2(temp1)
		paddw temp2,temp1   //temp2=SC2(temp1)+temp1
		movq SC3(rdx),temp2 //SC3=temp2				14+17
		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae initloop		//continue loop if less than width

		
		//rdx=index for SC0/SC1;

//********************************************************************************************************
//second line		

		add rcx,rbx			//soffset=next line
		
		
//setup SR0, SR1,SR2 and SR3
		movq SR0,[rcx+rbx-1*8]	//Move last qword into SR0 SR0:ok = 1
		
		movq SR1,[rcx+rbx-2*8]	//Move the second last qword into SR1= 1 0
		movq SR2,SR1			//SR2 = 1 0
		paddw SR1,SR0			//add SR0 to SR1 SR1:OK		1+ 1 0=1 1
		
		movq SR3,[rcx+rbx-3*8]	//SR3= 1 0 0
		paddw SR2,SR3	//add third last qword to the second last SR2= 1 0 + 1 0 0= 1 1 0
		
		paddw SR3,[rcx+rbx-4*8] //add 4. last qword to SR3= 1 0 0+ 1 0 0 0= 1 1 0 0
		paddw SR3,SR2				//add SR2+SR3= 1 1 0+1 1 0 0 =  1 2 1 0 
		
		paddw SR2,SR1				//add SR1+SR2  SR2:ok= 1 1 + 1 1 0 = 1 2 1
		paddw SR3,SR2				//add SR2+SR3 SR3:ok= 1 2 1+ 1 2 1 0 = 1 3 3 1


		xor rdx,rdx			//rdx =x=0
		
//first Qword	
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq pix1,temp1		  //save first pixel
		psrlq temp1,16		//shift first pixel out(to calculate the last)
		
		movq mm6,[rcx+rbx-1*8]		//move SR0 in to extract last pixel
		psrlq mm6,48		
		psllq mm6,48		
		paddw temp1,mm6		//add last pixel		1

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2
		
		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9


		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]	10

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		/*movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rcx+rdx-8],temp1		//dst=temp1*/
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17
		

		add rdx,8			//increase counter

//second qword
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq temp2,temp1     //save second pixel temporarely in temp2
		psrlq temp1,16		//shift first pixel out(to calculate the last)
		
		movq mm6,[rcx+rbx-1*8]		//move last qword in to extract last pixel
		psrlq mm6,48		
		psllq mm6,48		
		paddw temp1,mm6		//add last pixel		1
		movq pix2,temp2		//Move second pixel into mm6

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2

		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9


		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]	10

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		/*movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rcx+rdx-8],temp1		//dst=temp1*/
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17

		pand pix1,mask		//select first pixel     1 0
		pand pix2,mask		//select second pixel     1
		psllq SR0,16		//align
		paddw SR0,mm6		//SR0:ok				1

		psllq SR1,16		//align
		paddw mm6,mm7			//mm6 = 1 1

		paddw SR1,mm6		//SR1:ok
		

		psllq SR2,16		//align
		paddw mm7,mm7			//mm7 = 1 1 0
		paddw mm6,mm7			//mm6 = 1 2 1
		paddw SR2,mm6		//SR2:ok

		psllq SR3,16		//align
		paddw mm7,mm7			//mm7 = 1 2 1 0
		paddw mm6,mm7			//mm6 = 1 3 3 1
		paddw SR3,mm6		//SR2:ok

		add rdx,8			//increase counter




xloop1:
		movq temp1,[rcx+rdx]  //temp1=srcp[x]		1

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2

		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9


		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]	10

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		/*movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rcx+rdx-8],temp1		//dst=temp1*/
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17

		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae xloop1		//continue loop if less than width
		
		

//************************************************************************
//main loop





yloop:
		
		add rcx,rbx			//soffset=next line
		cmp ecx,endline
		jae lastline
//setup SR0, SR1,SR2 and SR3
		movq SR0,[rcx+rbx-1*8]	//Move last qword into SR0 SR0:ok = 1
		
		movq SR1,[rcx+rbx-2*8]	//Move the second last qword into SR1= 1 0
		movq SR2,SR1			//SR2 = 1 0
		paddw SR1,SR0			//add SR0 to SR1 SR1:OK		1+ 1 0=1 1
		
		movq SR3,[rcx+rbx-3*8]	//SR3= 1 0 0
		paddw SR2,SR3	//add third last qword to the second last SR2= 1 0 + 1 0 0= 1 1 0
		
		paddw SR3,[rcx+rbx-4*8] //add 4. last qword to SR3= 1 0 0+ 1 0 0 0= 1 1 0 0
		paddw SR3,SR2				//add SR2+SR3= 1 1 0+1 1 0 0 =  1 2 1 0 
		
		paddw SR2,SR1				//add SR1+SR2  SR2:ok= 1 1 + 1 1 0 = 1 2 1
		paddw SR3,SR2				//add SR2+SR3 SR3:ok= 1 2 1+ 1 2 1 0 = 1 3 3 1


		xor rdx,rdx			//rdx =x=0
		
//first Qword	
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq pix1,temp1		  //save first pixel
		psrlq temp1,16		//shift first pixel out(to calculate the last)
		
		movq mm6,[rcx+rbx-1*8]		//move SR0 in to extract last pixel
		psrlq mm6,48		
		psllq mm6,48		
		paddw temp1,mm6		//add last pixel		1

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2
		
		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9


		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]	10

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rbx-16],temp1		//dst=temp1
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17
		

		add rdx,8			//increase counter

//second qword
		movq temp1,[rcx+rdx]  //temp1=srcp[x]
		movq temp2,temp1     //save second pixel temporarely in temp2
		psrlq temp1,16		//shift first pixel out(to calculate the last)
		
		movq mm6,[rcx+rbx-1*8]		//move last qword in to extract last pixel
		psrlq mm6,48		
		psllq mm6,48		
		paddw temp1,mm6		//add last pixel		1
		movq pix2,temp2		//Move second pixel into mm6

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2

		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9


		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]	10

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rbx-8],temp1		//dst=temp1
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17

		pand pix1,mask		//select first pixel     1 0
		pand pix2,mask		//select second pixel     1
		psllq SR0,16		//align
		paddw SR0,mm6		//SR0:ok				1

		psllq SR1,16		//align
		paddw mm6,mm7			//mm6 = 1 1

		paddw SR1,mm6		//SR1:ok
		

		psllq SR2,16		//align
		paddw mm7,mm7			//mm7 = 1 1 0
		paddw mm6,mm7			//mm6 = 1 2 1
		paddw SR2,mm6		//SR2:ok

		psllq SR3,16		//align
		paddw mm7,mm7			//mm7 = 1 2 1 0
		paddw mm6,mm7			//mm6 = 1 3 3 1
		paddw SR3,mm6		//SR2:ok

		add rdx,8			//increase counter



		movq CONST128,const8
xloop:
		movq temp1,[rcx+rdx]  //temp1=srcp[x]		1

		movq temp2,SR0		//temp2 SR0
		paddw temp2,temp1		//temp2=temp1+SR0	2

		movq SR0,temp1		//SR0=temp1				3

		movq temp1,SR1		//temp1=SR1
		paddw temp1,temp2	//temp1=SR1+temp2		4

		movq SR1,temp2		//SR1=temp2				5

		movq temp2,SR2		//temp2=SR2
		paddw temp2,temp1   //temp2=SR2+temp1		6	

		movq SR2,temp1		//SR2=temp1				7

		movq temp1,SR3		//temp1=SR3
		paddw temp1,temp2   //temp1=SR3+temp2		8

		movq SR3,temp2		//SR3=temp2				9


		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp1	//temp2=temp1+SCO[x]	10

		movq SC0(rdx),temp1	//SC0[x]=temp1			11

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,CONST128		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC3[x]+temp2+8)/256
		movq [rax+rdx-16],temp1		//dst=temp1
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17

		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae xloop		//continue loop if less than width
		
		add rax,rbx		//increase doffset
		jmp yloop		//next line 

//*************************************************************************************************************
//two last lines

lastline:
		xor rdx,rdx			//rdx =x=0
				
//first qword
		
		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2	//temp2=2*SCO[x]		10

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rbx-16],temp1		//dst=temp1
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17
		

		add rdx,8			//increase counter

//second qword

		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2	//temp2=temp1+SCO[x]	10

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rbx-8],temp1		//dst=temp1
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17

		add rdx,8			//increase counter
xloope:
		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2	//temp2=temp1+SCO[x]	10

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq SC1(rdx),temp2 //SC1(x)=temp2			13

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq SC2(rdx),temp1	//SC2[x]=temp1			15

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rdx-16],temp1		//dst=temp1
		
		movq SC3(rdx),temp2		//SC3[x]=temp2		17

		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae xloope		//continue loop if less than width
		
		add rax,rbx		//increase doffset

//Last Line
		

		xor rdx,rdx			//rdx =x=0
				
//first qword
		
		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2	//temp2=temp1+SCO[x]	10

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rbx-16],temp1		//dst=temp1
		
		add rdx,8			//increase counter

//second qword

		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2	//temp2=temp1+SCO[x]	10

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14

		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rbx-8],temp1		//dst=temp1
		
		add rdx,8			//increase counter
xloopend:
		movq temp2,SC0(rdx)	//temp2=SC0[x]
		paddw temp2,temp2	//temp2=temp1+SCO[x]	10

		movq temp1,SC1(rdx)	//temp1=SC1[x]
		paddw temp1,temp2	//temp1=SC1[x]+temp2	12

		
		movq temp2,SC2(rdx)	//temp2=SC2[x]
		paddw temp2,temp1	//temp2=temp1+SC2[x]	14
	
		movq temp1,SC3(rdx) //temp1=SC3[x]
		paddw temp1,temp2  //temp1=SC3[x]+temp2
		paddw temp1,const8		//temp1=SC3[x]+temp2+128  (+128 because off rounding)
		psrlw temp1,8		//temp1=(SC1[x]+temp2+8)/256
		movq [rax+rdx-16],temp1		//dst=temp1

		add rdx,8			//increase counter
		cmp rdx,rbx			//compare with width
		jnae xloopend		//continue loop if less than width




		emms
		

	}
	#undef temp1
	#undef temp2
	#undef SR0
	#undef SR1
	#undef SR2
	#undef SR3

	#undef SC0
	#undef SC1
	#undef SC2
	#undef SC3

	#undef pix1
	#undef pix2
}