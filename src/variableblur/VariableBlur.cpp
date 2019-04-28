// VariableBlur.cpp : Defines the entry point for the DLL application.
//

#include "GaussianBlur.h"

AVSValue __cdecl Create_binomialblur(AVSValue args, void* user_data, IScriptEnvironment* env){
	unsigned int radY=args[1].AsFloat(1.5)*2+0.1;
	unsigned int radC=args[2].AsFloat(1.5)*2+0.1;
	return new BinomialBlur(args[0].AsClip(),
		radY,
		radC,
		args[3].AsInt(3),
		args[4].AsInt(2),
		args[5].AsInt(2),
		true,
		args[6].AsBool(true),
		args[6].Defined(),
		true,
		env);
}

AVSValue __cdecl Create_averageblur(AVSValue args, void* user_data, IScriptEnvironment* env){
	return new BinomialBlur(args[0].AsClip(),
		args[1].AsInt(3),
		args[2].AsInt(3),
		args[3].AsInt(3),
		args[4].AsInt(2),
		args[5].AsInt(2),
		false,
		true,
		true,
		true,
		env);
}


AVSValue __cdecl Create_gaussianblur(AVSValue args, void* user_data, IScriptEnvironment* env){
	double sdY = sqrt(args[1].AsFloat(1));
	double sdC = sqrt(args[2].AsFloat(1));
	return new GaussianBlur(args[0].AsClip(),
		sdY,
		sdC,
        args[3].AsInt(4),
		args[4].AsBool(false),
		0,
		args[5].AsInt(3),
		args[6].AsInt(2),
		args[7].AsInt(2),
		args[8].AsInt(-1),
		args[9].AsInt(-1),
		args[10].AsInt(0),
		args[11].AsInt(0),
		args[12].AsInt(1),
		env);
}

AVSValue __cdecl Create_unsharp(AVSValue args, void* user_data, IScriptEnvironment* env){
	double sdY = sqrt(args[1].AsFloat(1));
	double sdC = sqrt(args[2].AsFloat(1));
	return new GaussianBlur(args[0].AsClip(),
		sdY,
		sdC,
		args[4].AsInt(4),
		args[5].AsBool(false),
		args[3].AsFloat(0.7),
		args[6].AsInt(3),
		args[7].AsInt(2),
		args[8].AsInt(2),
		args[9].AsInt(-1),
		args[10].AsInt(-1),
		args[11].AsInt(0),
		args[12].AsInt(0),
		args[13].AsInt(1),
		env);
}

const AVS_Linkage *AVS_linkage;

extern "C" __declspec(dllexport)
const char * __stdcall AvisynthPluginInit3(IScriptEnvironment *env, const AVS_Linkage *const vectors)
{
	AVS_linkage = vectors;
	
env->AddFunction("BinomialBlur","c[varY]f[varC]f[Y]i[U]i[V]i[usemmx]b", Create_binomialblur, 0); 
env->AddFunction("AverageBlur","c[radY]i[radC]i[Y]i[U]i[V]i", Create_averageblur, 0); 
env->AddFunction("GaussianBlur","c[varY]f[varC]f[border]i[integrate]b[Y]i[U]i[V]i[gfunc]i[gfuncc]i[pcr]i[pcrc]i[nthreads]i", Create_gaussianblur, 0); 
env->AddFunction("Unsharp","c[varY]f[varC]f[strength]f[border]i[integrate]b[Y]i[U]i[V]i[gfunc]i[gfuncc]i[pcr]i[pcrc]i[nthreads]i", Create_unsharp, 0); 
return "VariableBlur version 0.7"; 
} 




