LLVM version: 07f1b21cb5c9ccda085f810c2c9d0d154e29f8c9
CMake version: 2.8.12.2
visual studio: 2013

Input File:
```
#!c++
clang -emit-llvm -c test.c
opt -mem2reg -instnamer -lcssa -indvars test.bc -o test1.bc
```
Check llvm IR:
```
#!c++
llvm-dis test1.bc
```
The method to compile multiple source codes:

Assume we have 2 source files, main.c and sum.c.
-----------
main.c:
```
#!c++
#include <stdio.h>
int sum(int x, int y);
int main() {
   int r = sum(3, 4);
   printf("r = %d\n", r);
   return 0;
}
```
sum.c:
```
#!c++
int sum(int x, int y) {
   return x+y;
}
```
Solution:

```
#!c++

clang -g -O3 -emit-llvm -c main.c -o main.bc
clang -g -O1 -emit-llvm -c sum.c -o sum.bc
llvm-link main.bc sum.bc -o sum.linked.bc
opt -mem2reg -instnamer -lcssa -indvars sum.linked.bc -o sum.final.bc
```

How to run lin-profile:

```
#!c++

When we use visual studio 2013 to run lin-profile, we need to specify the input arguments as:

input_bitcode.bc path_of_input_bitcode kernel_function_name kernel_subfunction_name

eg, for stencil3d benchmark:

E:\llvm\home\henry\llvm_3.5\input\stencil3d\stencil3d_app.bc E:\llvm\home\henry\llvm_3.5\input\stencil3d\ stencil3d
```
---------------------------------------------
Download boost graph library:

```
#!c++

wget http://softlayer-sng.dl.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.gz
```


---------------------------------------------
Install zlib library in linux:

```
#!c++

sudo apt-get install zlib1g-dev
```


Build zlib library in visual studio 2013:
1. Download zlib 1.2.8:  

```
#!c++

wget http://zlib.net/zlib-1.2.8.tar.gz
```
2. change your_path_to\zlib-1.2.8\contrib\masmx86\bld_ml32.bat:

```
#!c++

ml.exe /coff /Zi /c /Flmatch686.lst match686.asm
ml.exe /coff /Zi /c /Flinffas32.lst inffas32.asm
```
3. Make sure you set the path of system environment with: (this is used for finding ml.exe)

```
#!c++

PATH = %PATH:your_path_to\Microsoft Visual Studio 13.0\VC\bin
```
4. Open zlibvc.sln using vs2013 located at your_path_to\zlib-1.2.8\contrib\vstudio\vc11
5. Open properties of zlibstat project and remove "ZLIB_WINAPI;" from  "Configuration Properties -> C/C++ -> Preprocessor -> Preprocessor Definitions"
6. If we generate 64-bit lin-profiler, then we also need to generate 64-bit zlib library. To do this, we need to Open "Configuration Manager" and change "Win32" to "x64" under "Platform"
7. Build the static library. 
8. After successfully building the zlib static library, we need to specify its absolute path for ZLIB_INCLUDE_DIRS and ZLIB_LIBRARY variables.

Done!

This project uses Cmake to link library automatically, if users want to link libraries manually, you can check this link: [http://www.codeproject.com/Articles/85391/Microsoft-Visual-C-Static-and-Dynamic-Libraries](http://www.codeproject.com/Articles/85391/Microsoft-Visual-C-Static-and-Dynamic-Libraries)
---------------------------------------------

Before we use cmake to configure the environment, we need to edit the CMakeLists.txt located at $LLVM_SOURCE_DIR. Just add "set(LLVM_REQUIRES_RTTI 1)" to enable RTTI feature, because of boost graph library we used.

To build lin-profiler in Linux or Unix
Follow the steps in setup.h


```
#!c++
In build folder, using the following commands for building llvm with lin-profiler
cmake ../llvm-3.5.0.src/ -DENABLE_TESTSUITE=ON -DBOOST_INCLUDE_DIR=/home/guanwen1/boost_1_57_0 -DZLIB_INCLUDE_DIRS=/usr/include -DZLIB_LIBRARY=/usr/lib/x86_64-linux-gnu/libz.so
```

loop pipelining configuration format:
pipeline,function name,loop number,loop level

The architecture of Lin-analyzer:
![Lin-Analyzer_Architecture.jpg](https://bitbucket.org/repo/qokpMA/images/1538122950-Lin-Analyzer_Architecture.jpg)

Add -profiling-time support:
When specifying this option, Lin-anlyzer will only run profiling without estimating FPGA execution cycles. This is used to figure out how much time spend on profiling step.