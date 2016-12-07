# Lin-Analyzer

A high-level performance analysis tool for FPGA-based accelerator.

## Installation
### Linux Installation

Assume you have the following directory structure (windows setup is similar):
```
$HOME
   ~/llvm
   ~/llvm/tools/clang
   ~/llvm/tools/lin-analyzer
   ~/build
   ~/boost_1_57_0
```

* LLVM and clang 3.5: 

```
wget http://llvm.org/releases/3.5.0/llvm-3.5.0.src.tar.xz
wget http://llvm.org/releases/3.5.0/cfe-3.5.0.src.tar.xz
tar -xvf llvm-3.5.0.src.tar.xz
mv llvm-3.5.0.src llvm
tar -xvf cfe-3.5.0.src.tar.xz
mv cfe-3.5.0.src clang
cp -r clang llvm/tools
```

* Lin-Analyzer:

```
cd ~/llvm/tools
git clone git@github.com:zhguanw/lin-analyzer.git
open ~/llvm/tools/CMakeLists.txt and add "add_llvm_tool_subdirectory(lin-analyzer)" to it
open ~/llvm/CMakeLists.txt and add "set(LLVM_REQUIRES_RTTI 1)" to enable RTTI feature
```

* Boost graph library: 

```
cd ~
wget http://softlayer-sng.dl.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.gz
tar -xvf boost_1_57_0.tar.gz
```

* zlib library, version 1.2.8

```
sudo apt-get install zlib1g-dev
```

* cmake

```
sudo apt-get install cmake
```
(I use version 2.8.12.2, higher version should be fine)

* Install Lin-Analyzer
```
cd ~
mkdir build
cd ~/build
cmake ~/llvm -DBOOST_INCLUDE_DIR=/your-path-to/boost_1_57_0 -DZLIB_INCLUDE_DIRS=/usr/include -DZLIB_LIBRARY=/usr/lib/x86_64-linux-gnu/libz.so
```

### Windows Installation
I use cygwin, so the first three steps above are the same.

* Install Visual Studio 2013 (32-bit or 64-bit)
* zlib 1.2.8:
```
(a). Download zlib from http://zlib.net/zlib-1.2.8.tar.gz and uncompress it
(b). vim your_path_to\zlib-1.2.8\contrib\masmx86\bld_ml32.bat
     ml.exe /coff /Zi /c /Flmatch686.lst match686.asm
     ml.exe /coff /Zi /c /Flinffas32.lst inffas32.asm
(c). Set Visual Studio to the system environment:
     PATH = %PATH:your_path_to\Microsoft Visual Studio 13.0\VC\bin
(d). Open zlibvc.sln using vs2013 located at your_path_to\zlib-1.2.8\contrib\vstudio\vc11
(e). Open properties of zlibstat project and remove "ZLIB_WINAPI;" from  "Configuration Properties -> C/C++ -> Preprocessor -> Preprocessor Definitions"
(f). If you generate 64-bit lin-analyzer, then you also need to generate 64-bit zlib library. To do this, you need to Open "Configuration Manager" and change "Win32" to "x64" under "Platform"
(g). Build the static library. 
(h). After successfully building the zlib static library, you need to specify its absolute path for ZLIB_INCLUDE_DIRS and ZLIB_LIBRARY variables in cmake.
```

* Install cmake 2.8.12.2
```
Use cmake to compile llvm and setup Visual Studio project files
Will update detailed instructions soon...
```

## Getting started
* Lin-Analyzer options
```
Usage:	lin-analyzer [file.bc] -Ipath=[path] -config=[filename] [kernel-name] -Opath=[path] -TargetLoops=[index] [options]
Options:
	-h, --help               Help information.
	-profiling-time          Profile kernels without FPGA estimation.
	-mem-trace               Obtain memory trace for access pattern analysis without FPGA estimation. This
	                         option should be used with -profiling-time together.
	-no-trace                Disable dynamic trace generation.
	-TargetLoops             Specify target loops focused. Eg., -TargetLoops=2,3: only analyse loop 2 and 3.
	-cfg-detailed            Show CFG with detailed instructions.
	-cfg-only                Show CFG only with basic blocks.
	-dddg-bf-opt             Show DDDG before optimization. May slow down program if input size is large
	-dddg-af-opt             Show DDDG after optimization. May slow down program if input size is large
	-verbose                 Verbose mode, print more information.
	-dis-store-buffer        Disable store-buffer optimization.
	-shared-load-removal     Enable shared-load-removal optimization.
	-dis-shared-load-removal Disable shared-load-removal opt., even for completely unrolling config.
	-dis-rp-store-removal    Disable repeated-store-removal optimization.
	-THR-float               Enable tree height reduction optimization for floating point operations.
	-THR-integer             Enable tree height reduction optimization for integer operations.
	-memory-disambig         Enable memory disambiguation optimization.
	-dis-fp-threshold        Disable floating point unit threshold and area budget is unlimited.
	-en-extra-scalar         Sometimes, result might be shifted, this option is used to improve prediction.
	-en-rw-rw-memory         Use memory with two ports and each supports read and write operation. Default is
	                         read-only and read-write.
	-vc707                   Target for Xilinx Virtex7 VC707. Default is Xilinx Zedboard or ZC702
```

2. Design Space Exploration
```
export BOOST_ROOT=~/boost_1_57_0
export LD_LIBRARY_PATH=~/build/bin/lib
export LLVM_BIN_HOME=~/build/bin
export LLVM_SRC=~/llvm
export TESTBENCH_HOME=~/llvm/tools/lin-analyzer/testsuite/Ecobench
cd ~/llvm/tools/lin-analyzer/testsuite/Ecobench/scripts
python run_dse.py
```

3. Windows
```
Similar to Linux command, but use Visual Studio;
Will update it soon...
```

## License

Lin-Analyzer is licensed under the GPL-3.0 license.

Dynamic data graph generation (DDDG) and few of optimization functions in Lin-Analyzer are modified from Aladdin project and Lin-Analyzer is also followed Aladdin's license. For more information, please visit Aladdin's website: https://github.com/ysshao/ALADDIN

## Citation

If you use Lin-Analyzer in your research, please cite

Lin-Analyzer: A High-level Performance Analysis Tool for FPGA-based Accelerators,
Guanwen Zhong, Alok Prakash,Yun Liang, Tulika Mitra, Smail Niar,
53rd ACM/IEEE Design Automation Conference, June 2016

==========================
Guanwen (Henry) Zhong,

guanwen@comp.nus.edu.sg

National University of Singapore, 2016
