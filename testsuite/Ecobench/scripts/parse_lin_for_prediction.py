#!/usr/bin/env python

import numpy as np
import sys
import os
import os.path
import math
import matplotlib
import run_analyzer
import time

matplotlib.use('Agg')
import matplotlib.pyplot as plt


#unroll = [1, 2, 4, 8, 16, 32, 64]
#pipe = [0,1]

'''
input_sizes = {
'atax' : 128,
'bicg' : 256,
'convolution2d' : 128,
'convolution3d' : 32,
'gemm' : 128,
'gesummv' : 128,
'mm' : 128,
'mvt' : 256,
'syr2k' : 128,
'syrk' : 128,
}
'''

numLplevels = {
'atax' : 2,
'bicg' : 2,
'convolution2d' : 2,
'convolution3d' : 3,
'gemm' : 3,
'gesummv' : 2,
'mm' : 3,
'mvt' : 2,
'syr2k' : 3,
'syrk' : 3,
}

#'convolution2d' : [1, 2, 7, 14, 49, 98],
#'convolution3d' : [1, 2, 3, 10, 15, 30],

unrolls = {
'atax' : [1, 2, 4, 8, 16, 32, 64, 128],
'bicg' : [1, 2, 4, 8, 16, 32, 64, 128, 256],
'convolution2d' : [1, 2, 3, 6, 7, 9, 14, 18, 21, 42, 63, 126],
'convolution3d' : [1, 2, 3, 5, 6, 10, 15, 30],
'gemm' : [1, 2, 4, 8, 16, 32, 64, 128],
'gesummv' : [1, 2, 4, 8, 16, 32, 64, 128],
'mm' : [1, 2, 4, 8, 16, 32, 64, 128],
'mvt' : [1, 2, 4, 8, 16, 32, 64, 128, 256],
'syr2k' : [1, 2, 4, 8, 16, 32, 64, 128],
'syrk' : [1, 2, 4, 8, 16, 32, 64, 128],
}

# 0 means no loop pipelining applied; 1 means loop pipelining is applied to loop level 1
# 2 means loop pipelining is applied to loop level 2, ...
pipelines = {
'atax' : ['0', '1', '2'],
'bicg' : ['0', '1', '2'],
'convolution2d' : ['0', '1', '2'],
'convolution3d' : ['0', '2', '3'],
'gemm' : ['0', '2', '3'],
'gesummv' : ['0', '1', '2'],
'mm' : ['0', '2', '3'],
'mvt' : ['0', '1', '2'],
'syr2k' : ['0', '2', '3'],
'syrk' : ['0', '2', '3'],
}

partitions = [1, 2, 4, 8, 16]

#You can parallel this by running multiply jobs together
#In this case, llvm_coompile.py inside run_aladdin.py only need to run once
#to generate the dynamic instruction trace

#for f_unroll in unroll:
  #for f_part in part:
    #for f_pipe in pipe:
      #os.system('python run_aladdin.py %s %i %i %i' % (bench, f_part, f_unroll, f_pipe))

benchmarks = ['atax', 'bicg', 'convolution2d', 'convolution3d', 'gemm', 'gesummv', 'mm', 'mvt', 'syr2k', 'syrk']
#benchmarks = ['mm']
upperlevel = {
'atax' : '1',
'bicg' : '1',
'convolution2d' : '1',
'convolution3d' : '2',
'gemm' : '2',
'gesummv' : '1',
'mm' : '2',
'mvt' : '1',
'syr2k' : '2',
'syrk' : '2',
}

# profiling time unit: second
profiling_times = {
'atax' : 0,
'bicg' : 0,
'convolution2d' : 0,
'convolution3d' : 0,
'gemm' : 0,
'gesummv' : 0,
'mm' : 0,
'mvt' : 0,
'syr2k' : 0,
'syrk' : 0,
}

# exploration time unit: second
exploration_times = {
'atax' : 0,
'bicg' : 0,
'convolution2d' : 0,
'convolution3d' : 0,
'gemm' : 0,
'gesummv' : 0,
'mm' : 0,
'mvt' : 0,
'syr2k' : 0,
'syrk' : 0,
}

if not 'TESTBENCH_HOME' in os.environ:
  raise Exception('Set TESTBENCH_HOME directory as an environment variable')
	
TESTBENCH_HOME = os.getenv('TESTBENCH_HOME')
RESULT_DIRECTORY = TESTBENCH_HOME+'/scripts/results'

# Write profiling and exploration time into csv file
os.chdir(TESTBENCH_HOME)

ML_input_DIR = TESTBENCH_HOME + '/scripts/ML_input'
if not os.path.isdir(ML_input_DIR):
  os.mkdir(ML_input_DIR)

# XILINX ZC702 Resources
max_dsp = 220
max_bram18k = 280
max_ff = 106400
max_lut = 53200

for kernel in benchmarks:
  print 'Kernel name: ' + kernel
  
  fpga_cycles = []
  area = []
  dsp_used = []
  bram18k_used = []
  ff_used = []
  lut_used = []
  fadd_used = []
  fsub_used = []
  fmul_used = []
  fdiv_used = []
	
  asap_il = []
  ave_par = []
  maxAccPerBank = []
  maxAveLdPerBank = []
  maxAveStPerBank = []
	
  ld_inst = []
  st_inst = []
  fadd_inst = []
  fsub_inst = []
  fmul_inst = []
  fdiv_inst = []
  fcmp_inst = []
  int_inst = []
  bit_inst = []
  cont_inst = []
  br_inst = []
  
  KERNEL_HOME = TESTBENCH_HOME + '/' + kernel
  os.chdir(KERNEL_HOME+'/sim')
  dirs = os.listdir(KERNEL_HOME + '/sim/')
 
  predict_ff_lut_filename = ML_input_DIR + '/' + kernel + '_' + str(unrolls[kernel][-1]) + '_lin_predict_ff_lut.csv'
  predict_ff_lut = open(predict_ff_lut_filename, 'w')
#  predict_ff_lut.write('InputSize,NumberLpL,Partition,Unroll,PipeUpperLpL,PipeInnermostLpL,dsp,bram18k,usedFadd,usedFsub,usedFmul,usedFdiv\n')
  predict_ff_lut.write('InputSize,NumberLpL,Partition,Unroll,PipeUpperLpL,PipeInnermostLpL,\
dsp,bram18k,asapIL,avePar,maxAccPerBank,maxAveLdPerBank,maxAveStPerBank,ldInst,stInst,\
faddInst,fsubInst,fmulInst,fdivInst,fcmpInst,intInst,bitInst,contInst,brInst,LinCycles\n')
  for dir in dirs:
    config_name = dir.split('/')[-1]
    #print config_name
    sum_file_name = dir+'/'+kernel+'_summary.log'
    if os.path.isfile(sum_file_name):
      sum_file = open(sum_file_name, 'r')
      for line in sum_file:
        if 'cycles' in line:
          fpga_cycles.append(int(line.split(' ')[3]))
        elif 'DSP used' in line:
          dsp_used.append(int(line.split(' ')[2]))
        elif 'BRAM18K used' in line:
          bram18k_used.append(int(line.split(' ')[2]))
        elif 'FF used' in line:
          ff_used.append(int(line.split(' ')[2]))
        elif 'LUT used' in line:
          lut_used.append(int(line.split(' ')[2]))
        elif 'Fadd used' in line:
          fadd_used.append(int(line.split(' ')[2]))
        elif 'Fsub used' in line:
          fsub_used.append(int(line.split(' ')[2]))
        elif 'Fmul used' in line:
          fmul_used.append(int(line.split(' ')[2]))
        elif 'Fdiv used' in line:
          fdiv_used.append(int(line.split(' ')[2]))
        elif 'ideal iteration latency' in line:
          asap_il.append(int(line.split(':')[1]))
        elif 'average parallelism' in line:
          ave_par.append(float(line.split(':')[1]))
        elif 'maximum number of memory access to a memory bank' in line:
          maxAccPerBank.append(float(line.split(':')[1]))
        elif 'maximum average number of memory bank load access' in line:
          maxAveLdPerBank.append(float(line.split(':')[1]))
        elif 'maximum average number of memory bank store access' in line:
          maxAveStPerBank.append(float(line.split(':')[1]))
        elif 'load instruction num' in line:
          ld_inst.append(int(line.split(':')[1]))
        elif 'store instruction num' in line:
          st_inst.append(int(line.split(':')[1]))
        elif 'fadd instruction num' in line:
          fadd_inst.append(int(line.split(':')[1]))
        elif 'fsub instruction num' in line:
          fsub_inst.append(int(line.split(':')[1]))
        elif 'fmul instruction num' in line:
          fmul_inst.append(int(line.split(':')[1]))
        elif 'fdiv instruction num' in line:
          fdiv_inst.append(int(line.split(':')[1]))
        elif 'fcmp instruction num' in line:
          fcmp_inst.append(int(line.split(':')[1]))
        elif 'integer instruction num' in line:
          int_inst.append(int(line.split(':')[1]))
        elif 'bitwise instruction num' in line:
          bit_inst.append(int(line.split(':')[1]))
        elif 'control instruction num' in line:
          cont_inst.append(int(line.split(':')[1]))
        elif 'branch instruction num' in line:
          br_inst.append(int(line.split(':')[1]))
          #break
        else:
          continue
      sum_file.close() 
  	
  os.chdir(RESULT_DIRECTORY)
  
  area = [float(dsp)/max_dsp + float(bram18k)/max_bram18k for dsp, bram18k in zip(dsp_used, bram18k_used)]
	
  InputSize = unrolls[kernel][-1] #input_sizes[kernel]
  NumberLpL = numLplevels[kernel] 
  zipped = zip(dirs,dsp_used,bram18k_used,asap_il,ave_par,maxAccPerBank,maxAveLdPerBank,maxAveStPerBank,ld_inst,st_inst,fadd_inst,fsub_inst,fmul_inst,fdiv_inst,fcmp_inst,int_inst,bit_inst,cont_inst,br_inst,fpga_cycles)

  for dir_t, dsp_u, bram18k_u, asap_il_t, ave_par_t, maxAccPB_t, maxAveLdPB_t, maxAveStPB_t, ld, st, fadd, fsub, fmul, fdiv, fcmp, integer, bitwise, control, branch, fpga_cycle in zipped:
    config_t = dir_t.split('/')[-1]
    pf_str = config_t.split('_')[0]
    uf_str = config_t.split('_')[1]
    pipe_str = config_t.split('_')[2]
    pf = pf_str.split('p')[-1]
    uf = uf_str.split('u')[-1]
    pipe = pipe_str.split('P')[-1]
    if pipe == '0':
      PipeUpperLpL = 0
      PipeInnermostLpL = 0
    elif pipe == str(NumberLpL):
      PipeUpperLpL = 0
      PipeInnermostLpL = 1
    elif pipe == str(NumberLpL -1):
      PipeUpperLpL = 1
      PipeInnermostLpL = 0
    else:
      PipeUpperLpL = 0
      PipeInnermostLpL = 0
      print 'Current design space only consider loop pipelining applied at the innermost two loop levels\n'			
    #predict_ff_lut.write(config_t + ',' + str(area_t) + '\n')
    #InputSize,NumberLpL,Partition,Unroll,PipeUpperLpL,PipeInnermostLpL,dsp,bram18k,usedFadd,usedFsub,usedFmul,usedFdiv
    predict_ff_lut.write(str(InputSize)+','+str(NumberLpL)+','+pf+','+uf+','+\
		str(PipeUpperLpL)+','+str(PipeInnermostLpL)+','+str(dsp_u)+','+str(bram18k_u)\
		+','+str(asap_il_t)+','+str(ave_par_t)+','+str(maxAccPB_t)+','+str(maxAveLdPB_t)\
		+','+str(maxAveStPB_t)+','+str(ld)+','+str(st)+','+str(fadd)+','+str(fsub)\
		+','+str(fmul)+','+str(fdiv)+','+str(fcmp)+','+str(integer)+','+str(bitwise)\
		+','+str(control)+','+str(branch)+','+str(fpga_cycle)+'\n')
  predict_ff_lut.close()

print 'Done!'

'''	
  for dir, dsp_u, bram18k_u, ff_u, lut_u, area_t, fadd_u, fsub_u, fmul_u, fdiv_u in zip(dirs,dsp_used,bram18k_used,ff_used,lut_used,area,fadd_used,fsub_used,fmul_used,fdiv_used):
    config_t = dir.split('/')[-1]
    pf_str = config_t.split('_')[0]
    uf_str = config_t.split('_')[1]
    pipe_str = config_t.split('_')[2]
    pf = pf_str.split('p')[-1]
    uf = uf_str.split('u')[-1]
    pipe = pipe_str.split('P')[-1]
    if pipe == '0':
      PipeUpperLpL = 0
      PipeInnermostLpL = 0
    elif pipe == str(NumberLpL):
      PipeUpperLpL = 0
      PipeInnermostLpL = 1
    elif pipe == str(NumberLpL -1):
      PipeUpperLpL = 1
      PipeInnermostLpL = 0
    else:
      PipeUpperLpL = 0
      PipeInnermostLpL = 0
      print 'Current design space only consider loop pipelining applied at the innermost two loop levels\n'			
    #predict_ff_lut.write(config_t + ',' + str(area_t) + '\n')
    #InputSize,NumberLpL,Partition,Unroll,PipeUpperLpL,PipeInnermostLpL,dsp,bram18k,usedFadd,usedFsub,usedFmul,usedFdiv
    predict_ff_lut.write(str(InputSize)+','+str(NumberLpL)+','+pf+','+uf+','+\
		str(PipeUpperLpL)+','+str(PipeInnermostLpL)+','+str(dsp_u)+','+str(bram18k_u)\
		+','+str(fadd_u)+','+str(fsub_u)+','+str(fmul_u)+','+str(fdiv_u)+'\n')
'''

