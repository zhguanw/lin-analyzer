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
#benchmarks = ['convolution2d']
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

designPointNum = {
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
for kernel in benchmarks:
  run_analyzer.gen_bitcode(TESTBENCH_HOME, kernel)


## Copy input bitcode files to Ecobench folder and generate dynamic trace
#BITCODE_HOME = os.getenv('ECOBENCH_BC_PATH')
for kernel in benchmarks:
  os.chdir(TESTBENCH_HOME+'/'+kernel)
  if os.path.isdir(TESTBENCH_HOME+'/'+kernel+'/sim'):
    os.system('rm -rf sim')
  
  #bitcode_file = BITCODE_HOME+'/'+kernel+'/'+kernel+'_linked_opt.bc'
  #os.system('cp %s %s' % (bitcode_file, TESTBENCH_HOME+'/'+kernel))
  #Usage:  lin-analyzer [file.bc] -Ipath=[path] -config=[filename] [kernel-name] -Opath=[path] [options]
  input_bitcode = TESTBENCH_HOME+'/'+kernel+'/'+kernel+'_linked_opt.bc'
  input_path = TESTBENCH_HOME+'/'+kernel
  start_time = time.time()
  os.system('lin-analyzer %s %s %s %s %s' % (input_bitcode, '-Ipath='+input_path, '-config=no', kernel, '-profiling-time'))
  end_time = time.time()
  profiling_times[kernel] = end_time - start_time
  

## Run design space exploration
for kernel in benchmarks:
  counter = 0
  explo_start = time.time()
  input_size = unrolls[kernel][-1]
  for pipeline in pipelines[kernel]:
    if pipeline == upperlevel[kernel]:
      for part in partitions:
        counter += 1
        os.chdir(TESTBENCH_HOME+'/scripts')
        os.system('python run_analyzer.py %s %i %i %s %i' % (kernel, part, 1, pipeline, input_size))
        #print 'kernel name = ' + kernel + ', unroll = 1, ' + 'partition = ' + str(part) + ', pipeline = ' + pipeline + '\n'
    else:
      for unroll in unrolls[kernel]:
        for part in partitions:
          counter += 1
          os.chdir(TESTBENCH_HOME+'/scripts')
          os.system('python run_analyzer.py %s %i %i %s %i' % (kernel, part, unroll, pipeline, input_size))
          #if pipeline == '0':
            #print 'kernel name = ' + kernel + ', unroll = ' + str(unroll) + \
            #', partition = ' + str(part) + ', pipeline = No\n'
          #else:
            #print 'kernel name = ' + kernel + ', unroll = ' + str(unroll) + \
            #', partition = ' + str(part) + ', pipeline = ' + pipeline + '\n'
  designPointNum[kernel] = counter
  explo_end = time.time()
  exploration_times[kernel] = explo_end - explo_start
	

# Write profiling and exploration time into csv file
os.chdir(TESTBENCH_HOME)
RESULT_DIRECTORY = TESTBENCH_HOME+'/scripts/results'
profilingTime_file = RESULT_DIRECTORY + '/profiling_time.csv'
explorationTime_file = RESULT_DIRECTORY + '/exploration_time.csv'
designPointNum_file = RESULT_DIRECTORY + '/design_point_num.csv'
if not os.path.isdir(RESULT_DIRECTORY):
  os.mkdir(RESULT_DIRECTORY)
	
print 'Writing to profiling_time.csv, exploration_time.csv and design_point_num.csv'
prof_time_fp = open(profilingTime_file, 'w')
expl_time_fp = open(explorationTime_file, 'w')
dp_num_fp = open(designPointNum_file, 'w')

print '\n'
for kernel in benchmarks:
  prof_time_fp.write(kernel + ',' + str(profiling_times[kernel])+'\n')
  print 'Profiling time of kernel, ' + kernel + ', ' + str(profiling_times[kernel]) + ' (s)'

prof_time_fp.close()
print '\n'

for kernel in benchmarks:
  expl_time_fp.write(kernel + ',' + str(exploration_times[kernel]) + '\n')
  print 'Exploration time of kernel, ' + kernel + ', ' + str(exploration_times[kernel]) + ' (s)'
	
expl_time_fp.close()
print '\n'
	
for kernel in benchmarks:
  dp_num_fp.write(kernel + ',' + str(designPointNum[kernel]) + '\n')
  print 'Number of design points in kernel, ' + kernel + ', ' + str(designPointNum[kernel])

dp_num_fp.close()
print '\n'

CONFIG2CYCLES_DIR = TESTBENCH_HOME + '/scripts/config2cycles' 
if not os.path.isdir(CONFIG2CYCLES_DIR):
  os.mkdir(CONFIG2CYCLES_DIR)

AREA_DIR = TESTBENCH_HOME + '/scripts/config2area'
if not os.path.isdir(AREA_DIR):
  os.mkdir(AREA_DIR)

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
  fadd_used = []
  fsub_used = []
  fmul_used = []
  fdiv_used = []
  
  KERNEL_HOME = TESTBENCH_HOME + '/' + kernel
  os.chdir(KERNEL_HOME+'/sim')
  dirs = os.listdir(KERNEL_HOME + '/sim/')
 
  config2cycles_filename = CONFIG2CYCLES_DIR + '/' + kernel + '_config2cycles.csv'
  config2cycles = open(config2cycles_filename, 'w')
  config2area_filename = AREA_DIR + '/' + kernel + '_config2area.csv'
  config2area = open(config2area_filename, 'w')
 
  for dir in dirs:
    config_name = dir.split('/')[-1]
    #print config_name
    sum_file_name = dir+'/'+kernel+'_summary.log'
    if os.path.isfile(sum_file_name):
      sum_file = open(sum_file_name, 'r')
      for line in sum_file:
        if 'cycles' in line:
          fpga_cycles.append(int(line.split(' ')[3]))
          config2cycles.write(config_name + ',' + str(line.split(' ')[3]))
        elif 'DSP used' in line:
          dsp_used.append(int(line.split(' ')[2]))
        elif 'BRAM18K used' in line:
          bram18k_used.append(int(line.split(' ')[2]))
        elif 'Fadd used' in line:
          fadd_used.append(int(line.split(' ')[2]))
        elif 'Fsub used' in line:
          fsub_used.append(int(line.split(' ')[2]))
        elif 'Fmul used' in line:
          fmul_used.append(int(line.split(' ')[2]))
        elif 'Fdiv used' in line:
          fdiv_used.append(int(line.split(' ')[2]))
          break
        else:
          continue
      sum_file.close()
  
  config2cycles.close()
  	
  os.chdir(RESULT_DIRECTORY)
  fig = plt.figure()
  fig.suptitle(kernel + ' Execution Time Design Space')
  curr_plot = fig.add_subplot(111)
  #curr_plot.scatter(range(designPointNum[kernel]), fpga_cycles)
  curr_plot.scatter(range(len(fpga_cycles)), fpga_cycles)
  curr_plot.set_xlabel('Configurations')
  curr_plot.set_ylabel('Execution Time (cycles)')
  curr_plot.grid(True)
  plt.savefig(kernel + '-cycles.pdf')
  if len(fpga_cycles) != designPointNum[kernel]:
    raise Exception('Length of fpga_cycles and designPointNum[kernel] mismatched!')
  
  area = [float(dsp)/max_dsp + float(bram18k)/max_bram18k for dsp, bram18k in zip(dsp_used, bram18k_used)]
  fig = plt.figure()
  fig.suptitle(kernel + 'Execution Time - Area Design Space')
  curr_plot = fig.add_subplot(111)
  curr_plot.scatter(area, fpga_cycles)
  curr_plot.set_xlabel('Area')
  curr_plot.set_ylabel('Execution Time (cycles)')
  curr_plot.grid(True)
  plt.savefig(kernel + '-cycles-area.pdf')
  if len(area) != designPointNum[kernel]:
    raise Exception('Lenght of area and designPointNum[kernel] mismatched!')
  
  for dir, area_t in zip(dirs,area):
    config_t = dir.split('/')[-1]
    config2area.write(config_t + ',' + str(area_t) + '\n')
  config2area.close()

