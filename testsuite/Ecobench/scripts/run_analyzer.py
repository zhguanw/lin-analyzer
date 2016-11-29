  #!/usr/bin/env python
import sys
import os
import os.path
import shutil

def gen_bitcode(directory, kernel):
  if not 'LLVM_BIN_HOME' in os.environ:
	  raise Exception('Set LLVM_BIN_HOME directory as an environment variable')

  KERNEL_HOME = directory + '/' + kernel
  os.chdir(KERNEL_HOME)
  clang_arg = ' -g -O1 -emit-llvm -c '
  llvm_link_arg = ' main.bc abstract_kernel.bc ' + kernel + '.bc -o ' + kernel + '_linked.bc'
  opt_arg = ' -mem2reg -instnamer -lcssa -indvars ' + kernel + '_linked.bc -o ' + kernel + '_linked_opt.bc'
  os.system('clang %s' % (clang_arg + 'main.c -o main.bc'))
  os.system('clang %s' % (clang_arg + 'abstract_kernel.c -o abstract_kernel.bc'))
  os.system('clang %s' % (clang_arg + kernel + '.c -o ' + kernel + '.bc'))
  os.system('llvm-link %s' % (llvm_link_arg))
  os.system('opt %s' % (opt_arg))

def gen_config(directory, kernel, part, unroll, pipe, input_size):

  print '--Running config.main()'

  d = 'p%s_u%s_P%s' % (part, unroll, pipe)

  print 'Kernel = %s, Part = %s, unroll = %s' % (kernel, part, unroll)
  
  array_names = {
  'atax' : ['A','x','y','tmp'],
  'bicg' : ['A','r','s','p','q'],
  'convolution2d' : ['A','B'],
  'convolution3d' : ['A','B'],
  'gemm' : ['A','B','C'],
  'gesummv' : ['A','B','x','y','tmp'],
  'mm' : ['A','B','C'],
  'mvt' : ['a','x1','y1'],
  'syr2k' : ['A','B','C'],
  'syrk' : ['A','C'],
  }
  array_partition_type = {
  'atax' : ['cyclic','complete','complete','complete'],
  'bicg' : ['cyclic','complete','complete','complete','complete'],
  'convolution2d' : ['cyclic','cyclic'],
  'convolution3d' : ['cyclic','cyclic'],
  'gemm' : ['cyclic','cyclic','cyclic'],
  'gesummv' : ['cyclic','cyclic','complete','complete','complete'],
  'mm' : ['cyclic','cyclic','cyclic'],
  'mvt' : ['cyclic','complete','complete'],
  'syr2k' : ['cyclic','cyclic','cyclic'],
  'syrk' : ['cyclic','cyclic'],
  }
  array_size = {
  'atax_128' : ['16384','128','128','128'],
  'atax_64' : ['4096','64','64','64'],
  'atax_32' : ['1024','32','32','32'],
  'atax_16' : ['256','16','16','16'],
  'bicg_256' : ['65536','256','256','256','256'],
  'bicg_128' : ['16384','128','128','128','128'],
  'bicg_64' : ['4096','64','64','64','64'],
  'bicg_32' : ['1024','32','32','32','32'],
  'bicg_16' : ['256','16','16','16','16'],
  'convolution2d_126' : ['16384', '16384'],
  'convolution3d_30' : ['32768','32768'],
  'gemm_128' : ['16384','16384','16384'],
  'gemm_64' : ['4096','4096','4096'],
  'gemm_32' : ['1024','1024','1024'],
  'gemm_16' : ['256','256','256'],
  'gesummv_128' : ['16384','16384','128','128','128'],
  'gesummv_64' : ['4096','4096','64','64','64'],
  'gesummv_32' : ['1024','1024','32','32','32'],
  'gesummv_16' : ['256','256','16','16','16'],
  'mm_128' : ['16384','16384','16384'],
  'mm_64' : ['4096','4096','4096'],
  'mm_32' : ['1024','1024','1024'],
  'mm_16' : ['256','256','256'],
  'mvt_256' : ['65536','256','256'],
  'mvt_128' : ['16384','128','128'],
  'mvt_64' : ['4096','64','64'],
  'mvt_32' : ['1024','32','32'],
  'mvt_16' : ['256','16','16'],
  'syr2k_128' : ['16384','16384','16384'],
  'syr2k_64' : ['4096','4096','4096'],
  'syr2k_32' : ['1024','1024','1024'],
  'syr2k_16' : ['256','256','256'],
  'syrk_128' : ['16384','16384'],
  'syrk_64' : ['4096','4096'],
  'syrk_32' : ['1024','1024'],
  'syrk_16' : ['256','256'],
  }
  #wordsize in bytes
  #sizeof(float) = 4
  array_wordsize = {
  'atax' : ['4','4','4','4'],
  'bicg' : ['4','4','4','4','4'],
  'convolution2d' : ['4','4'],
  'convolution3d' : ['4','4'],
  'gemm' : ['4','4','4'],
  'gesummv' : ['4','4','4','4','4'],
  'mm' : ['4','4','4'],
  'mvt' : ['4','4','4'],
  'syr2k' : ['4','4','4'],
  'syrk' : ['4','4'],
  }
  
  BaseFile = directory
  os.chdir(BaseFile)

  if not os.path.isdir(BaseFile + '/sim/'):
    os.mkdir(BaseFile + '/sim/')

  if os.path.isdir(BaseFile + '/sim/' + d):
    shutil.rmtree(BaseFile + '/sim/' + d)

  if not os.path.isdir(BaseFile + '/sim/' + d):
    os.mkdir(BaseFile + '/sim/' + d)

  print 'Writing config file'
  config = open(BaseFile + '/sim/' + d + '/config_' + d, 'w')
  
  if pipe != '0':
    config.write('pipeline,' + kernel + ',0,' + str(pipe) + "\n")

  #memory partition
  names = array_names[kernel]
  types = array_partition_type[kernel]
  sizes = array_size[kernel+'_'+str(input_size)]
  wordsizes = array_wordsize[kernel]
  assert (len(names) == len(types) and len(names) == len(sizes))
  for name,type_t,size,wordsize in zip(names, types, sizes, wordsizes):
    config.write('array,' + name + ',' + str(int(size)*int(wordsize)) + \
    ',' + str(wordsize) + "\n")
    if type_t == 'complete':
      config.write('partition,'+ type_t + ',' + name + ',' + \
      str(int(size)*int(wordsize)) + "\n")
    elif type_t == 'block' or type_t == 'cyclic':
      config.write('partition,'+ type_t + ',' + name + ',' + \
      str(int(size)*int(wordsize)) + ',' + str(wordsize) + ',' + str(part) + "\n")
    else:
      print "Unknown partition type: " + type_t
      sys.exit(0)

  #config.close()

  #loop unrolling
  if kernel == 'atax':
    config.write('unrolling,atax,0,1,7,1\n')
    config.write('unrolling,atax,0,2,9,%s\n' %(unroll))

  elif kernel == 'bicg':
    config.write('unrolling,bicg,0,1,11,1\n')
    config.write('unrolling,bicg,0,2,13,%s\n' %(unroll))

  elif kernel == 'convolution2d':
    config.write('unrolling,convolution2d,0,1,17,1\n')
    config.write('unrolling,convolution2d,0,2,18,%s\n' %(unroll))

  elif kernel == 'convolution3d':
    config.write('unrolling,convolution3d,0,1,16,1\n')
    config.write('unrolling,convolution3d,0,2,17,1\n')
    config.write('unrolling,convolution3d,0,3,18,%s\n' %(unroll))

  elif kernel == 'gemm':
    config.write('unrolling,gemm,0,1,7,1\n')
    config.write('unrolling,gemm,0,2,8,1\n')
    config.write('unrolling,gemm,0,3,11,%s\n' %(unroll))

  elif kernel == 'gesummv':
    config.write('unrolling,gesummv,0,1,7,1\n')
    config.write('unrolling,gesummv,0,2,10,%s\n' %(unroll))

  elif kernel == 'mm':
    config.write('unrolling,mm,0,1,7,1\n')
    config.write('unrolling,mm,0,2,8,1\n')
    config.write('unrolling,mm,0,3,10,%s\n' %(unroll))

  elif kernel == 'mvt':
    config.write('unrolling,mvt,0,1,7,1\n')
    config.write('unrolling,mvt,0,2,8,%s\n' %(unroll))
	
  elif kernel == 'syr2k':
    config.write('unrolling,syr2k,0,1,13,1\n')
    config.write('unrolling,syr2k,0,2,14,1\n')
    config.write('unrolling,syr2k,0,3,16,%s\n' %(unroll))

  elif kernel == 'syrk':
    config.write('unrolling,syrk,0,1,14,1\n')
    config.write('unrolling,syrk,0,2,15,1\n')
    config.write('unrolling,syrk,0,3,17,%s\n' %(unroll))

  config.close()

def main(kernel, part, unroll, pipe, size):

  if not 'TESTBENCH_HOME' in os.environ:
    raise Exception('Set TESTBENCH_HOME directory as an environment variable')

  TESTBENCH_HOME = os.getenv('TESTBENCH_HOME')
  KERNEL_HOME = TESTBENCH_HOME + '/' + kernel

  os.chdir(KERNEL_HOME)
  d = 'p%s_u%s_P%s' % (part, unroll, pipe)

  #Generate accelerator design config file
  gen_config(KERNEL_HOME, kernel, part, unroll, pipe, size)

  print 'Start Lin-Analyzer'
  config_file = 'config_' + d

	# Trace input store in KERNEL_HOME folder
  input_path = KERNEL_HOME
  output_path = os.path.join(KERNEL_HOME, 'sim', d)
  print 'Changing directory to %s' % output_path
  os.chdir(output_path)
	
  #Usage:  lin-analyzer [file.bc] -Ipath=[path] -config=[filename] [kernel-name] -Opath=[path] [options]
  os.system('lin-analyzer %s %s %s %s %s %s %s' % (input_path+'/'+kernel+'_linked_opt.bc', '-Ipath='+input_path, \
	'-config='+config_file, kernel, '-Opath='+output_path, '-no-trace', '-verbose'))
  #os.system('lin-analyzer %s %s %s %s %s %s' % (input_path+'/'+kernel+'_linked_opt.bc', '-Ipath='+input_path, \
	#'-config='+config_file, kernel, '-Opath='+output_path, '-no-trace'))
  
  #delete *.gz files to save disk space
  os.system('rm %s' % ('*.gz'))

if __name__ == '__main__':
  kernel = sys.argv[1]
  part = sys.argv[2]
  unroll = sys.argv[3]
  pipe = sys.argv[4]
  size = sys.argv[5]
  main(kernel, part, unroll, pipe, size)
