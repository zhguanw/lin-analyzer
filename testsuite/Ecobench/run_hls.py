import os
import sys
import commands
import subprocess

root = os.getcwd()
polybench = ['atax', 'bicg', 'convolution2d', 'convolution3d', 'correlation', 'covariance', 'fdtd2d', 'gemm', 'gesummv', 'gramschmidt', 'mm2', 'mm3', 'mvt', 'syr2k', 'syrk']
#polybench = ['atax']
reports = "report_repository"
report_dir = root + '/' + reports
subprocess.call(["mkdir", "-p", reports])

for name in polybench:
	print 'Working on kernel ' + name
	kernel_path = root + '/' + name
	#print kernel_path
	os.chdir(kernel_path)
	subprocess.call(["vivado_hls", "-f", "script.tcl"])
	print 'Copy synthesized report into root folder'
	syn_rpt_name = kernel_path+'/'+name+'/solution1/syn/report/'+name+'_csynth.rpt'
	subprocess.call(["cp", syn_rpt_name, report_dir])
	print 'Clean Vivado_HLS project for kernel ' + name
	subprocess.call(["rm", "-r", name])
	os.chdir(root)
	print 'Kernel ' + name + ' Done!'

print '\nParsing synthesized reports to obtain latency and area of each kernel'
os.chdir(report_dir)
temp_filenames = commands.getstatusoutput("ls .")
filenames = temp_filenames[1].split('\n')
#print temp_filenames
result_dir = root+'/results'
subprocess.call(["mkdir", "-p", result_dir])
clock_check_str = "    |  Clock  | Target| Estimated| Uncertainty|\n"
latency_check_str = "    |   Latency   |   Interval  | Pipeline|\n"
resources_check_str = "|       Name      | BRAM_18K| DSP48E|   FF   |  LUT  |\n"

ecobenchResult = open(result_dir+'/ecobenchResult.csv', 'w')
title_str = "Benchmark,Target Clock,Estimated Clock,Uncertainty,min Latency,max Latency,min Interval,max Interval,BRAM_18K,DSP48E,FF,LUT"
ecobenchResult.write(title_str + '\n')
print title_str
for filename in filenames:
	f = open(report_dir + '/' + filename,'r')
	
	clock_str_enc = False
	clock_ct = 0	

	lat_str_enc = False
	lat_ct = 0
	lat_found = False

	res_str_enc = False
	res_ct = 0

	out = ''
	
	kernel = filename.replace('_', ' ').split()
	kernel_name = kernel[0]
	for line in f:
		y = line.replace('|',' ').split()

		if(len(y)>=4 and y[0] == "Clock" and y[1] == "Target" and y[2] == "Estimated" and y[3] == "Uncertainty"):
			clock_str_enc = True
		if(clock_str_enc == True):
			clock_ct += 1
		if(clock_ct == 3):
			clock_ct = 0
			clock_str_enc = False
			x = line.replace('|', ' ').split()
			out += ',' + str(x[1]) + ',' + str(x[2]) + ',' + str(x[3])

		if(lat_found == False and len(y)>=3 and y[0] == "Latency" and y[1] == "Interval" and y[2] == "Pipeline"):
			lat_str_enc = True
		if(lat_str_enc == True):
			lat_ct += 1
		if(lat_ct == 4):
			lat_ct = 0
			lat_str_enc = False
			lat_found = True

			x = line.replace('|',' ').split()
			out += ',' + str(long(x[0])) + ',' + str(long(x[1])) + ',' + str(long(x[2])) + ',' + str(long(x[3]))
			#out += ';' + str( (long(x[1]) + long(x[0]) )/2.0)
			#out += ';' + str( (long(x[3]) + long(x[2]) )/2.0)

		if(line == resources_check_str):
			res_str_enc = True
		if(res_str_enc == True):
			res_ct += 1
		if(res_ct == 11):
			res_ct = 0
			res_str_enc = False

			x = line.replace('|',' ').split()
			out += ',' +  x[1] + ',' + x[2] + ',' + x[3] + ',' + x[4]
			break
	f.close()

	print kernel_name + out
	ecobenchResult.write(kernel_name + out + '\n')

ecobenchResult.close()
os.chdir(root)
#print root
print "Finish collecting ecobench FPGA results!\n"
