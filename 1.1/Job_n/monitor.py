import os

def monitor_print(texto,directorio_actual= os.path.dirname(os.path.realpath(__file__)),new=False):
	monitor_output = open(directorio_actual+'/real_time_monitor.txt','a')
	if new:
		monitor_output.write("----------------------------------Nueva sesión de cálculos (inicia programa)------------------------------"+'\n')
	print(texto)
	monitor_output.write(texto+'\n')
	monitor_output.close()