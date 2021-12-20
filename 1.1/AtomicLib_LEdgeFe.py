def AtomicLib_LEdgeFe(Config):
	AP = {'Gd':{'SOC3d':0,'Fdd2':0,'Fdd4':0}, 'Ex':{'SOC2p':0, 'SOC3d':0,'Fdd2':0,'Fdd4':0,'Fpd':0,'Gpd1':0,'Gpd3':0}, 'Shift':0}
	if Config==7: #Fe+
		AP['Gd']['SOC3d']=0.046
		AP['Gd']['Fdd2']=9.762
		AP['Gd']['Fdd4']=6.018

		AP['Ex']['SOC2p']=8.202
		AP['Ex']['SOC3d']=0.059
		AP['Ex']['Fdd2']=10.623
		AP['Ex']['Fdd4']=6.560
		AP['Ex']['Fpd']=6.143
		AP['Ex']['Gpd1']=4.467
		AP['Ex']['Gpd3']=2.538

		AP['Shift'] = 709.7

	elif Config==6: #Fe2+
		AP['Gd']['SOC3d']=0.052
		AP['Gd']['Fdd2']=10.966
		AP['Gd']['Fdd4']=6.815

		AP['Ex']['SOC2p']=8.200
		AP['Ex']['SOC3d']=0.067
		AP['Ex']['Fdd2']=11.779
		AP['Ex']['Fdd4']=7.327
		AP['Ex']['Fpd']=6.793
		AP['Ex']['Gpd1']=5.004
		AP['Ex']['Gpd3']=2.844

		AP['Shift'] = 711.2

	elif Config==5: #Fe3+
		AP['Gd']['SOC3d']=0.059
		AP['Gd']['Fdd2']=12.043
		AP['Gd']['Fdd4']=7.535

		AP['Ex']['SOC2p']=8.199
		AP['Ex']['SOC3d']=0.074
		AP['Ex']['Fdd2']=12.8189
		AP['Ex']['Fdd4']=8.023
		AP['Ex']['Fpd']=7.446
		AP['Ex']['Gpd1']=5.566
		AP['Ex']['Gpd3']=3.166

		AP['Shift'] = 712.7
	elif Config==4: #Fe4+
		AP['Gd']['SOC3d']=0.066
		AP['Gd']['Fdd2']=13.030
		AP['Gd']['Fdd4']=8.198

		AP['Ex']['SOC2p']=8.199
		AP['Ex']['SOC3d']=0.082
		AP['Ex']['Fdd2']=13.776
		AP['Ex']['Fdd4']=8.668
		AP['Ex']['Fpd']=8.103
		AP['Ex']['Gpd1']=6.153
		AP['Ex']['Gpd3']=3.503

		AP['Shift'] = 714.2

	return AP