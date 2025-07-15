import matplotlib.pyplot as plt
import lssd as ld
import os

os.mkdir('outputs') if os.path.exists('outputs')==False else None
file='datatest.binx'

data=ld.LumiData(file) #read the .binx
mdf = data.maindataframe #maindataframe
data.integrate_OSL([1,2],dt=0.1) #integrates between 1 and 2 seconds

for i in range(1,4):
	proc = data.maindataframe.loc[data.maindataframe.id==i] #select the register on dataframe
	data.plot(i, **{'color':'purple'}) #plot the signal
	plt.savefig(f'outputs/id{i}_sample{proc.sample_name.item()}.png')
	plt.close()
	
	
#Updating sample names by the carousel position
#Sample A, B, C, D, E
#Each sample has 3 aliquots
sample_dict = {'A':[1,2,3],'B':[4,5,6],'C':[7,8,9],'D':[10,11,12],'E':[13,14,15]}

for sample in list(sample_dict.keys()):
    mdf.loc[mdf.carousel_pos.between(sample_dict[sample][0],sample_dict[sample][-1]),'sample_name'] = sample   
    print(f'Aliquots from sample {sample}: {sample_dict[sample]}')
    
#Acessing the signals
signal = data.measures[5]
plt.plot(signal[0], signal[1], '-r')
print(len(signal[1]))
plt.show()
