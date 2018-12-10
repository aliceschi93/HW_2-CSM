# libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ast
 
correlazione = []
with open('autocorrAnalisi.0.125.txt', 'U') as f: #apre file
	data=f.read() #copia in una stringa tutto il documento
	#data=data.split('\n') #divide per riga in una lista di stringhe
	#for numero in data[:-1]:
		#print numero
		#energie.append(float(numero))
	correlazione=ast.literal_eval(data)


# Create a dataset:
#df=pd.DataFrame({'k': range(1,101), 'sigma': np.random.randn(100)*15+range(1,101) })
df=pd.DataFrame({'k': range(0,700), 'C(k)': np.array(correlazione)})
# plot
plt.plot('k', 'C(k)', data=df, linewidth=0.3)
plt.xlabel('Value of k')
plt.ylabel('Value of C(k)')
plt.show()