import ast
g2 = []
with open('sizeeffects\\valoriG2.txt', 'U') as f: #apre file
	data=f.read() #copia in una stringa tutto il documento
	#data=data.split('\n') #divide per riga in una lista di stringhe
	#for numero in data[:-1]:
		#print numero
		#energie.append(float(numero))
	g2=ast.literal_eval(data)
EPS = 0.05

risultato = []
for indice,puntoG in enumerate(g2):
	risultato.append(str(EPS * indice) + '\t' + str(puntoG) +'\n')
with open ('gnuPlot.txt', 'w') as f:
	for puntoG in risultato:
		f.write(puntoG)
