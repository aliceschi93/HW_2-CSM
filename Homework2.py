import random
import math
import numpy as np
import time
from multiprocessing import Process



def calcolaL(numeroMolecole):
	L = (2*numeroMolecole)**(1.0/3.0) #L = L/sigma
	return L

def generaCubo(numeroMolecole, L):
	insiemeCoordinate= set()
	while len(insiemeCoordinate)<numeroMolecole:
		x= random.random()*L
		y= random.random()*L
		z= random.random()*L
		insiemeCoordinate.add((x,y,z))
	return list(insiemeCoordinate)
#print generaCubo(2, 13.5)

def controllo(coordinata, L):
	if coordinata < 0:
		return coordinata + L
	if coordinata > L :
		return coordinata - L
	return coordinata

def propostaMossa(molecola, L, delta):
	r1= random.random()
	r2= random.random()
	r3= random.random()
	deltar1= delta*(r1-0.5)
	deltar2= delta*(r2-0.5)
	deltar3= delta*(r3-0.5)
	X=molecola[0]+deltar1
	Y=molecola[1]+deltar2
	Z=molecola[2]+deltar3
	X = controllo(X,L)
	Y = controllo(Y,L)
	Z = controllo(Z,L)
	return (X,Y,Z)

def metropolis(deltaE):
	if deltaE < 0:
		return True
	trovato = False
	p = 0.0
	while not trovato:
		p = random.random()
		if p > 0:
			trovato = True
	if (- math.log(p) > deltaE):
		return True
	return False


def spostaCubo(cubo, delta, L, energia, tabellaEnergie, iterazione, accettazioni, energiaParticella, rc, NMolecole):
	for indice, molecola in enumerate(cubo):
		X,Y,Z = propostaMossa(molecola, L, delta) 
		nuovoCubo = cubo[:]
		nuovoCubo[indice] = (X,Y,Z)
		nuovaEnergia, nuovaTabellaEnergie = aggiornaEnergia(nuovoCubo, L, tabellaEnergie, energia,indice, rc)
		deltaE = nuovaEnergia - energia
		if metropolis(deltaE):
			cubo[indice] = (X,Y,Z)
			energia = nuovaEnergia
			tabellaEnergie = nuovaTabellaEnergie
			if iterazione > 99:
				accettazioni +=1
	energiaParticella.append(energia/NMolecole)
	return (cubo, energia, tabellaEnergie, accettazioni)


def distanza(molecola1, molecola2, L):
	distanzaXYZ=[]
	for indice in range(0,3):
		di= molecola2[indice]-molecola1[indice]
		di-= L*round(di/L)
		distanzaXYZ.append(di)
	sommaQuadrati=0.0
	for coordinata in distanzaXYZ:
		sommaQuadrati+= coordinata**2
	distanza12=math.sqrt(sommaQuadrati)
	return distanza12

def generaEnergie(cubo, L, rc):
	energia = 0.0
	tabellaEnergie = []
	energiaMolecola = []
	for indice, molecola in enumerate(cubo):
		lineaEnergia = []
		for molecola2 in cubo[indice+1:]:
			dist = distanza(molecola, molecola2, L)
			if dist > rc: 
				singolaEnergia = 0
			else:
				singolaEnergia = 0.5*(dist**(-12.0) - dist**(-6.0))
			energia += singolaEnergia
			lineaEnergia.append(singolaEnergia)
		tabellaEnergie.append(lineaEnergia)
	del tabellaEnergie[-1]
	return (energia, tabellaEnergie) 

def aggiornaEnergia(cubo, L, tabellaEnergie, energia, indice, rc):
	nuovaTabellaEnergie = copiaTabella(tabellaEnergie)
	molecolaNuova = cubo[indice]
	for indiMol, molecola in enumerate(cubo):
		if indice == indiMol:
			continue
		dist = distanza(molecolaNuova, molecola, L)
		if dist > rc:
			singolaEnergia = 0
		else:
			singolaEnergia = 0.5 * (dist**(-12.0) - dist**(-6.0))
		if indice < indiMol:
			energia -= tabellaEnergie[indice][indiMol - indice - 1]
			energia += singolaEnergia
			nuovaTabellaEnergie[indice][indiMol - indice - 1] = singolaEnergia
		else:
			energia -= tabellaEnergie[indiMol][indice - indiMol - 1]
			energia += singolaEnergia
			nuovaTabellaEnergie[indiMol][indice - indiMol - 1] = singolaEnergia
	return (energia, nuovaTabellaEnergie)

def energiaTail(rc, NMolecole, L):
	enTail = math.pi * NMolecole * (1.0/L)**(3.0) * (1.0 - 3 * rc**6)/(9 * rc**9)
	return enTail

def pressioneTail(NMolecole, L, rc):
	pressTail = -2 * math.pi * (NMolecole / (L)**(3.0))**2 * (3* rc**(6.0) - 2)/(9 * rc**(9.0))
	return pressTail

def pressioneEx(cubo, NMolecole, L, rc):
	pressEx = 0.0
	for indice, molecola in enumerate(cubo):
		for molecola2 in cubo[indice+1:]:
			dist=distanza(molecola, molecola2, L)
			if dist > rc:
				pressEx += 0.0
			else:
				pressEx += -((1.0/L)**3 * (dist**(-6.0) - 2 * dist**(-12)))
	return pressEx

def copiaTabella(tabellaEnergie):
	nuovaTabellaEnergie = []
	for listaContributi in tabellaEnergie:
		nuovaTabellaEnergie.append(listaContributi[:])
	return nuovaTabellaEnergie

def autocorrelazione(energiaParticella, k, media):
	autoCorr=0.0
	for j in range (0,len(energiaParticella)-k):
		autoCorr+= (energiaParticella[j+k]-media)*(energiaParticella[j]-media)
	C = (1.0/(len(energiaParticella)-k))*autoCorr
	return C

def autocorrAnalisi(energiaParticella, media, delta):
	k = 0
	Trovato = False
	CK = []
	while ((not Trovato) or (k<50)):
		ktemp= autocorrelazione(energiaParticella, k, media)
		CK.append(ktemp)
		if ktemp < 0 and not Trovato:
			k1= k-1
			Trovato = True
		k +=1
	C0 = CK[0]
	tauint = 0.0
	for k in range(1, k1):
		tauint  += CK[k] / C0
	tauint += 0.5 
	sigma = math.sqrt((C0/len(energiaParticella))*2*tauint)
	return(k1, tauint, sigma, CK)

def pairDistr(cubo, NMolecole, L, EPS):
	M = int(L / EPS)
	V = L**3
	g2 = [0] * M
	h = [0] * M
	for indice,molecola in enumerate(cubo[:-2]):
		for molecola2 in cubo[indice+1:-1]:
			dist = distanza(molecola, molecola2, L)
			n = int(dist / EPS)
			Vs = 4.0/3.0 * math.pi * ((n * EPS + EPS)**3 - (n * EPS)**3)
			h[n] = h[n] + 1.0
			g2[n] = (2.0 /(NMolecole * (NMolecole - 1.0))) * (V/Vs) * h[n]
	return g2

def RoleOfDelta(NMolecole, listaDelta, rc):
	L = calcolaL(NMolecole)
	iterazioni = 5000
	for delta in listaDelta:
		box=generaCubo(NMolecole, L)
		energia, tabellaEnergie=generaEnergie(box, L, rc)
		energieDecrescono = []
		accettazioni = 0
		energiaParticella = []
		pressEx = []
		rc = L/2.0
		for indice in range (0,iterazioni):
			box, energia, tabellaEnergie, accettazioni=spostaCubo(box, delta, L, energia, tabellaEnergie, indice, accettazioni, energiaParticella, rc, NMolecole)
			energieDecrescono.append(energia)
			if indice > 99:
				pressEx.append(pressioneEx(box, NMolecole, L, rc))

		with open('roleofdelta\\energiaIterazioni' + str(delta) +'.txt', 'w') as the_file: 
			for ener in energieDecrescono:
				the_file.write(str(ener) + '\n')

		with open('roleofdelta\\pressioneIterazioni' + str(delta) +'.txt', 'w') as the_file: 
			for pressione in pressEx:
				the_file.write(str(pressione) + '\n')

		with open('roleofdelta\\energiaParticella' + str(delta) +'.txt', 'w') as the_file: 
			the_file.write(str(energiaParticella))

		media = 0.0
		for energia in energiaParticella[100:]:
			media += energia

		media /= (len(energiaParticella) - 100)
		k, tau, sigma, CK=autocorrAnalisi(energiaParticella[100:], media, delta)


		with open('roleofdelta\\EnergiaAutocorrAnalisi.' + str(delta) +'.txt', 'w') as the_file:
			the_file.write(str(CK))

		with open('roleofdelta\\EnergiakTauMediaSigma.' + str(delta) +'.txt', 'w') as the_file:
			the_file.write(str(k))
			the_file.write('\n')
			the_file.write(str(tau))
			the_file.write('\n')
			the_file.write(str(media))
			the_file.write('\n')
			the_file.write(str(sigma))

		with open('roleofdelta\\accettazioni' + str(delta) +'.txt', 'w') as the_file: 
			the_file.write(str(accettazioni))
			the_file.write('\n')
			the_file.write(str(NMolecole * (iterazioni - 100)))

		for indice in range(0, len(pressEx)):
			pressEx[indice] += 0.5

		mediaP = 0.0
		for pressione in pressEx:
			mediaP += pressione
		mediaP /= (len(pressEx))

		kP, tauP, sigmaP, CKP=autocorrAnalisi(pressEx, mediaP, delta)

		with open('roleofdelta\\PressioneAutocorrAnalisi.' + str(delta) +'.txt', 'w') as the_file:
			the_file.write(str(CKP))

		with open('roleofdelta\\PressionekTauMediaSigma.' + str(delta) +'.txt', 'w') as the_file:
			the_file.write(str(kP))
			the_file.write('\n')
			the_file.write(str(tauP))
			the_file.write('\n')
			the_file.write(str(mediaP))
			the_file.write('\n')
			the_file.write(str(sigmaP))

		enTail=energiaTail(rc, NMolecole, L)
		pressTail=pressioneTail(NMolecole, L, rc)

		with open('roleofdelta\\valoriEnergiaPressione' + str(delta) +'.txt', 'w') as the_file: 
			the_file.write('Energia media + tail + errore\n')
			the_file.write(str(media) + '\t' + str(enTail) + '\t' + str(sigma))
			the_file.write('\n')
			the_file.write('Pressione media + tail + errore\n')
			the_file.write(str(mediaP) + '\t' + str(pressTail) + '\t' + str(sigmaP))


def RoleOfCuteOff(listarc):
	NMolecole = 60
	iterazioni = 5000
	delta = 1
	L = calcolaL(NMolecole)
	for rc in listarc:
		box=generaCubo(NMolecole, L)
		energia, tabellaEnergie=generaEnergie(box, L, rc)
		energieDecrescono = []
		energiaParticella = []
		pressEx = []
		for indice in range (0,iterazioni):
			box, energia, tabellaEnergie, accettazioni=spostaCubo(box, delta, L, energia, tabellaEnergie, indice, 0, energiaParticella, rc, NMolecole)
			if indice > 99:
				#energieDecrescono.append(energia)
				pressEx.append(pressioneEx(box, NMolecole, L, rc))
		media = 0.0
		for energia in energiaParticella[100:]:
			media += energia
		media /= (len(energiaParticella)-100)


		for indice in range(0, len(pressEx)):
			pressEx[indice] += 0.5

		mediaP = 0.0
		for pressione in pressEx:
			mediaP += pressione
		mediaP /= (len(pressEx))

		enTail=energiaTail(rc, NMolecole, L)
		pressTail=pressioneTail(NMolecole, L, rc)

		with open('roleofcutoff\\valoriEnergiaPressioneCO' + str(rc) +'.txt', 'w') as the_file:
			the_file.write('Energia media + tail energia\n') 
			the_file.write(str(media) + '\t' + str(enTail))
			the_file.write('\n') 
			the_file.write('Pressione media + tail pressione\n') 
			the_file.write(str(mediaP) + '\t' + str(pressTail))

def SizeEffect(delta, listaNMolecole):
	molecoleColG = 200
	for NMolecole in listaNMolecole:
		L = calcolaL(NMolecole)
		box=generaCubo(NMolecole, L)
		EPS = 0.05
		M = int(L / EPS)
		rc= L/2.0
		iterazioni = 5000
		delta = 1.0/2.0
		energia, tabellaEnergie=generaEnergie(box, L, rc)
		energieDecrescono = []
		accettazioni = 0
		energiaParticella = []
		pressEx = []
		gmean = [0] * M
		for indice in range (0,iterazioni):
			box, energia, tabellaEnergie, accettazioni=spostaCubo(box, delta, L, energia, tabellaEnergie, indice, accettazioni, energiaParticella, rc, NMolecole)
			energieDecrescono.append(energia)
			if indice > 99:
				pressEx.append(pressioneEx(box, NMolecole, L, rc))
				if NMolecole == molecoleColG:
					g = pairDistr(box, NMolecole, L, EPS)
					for indice in range(0, len(g)):
						gmean[indice] += g[indice]

		if NMolecole == molecoleColG:
			for indice in range(0,len(gmean)):
				gmean[indice] /= (iterazioni - 100)
			with open('sizeeffects\\valoriG2.txt', 'w') as the_file: 
				the_file.write(str(gmean))

		media = 0.0
		for energia in energiaParticella[100:]:
			media += energia
		media /= (len(energiaParticella)-100)

		for indice in range(0, len(pressEx)):
			pressEx[indice] += 0.5

		mediaP = 0.0
		for pressione in pressEx[100:]:
			mediaP += pressione
		mediaP /= (len(pressEx)-100)

		enTail=energiaTail(rc, NMolecole, L)
		pressTail=pressioneTail(NMolecole, L, rc)

		k, tau, sigma, CK=autocorrAnalisi(energiaParticella[100:], media, delta)


		with open('sizeeffects\\EnergiaAutocorrAnalisi.' + str(NMolecole) +'.txt', 'w') as the_file:
			the_file.write(str(CK))

		with open('sizeeffects\\EnergiakTauMediaSigma.' + str(NMolecole) +'.txt', 'w') as the_file:
			the_file.write(str(k))
			the_file.write('\n')
			the_file.write(str(tau))
			the_file.write('\n')
			the_file.write(str(media))
			the_file.write('\n')
			the_file.write(str(sigma))

		with open('sizeeffects\\valoriEnergiaPressioneSE' + str(NMolecole) +'.txt', 'w') as the_file:
			the_file.write('Energia media (SizeEffect) + tail energia\n') 
			the_file.write(str(media) + '\t' + str(enTail))
			the_file.write('\n') 
			the_file.write('Pressione media (SizeEffect) + tail pressione\n') 
			the_file.write(str(mediaP) + '\t' + str(pressTail))

		kP, tauP, sigmaP, CKP=autocorrAnalisi(pressEx, mediaP, delta)

		with open('sizeeffects\\PressioneAutocorrAnalisi.' + str(NMolecole) +'.txt', 'w') as the_file:
			the_file.write(str(CKP))

		with open('sizeeffects\\PressionekTauMediaSigma.' + str(NMolecole) +'.txt', 'w') as the_file:
			the_file.write(str(kP))
			the_file.write('\n')
			the_file.write(str(tauP))
			the_file.write('\n')
			the_file.write(str(mediaP))
			the_file.write('\n')
			the_file.write(str(sigmaP))

#main

def main():

	try:
		p = Process(target=processoCinghiale, args=( ))
		p2 = Process(target=processoCervo, args=( ))
		p3 = Process(target=processoVolpe, args=())
		p4 = Process(target=processoLepre, args=())
		p.start()
		p2.start()
		p3.start()
		p4.start()
		p.join()
		p2.join()
		p3.join()
		p4.join()
	except Exception as errtxt:
		a=8



def processoCinghiale():
	SizeEffect(1.0, [200])


def processoCervo():
	SizeEffect(1.0, [100])


def processoVolpe():
	NMolecole = 60
	L=calcolaL(NMolecole)
	listaDelta = [1.0/8.0, 1.0 / 4.0, 1.0/ 2.0, 1.0]
	RoleOfDelta(NMolecole, listaDelta, L/2.0)
def processoLepre():
	NMolecole = 60
	L=calcolaL(NMolecole)
	listaDelta = [2.0, L / 2.0]
	RoleOfDelta(NMolecole, listaDelta, L/2.0)
	listarc = [3.0/8.0*L, L / 4.0]
	RoleOfCuteOff(listarc)
	

	#listaDelta = [1.0/8.0, 1.0 / 4.0, 1.0/ 2.0, 1.0, 2.0, L / 2.0]
	#RoleOfDelta(NMolecole, listaDelta, L/2.0)
	#listarc = [3.0/8.0*L, L / 4.0]
	#RoleOfCuteOff(listarc)
	#SizeEffect(1.0)

if __name__ == '__main__':
	main()
