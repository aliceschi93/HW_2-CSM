import random
import math
import numpy as np
import time
from multiprocessing import Process


def energiaTail(rc, NMolecole, L):
	enTail = math.pi * NMolecole * (1.0/L)**(3.0) * (1.0 - 3 * rc**6)/(9 * rc**9)
	return enTail

def pressioneTail(NMolecole, L, rc):
	pressTail = -2 * math.pi * (NMolecole / (L)**(3.0))**2 * (3* rc**(6.0) - 2)/(9 * rc**(9.0))
	return pressTail


def calcolaL(numeroMolecole):
	L = (2*numeroMolecole)**(1.0/3.0) #L = L/sigma
	return L

NMolecole = 60
L = calcolaL(NMolecole)
pTailROD=pressioneTail(NMolecole, L, L / 2.0)
pTailROC38=pressioneTail(NMolecole, L, 3*L/8.0) 
pTailROC14=pressioneTail(NMolecole, L, L / 4.0)
L = calcolaL(100)
pTailSE100=pressioneTail(100, L, L / 2.0)
L = calcolaL(200)
pTailSE200=pressioneTail(200, L, L / 2.0)

print pTailROD
print pTailROC38
print pTailROC14
print pTailSE100
print pTailSE200
