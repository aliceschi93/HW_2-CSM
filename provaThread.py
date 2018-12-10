from multiprocessing import Process

def myfunction(mystring):
    for ele in mystring:
    	print ele



def main():
	try:
		p = Process(target=myfunction, args=(['MyStringHereimimimimimimimimimim',1,2,3,4,5,6,7],))
		#t = Thread(None,myfunction,None,(['MyStringHereimimimimimimimimimim',1,2,3,4,5,6,7],1))
		#t2 = Thread(None,myfunction,None,(['asjkjskjsksjksjbbbbbbbbbmmmmmmmmmmmmmddddddddddddddtrgrgrg','qqqq','opopopo','1pomn!'],1))
		p2 = Process(target=myfunction, args=(['asjkjskjsksjksjbbbbbbbbbmmmmmmmmmmmmmddddddddddddddtrgrgrg','qqqq','opopopo','1pomn!'],))
		p2.start()
		p.start()
		p.join()
		p2.join()
	except Exception as errtxt:
		print errtxt



if __name__ == '__main__':
	main()

