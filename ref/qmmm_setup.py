#!/usr/bin/env python

import time
def qmkind(f,foutput):
        fin=open(f).readlines()
        qmkinds=[]

        for i,line in enumerate(fin):
		if line.find("&QM_KIND") >= 0:
			q=line.split()[1]
			vals=[]			

			for k in range(10):
	                        if fin[i+k].find("&END QM_KIND")>=0:
        		                break

				else:
					for j in range(len(fin[i+k].split())):
						try:
							vals.append(int(fin[i+k].split()[j]))
						except ValueError:
							print "Oops! That's not an integer"

			qmkinds.append([q,vals])

	fout=open(foutput,"a")
	for l in range(len(qmkinds)):
		for m in range(len(qmkinds[l][1])):
			if m%5==0:
			    fout.write("QMKIND %s\t"%(qmkinds[l][0]))
	
			fout.write(" %6s"%(str(qmkinds[l][1][m]).ljust(6)))

			if m%5==4:
			    if m==(len(qmkinds[l][1])-1):
				continue

                            fout.write("\n")

		fout.write("\n\n")

	fout.close()				

def mmkind(f,foutput):
	fin=open(f).readlines()
	mmkinds=[]

	for i,line in enumerate(fin):
		if line.find("&MM_KIND") >= 0:
			m=line.split()[1]
			r=fin[i+1].split()[1]

			tmp=[m,r]
			mmkinds.append(tmp)

	fout=open(foutput,"a")
	for j in range(len(mmkinds)):
		fout.write("MMKIND %s\t %s\n"%(mmkinds[j][0],mmkinds[j][1].rjust(7)))

	fout.write("\n")
	fout.close()

def link(f,foutput):
	fin=open(f).readlines()
	links=[]
	
	for i,line in enumerate(fin):
		if line.find("&LINK") >= 0:
			qindex=kind=mindex=link=alpha=""
			for k in range(10):
				if fin[i+k].find("&END LINK")>=0:
					break

				elif fin[i+k].find("QM_INDEX")>=0:
					qindex=fin[i+k].split()[1]
				elif fin[i+k].find("QM_KIND")>=0:
					kind=fin[i+k].split()[1]
				elif fin[i+k].find("MM_INDEX")>=0:
					mindex=fin[i+k].split()[1]
				elif fin[i+k].find("LINK_TYPE")>=0:
					link=fin[i+k].split()[1]
				elif fin[i+k].find("ALPHA_IMOMM")>=0:
					alpha=fin[i+k].split()[1]

			tmp=[qindex,kind,mindex,link,alpha]
			links.append(tmp)

	fout=open(foutput,"a")
	for j in range(len(links)):
		fout.write("LINK %s\t %s\t %s\t %s\t %s\n"%(links[j][0],links[j][1],links[j][2],links[j][3],links[j][4]))

	fout.write("\n")
        fout.close()

def g3x3(f,foutput):
        fin=open(f).readlines()
        g3s=[]

        for i,line in enumerate(fin):
		if line.find("&G3X3") >= 0:
			mol=ex=""
			atoms=[]
			distances=[]
		
			for k in range(10):
				if fin[i+k].find("&END G3X3")>=0:
					break
				
				elif fin[i+k].find("MOLNAME")>=0:
					mol=fin[i+k].split()[1]
                                elif fin[i+k].find("EXCLUDE_QM")>=0:        
                                        ex=fin[i+k].split()[1]		 
				elif fin[i+k].find("ATOMS")>=0:
					for l in range(1,len(fin[i+k].split())):
						atoms.append(fin[i+k].split()[l])
                                elif fin[i+k].find("DISTANCES")>=0:
                                        for l in range(1,len(fin[i+k].split())):
                                                distances.append(fin[i+k].split()[l])				

			tmp=[mol,atoms,distances,ex]
			g3s.append(tmp)

	fout=open(foutput,"a")
	for j in range(len(g3s)):
		fout.write("G3X3 %s  \t %s\n"%(g3s[j][0],g3s[j][3]))

	fout.write("\n")
	fout.close()
				
if __name__=="__main__":
	
	foutput="qmmm."+str(time.time())+".init"

	finput="path.1.2.0.0.constr.inp"
	qmkind(finput,foutput)	

	finput="l.init"
	link(finput,foutput)

	finput="same.txt"
	mmkind(finput,foutput)
	
	finput="path.1.2.0.0.constr.inp"
	g3x3(finput,foutput)
