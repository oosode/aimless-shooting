#!/usr/bin/env python

import sys,math,heapq

#parameters
p1 = 5871
p2 = 5871
waterstart = 5876
waterend   = 55750


def read_input(org,pdb,out):

  p1 = 5871
  p2 = 5871
  waterstart = 5876
  waterend   = 55750

  cut  = 8
  cut2 = cut*cut

  cell = [92.0,70.0,90.0]

  nlines = 0
  natoms = 0
  coordinationNumber = 0
  steps = 0

  coord=[]
  coordFlag = False
  coordLine = 1e10

  coordinationLine = []
  coordinationStopLine = []

  forg = open(org).readlines()
  
  for i,line in enumerate(forg):
    
    if line.find("&QM_KIND Oqm")>=0:
        oqmStartLine = i
    elif line.find("&QM_KIND Nqm")>=0:
        oqmStopLine = i-1
    elif line.find("&QM_KIND Hqm")>=0:
	hqmStartLine = i
    elif line.find("&QM_KIND Pqm")>=0:
	hqmStopLine = i-1
 
    nlines+=1

   
  fpdb = open(pdb).readlines()
 
  for i,line in enumerate(fpdb):
  
    if i>0 and line.find("ATOM")>=0 and i<len(fpdb)-1:

        natoms+=1

        atom = line[11:17]
	x    = float(line[31:38])
	y    = float(line[39:46])
	z    = float(line[47:54])
	
	coord.append([atom,x,y,z])

  nCloseWats = 0
  closeWat = []

  for atom in range(waterstart-1,waterend-1,3):

    dist1 = 0.0
    dist2 = 0.0

    for k in range(0,3,1):
	tempdist = coord[p1-1][k+1]-coord[atom][k+1]
        if tempdist < -cell[k]/2.0:
	    tempdist += cell[k]
        elif tempdist > cell[k]/2.0:
            tempdist -= cell[k]

        dist1 += tempdist*tempdist

        tempdist = coord[p2-1][k+1]-coord[atom][k+1]
        if tempdist < -cell[k]/2.0:
            tempdist += cell[k]
        elif tempdist > cell[k]/2.0:
            tempdist -= cell[k]

        dist2 += tempdist*tempdist

    if dist1 < cut2 or dist2 < cut2:
        #make sure water is in box
        for k in range(0,3,1):
            if coord[atom][k+1] < -cell[k]/2.0:
		coord[atom+0][k+1] += cell[k]
		coord[atom+1][k+1] += cell[k]
                coord[atom+2][k+1] += cell[k]      
	    elif coord[atom][k+1] > cell[k]/2.0:
                coord[atom+0][k+1] -= cell[k]
                coord[atom+1][k+1] -= cell[k]
                coord[atom+2][k+1] -= cell[k]

    
        closeWat.append(atom+1)
        #print "%d"%(atom),
        nCloseWats+=1      
  
  fout = open(out,"w")
  
  for i,line in enumerate(forg):
    
    if i == oqmStartLine+3:
        fout.write("       MM_INDEX")
        for j in range(0,nCloseWats,1):
            fout.write(" %d"%(closeWat[j]))
	    if (j+1)%9==0 and (j+1)!=nCloseWats:
	        fout.write(" \\\n       ")
        fout.write("\n")
    elif i > oqmStartLine+3 and i < oqmStopLine:
        pass
    elif i == hqmStartLine+4:
        fout.write("       MM_INDEX")
	for j in range(0,nCloseWats,1):
	    fout.write(" %d %d"%(closeWat[j]+1,closeWat[j]+2))
	    if (j+1)%5==0 and (j+1)!=nCloseWats:
	        fout.write(" \\\n       ")
        
        fout.write("\n")
    elif i > hqmStartLine+4 and i < hqmStopLine:
        pass

    else:
        fout.write(line)

  fout.close()


def distance(a,b):
    
    x = a[0]-b[0]
    y = a[1]-b[1]
    z = a[2]-b[2]

    d = x*x + y*y + z*z

    return math.sqrt(d)

def geomean(vals):
    
    n=1.0
    for val in vals:
	n *= val

    return math.pow(n,1.0/float(len(vals)))

def CoordinationNumber(a,b,r0):

    r = distance(a,b)
    #from plumed coordination number
    num = 1 - pow(r/r0,6)
    den = 1 - pow(r/r0,12)

    return num/den

def checkATP(atoms):

    # gamma phosphorous
    Pgamma=[[i for i in atoms if i[0] == 5871][0][2], \
            [i for i in atoms if i[0] == 5871][0][3], \
            [i for i in atoms if i[0] == 5871][0][4]]

    Obeta=[] # beta oxygens
    for n in [5868,5869,5870]:
        Obeta.append([[i for i in atoms if i[0] == n][0][2], \
                      [i for i in atoms if i[0] == n][0][3], \
                      [i for i in atoms if i[0] == n][0][4]])
    
    Owater=[] # all other oxygens
    for atom in atoms:
	if atom[1]=="O" and atom[0]!=5868 and atom[0]!=5869 and atom[0]!=5870:
            Owater.append([atom[2],atom[3],atom[4]])

    c1=[]
    for o in Obeta:
        c1.append(CoordinationNumber(o,Pgamma,2.38))

    c2=[] 
    for o in Owater:
        c2.append(CoordinationNumber(o,Pgamma,2.38))
 #   print heapq.nlargest(5,c2)

    isATP=0
    if sum(c1)>0.5 and sum(c2)<3.7:
        isATP=1

    return [isATP,sum(c1),sum(c2)] 

def checkADP(atoms):
    # gamma phosphorous
    Pgamma=[[i for i in atoms if i[0] == 5871][0][2], \
            [i for i in atoms if i[0] == 5871][0][3], \
            [i for i in atoms if i[0] == 5871][0][4]]

    Obeta=[] # beta oxygens
    for n in [5868,5869,5870]:
        Obeta.append([[i for i in atoms if i[0] == n][0][2], \
                      [i for i in atoms if i[0] == n][0][3], \
                      [i for i in atoms if i[0] == n][0][4]])

    Owater=[] # all other oxygens
    for atom in atoms:
        if atom[1]=="O" and atom[0]!=5868 and atom[0]!=5869 and atom[0]!=5870:
            Owater.append([atom[2],atom[3],atom[4]])

    Hwater=[] # all hydrogens
    for atom in atoms:
        if atom[1]=="H":
	    Hwater.append([atom[2],atom[3],atom[4]])

    # find oxygen closest to the gamma phosphorous 
    # excluding the gamma oxygens for calculation
    # of the coordination number of gamma oxygens
    minimum=[0,1e10] # minimum distance [id,distance]
    for atom in atoms:
	if atom[1]=="O" and atom[0]!=5872 and atom[0]!=5873 and atom[0]!=5874:
            tmp=distance(Pgamma,[atom[2],atom[3],atom[4]])
	    if tmp<minimum[1]:
	        minimum[1]=tmp
	        minimum[0]=atom[0]

    Ogamma=[]
    for n in [5872,5873,5874,minimum[0]]:
        Ogamma.append([[i for i in atoms if i[0] == n][0][2], \
                       [i for i in atoms if i[0] == n][0][3], \
                       [i for i in atoms if i[0] == n][0][4]])
             
    # determine parameters 

    # coordination number between gamma phosphorous
    # and beta oxygens.
    c1=[]
    for o in Obeta:
        c1.append(CoordinationNumber(o,Pgamma,2.38))

    # coordination number between gamma phosphorous
    # and other oxygens.
    c2=[]
    for o in Owater:
        c2.append(CoordinationNumber(o,Pgamma,2.38))
    
    # coordination number of gamma oxygens and all
    # other hydrogens
    c3=[[]for i in range(4)]
    for i,o in enumerate(Ogamma):
        for h in Hwater:
	    c3[i].append(CoordinationNumber(h,o,1.50))

    tmp=[]	    
    for c in c3:
        tmp.append(sum(c))
    #print tmp
    tmp=heapq.nlargest(2,tmp)

    isADP=0
    if sum(c1)<0.5 and sum(c2)>3.7: # and geomean(tmp)>0.75:
        isADP=1

#    return [isADP,sum(c1),sum(c2),geomean(tmp)]
    return [isADP,sum(c1),sum(c2)]

def basins_xyz(fpdb,bpdb):
    
    """READ IN ATOM COORDINATES"""
    #forward trajectory
    forward = []
    fin = open(fpdb).readlines()
    gammaP = [float(fin[5872].split()[1]),float(fin[5872].split()[2]),float(fin[5872].split()[3])]
    
    for i,line in enumerate(fin):
        if i<2:
            pass
        else:
            atom=[float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
            if distance(gammaP,atom) < 8:
                forward.append([i-1,line.split()[0],
                                float(line.split()[1]),
                                float(line.split()[2]),
                                float(line.split()[3])])
    #backward trajectory
    backward = []
    fin = open(bpdb).readlines()
    gammaP = [float(fin[5872].split()[1]),float(fin[5872].split()[2]),float(fin[5872].split()[3])]
    
    for i,line in enumerate(fin):
        if i<2:
            pass
        else:
            atom=[float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
            if distance(gammaP,atom) < 8:
                backward.append([i-1,line.split()[0],
                                 float(line.split()[1]),
                                 float(line.split()[2]),
                                 float(line.split()[3])])


    """CHECK BASINS"""
    
    haf=0
    hbf=0
    hab=0
    hbb=0
    
    #forward trajectory
    haf=checkATP(forward)
    hbf=checkADP(forward)

    #backward trajectory
    hab=checkATP(backward)
    hbb=checkADP(backward)

    print "HAF", haf
    print "HBF", hbf
    print "HAB", hab
    print "HBB", hbb
    
    basin=open("basin_evals.txt","a")
    # This test is different than in the master h.f script from reference.
    # I'm thinking it is a typo there, but that this should be correct.
    conclusive = (haf[0] + hbf[0]) * (hab[0] + hbb[0])
    
    if conclusive == 0:
        basin.write("Inconclusive %d %d %d %d\n"%(haf[0],hbf[0],hab[0],hbb[0]))
    elif conclusive == 2:
        basin.write("Overlapping %d %d %d %d\n"%(haf[0],hbf[0],hab[0],hbb[0]))
    elif conclusive == 1:
        basin.write("Conclusive %d %d %d %d\n"%(haf[0],hbf[0],hab[0],hbb[0]))
    else:
        print "That's wierd. Exiting..."
        exit(0)

    basin.close()
    
    af=open("haf","w")
    ab=open("hab","w")
    bf=open("hbf","w")
    bb=open("hbb","w")
    
    af.write("%d"%(haf[0]))
    ab.write("%d"%(hab[0]))
    bf.write("%d"%(hbf[0]))
    bb.write("%d"%(hbb[0]))
    
    af.close()
    ab.close()
    bf.close()
    bb.close()
    
    return

def generate(inp,iteration,lsteps,ssteps):

    fin  = open(inp).readlines()
    
    fout = open("forw.org","w") # forward 
    bout = open("back.org","w") # backword 
    dout = open("dt.org","w")   # dt

    velFlag   = False
    printFlag = False
    dumpFlag  = False
    trajFlag  = False
  
    for i,line in enumerate(fin):
        if line.find("STEPS")>=0 and i < 50:
            fout.write("     STEPS %d\n"%(lsteps))
	    bout.write("     STEPS %d\n"%(lsteps))
	    dout.write("     STEPS %d\n"%(ssteps))
	elif line.find("PROJECT_NAME")>=0:
	    fout.write("   PROJECT_NAME FORWARD.%d\n"%(iteration))
            bout.write("   PROJECT_NAME BACKWARD.%d\n"%(iteration))
	    dout.write("   PROJECT_NAME DT.%d\n"%(iteration))	
	elif line.find("&VELOCITY")>=0 and i > 100:
	    velFlag=True
            fout.write(line)
            bout.write(line)
            dout.write(line)
	elif line.find("&END VELOCITY")>=0 and i > 100:
	    velFlag=False
            fout.write(line)
            bout.write(line)
            dout.write(line)
        elif line.find("&PRINT")>=0 and i < 200:
            printFlag=True
            fout.write(line)
            bout.write(line)
            dout.write(line)        
        elif line.find("&END PRINT")>=0 and i <200:
            printFlag=False
            fout.write(line)
            bout.write(line)
            dout.write(line)
        elif line.find("&DUMP_PDB")>=0:
            dumpFlag=True
            fout.write(line)
            bout.write(line)
            dout.write(line)
        elif line.find("&END DUMP_PDB")>=0:
            dumpFlag=False
            fout.write(line)
            bout.write(line)
            dout.write(line)
        elif line.find("&TRAJECTORY")>=0:
            trajFlag=True
            fout.write(line)
            bout.write(line)
            dout.write(line)
        elif line.find("&END TRAJECTORY")>=0:
            trajFlag=False
            fout.write(line)
            bout.write(line)
            dout.write(line)
	else:
	    if velFlag:
		x=float(line.split()[0])
		y=float(line.split()[1])
		z=float(line.split()[2])
	        fout.write("  %25.16E %25.16E %25.16E\n"%(x,y,z))
		bout.write("  %25.16E %25.16E %25.16E\n"%(-x,-y,-z))
		dout.write("  %25.16E %25.16E %25.16E\n"%(x,y,z))
            elif printFlag:
                if line.find("MD")>=0:
                    if trajFlag:
                        fout.write("         MD  1\n")
                        bout.write("         MD  1\n")
                        dout.write("         MD  1\n")
                    else:
                        fout.write("         MD  %d\n"%(lsteps))
                        bout.write("         MD  %d\n"%(lsteps))
                        dout.write("         MD  %d\n"%(ssteps))
                else:
                    fout.write(line)
                    bout.write(line)
                    dout.write(line)
            elif dumpFlag:
                if line.find("MD")>=0:
#                    fout.write("         MD  %d\n"%(lsteps))
#                    bout.write("         MD  %d\n"%(lsteps))
#                    dout.write("         MD  %d\n"%(ssteps))
                    fout.write("         MD  1\n")
                    bout.write("         MD  1\n")
                    dout.write("         MD  1\n")
                else:
                    fout.write(line)
                    bout.write(line)
                    dout.write(line)
	    else:
 		fout.write(line)
                bout.write(line)
                dout.write(line)

    fout.close()
    bout.close()
    dout.close()

    fout=open("forw.vmd","w")
    bout=open("back.vmd","w")
    dout=open("dt.vmd","w")

    fin=open("input/xyz2pdb.vmd").readlines() 

    for i,line in enumerate(fin):
        if line.find("addfile")>=0:
	    fout.write("mol addfile FORWARD.xyz type xyz first -1 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n")
	    bout.write("mol addfile BACKWARD.xyz type xyz first -1 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n")
	    dout.write("mol addfile DT.xyz type xyz first -1 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n")
	elif line.find("writepdb")>=0:
	    fout.write("$sel writepdb %d.forw.coor\n"%(iteration))
	    bout.write("$sel writepdb %d.back.coor\n"%(iteration))
	    dout.write("$sel writepdb %d.dt.coor\n"%(iteration))
	else:
	    fout.write(line)
	    bout.write(line)
	    dout.write(line)

    fout.close()
    bout.close()
    dout.close()

    return

def reseed(inp,out,seed,coord=False):

    read_input(inp,coord,out)

    fin  = open(out).readlines()
    fout = open(out,"w")

    for i,line in enumerate(fin):
	if line.find("SEED")>=0:
            fout.write("   SEED  %04d\n"%(seed))
	elif line.find("COORD_FILE_NAME")>=0:
	    if coord:
  	        fout.write("       COORD_FILE_NAME %s\n"%(coord))
	    else:
	        fout.write(line)
        else:
	    fout.write(line)
		
    return


if __name__=="__main__":

#    basins("SETUP-pos-1.pdb","SETUP-pos-1.pdb")
#    exit(0)
   
    cmd=""
    finput=""
    foutput=""
    coord=False
    vel=False
    seed=0000
    iteration=0
    
    if len(sys.argv)<2:

        print("\nusage: python as.py -s -g -c -p -v -i 'input' -o 'output'\n")

	print("-c -- Count the convergence into the wells.")
	print("-g -- Generate aimless shooting throws (forw,back,dt).") # include iteration number
        print("-s -- Seed simulation.") # include seed number 

        print("-i -- Input file.")
	print("-o -- Output file.")
        print("-p -- PDB file.")
	print("-v -- Velocity file.")
 
        print 
	exit(0)

    for i,val in enumerate(sys.argv):
	if val=="-i":
	    finput=sys.argv[i+1]
	elif val=="-o":
	    foutput=sys.argv[i+1]
	elif val=="-p":
	    coord=sys.argv[i+1]
	elif val=="-v":
	    vel=sys.argv[i+1]
	elif val=="-c":
	    cmd="count"
	    forw=sys.argv[i+1]
	    back=sys.argv[i+2]
	elif val=="-g":
	    cmd="generate"
	    iteration=int(sys.argv[i+1])
            lsteps=int(sys.argv[i+2])
            ssteps=int(sys.argv[i+3])
	elif val=="-s":
	    cmd="seed"
	    seed=int(sys.argv[i+1])

    if cmd=="count":
	basins_xyz(forw,back)

    elif cmd=="generate":
	generate(finput,iteration,lsteps,ssteps)

    elif cmd=="seed":
	reseed(finput,foutput,seed,coord)
