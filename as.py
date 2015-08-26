#!/usr/bin/env python

import sys,math,heapq

#parameters
p1 = 5871
p2 = 5871
waterstart = 5876
waterend   = 55750

#inp  = sys.argv[1] # template restart file
#out  = sys.argv[2] # output restart file
#dat  = sys.argv[3] # template plumed file
#plu  = sys.argv[4] # output plumed file
#cut  = int(sys.argv[5]) # radius (angstrom)

#cut2 = cut*cut

#cell = [92.0,70.0,90.0]

def read_template(inp):

    fin = open(inp).readlines()

    nlines = 0
    natoms = 0
    coordinationNumber = 0
    steps = 0

    coord=[]
    coordFlag = False
    coordLine = 1e10

    vel=[]
    velFlag = False
    velLine = 1e10

    coordinationLine = []
    coordinationStopLine = []

    for i,line in enumerate(fin):
        
        if line.find("&COORD")>=0 and i > 100:
            coordFlag = True
            coordLine = i
        elif line.find("&END COORD")>=0 and i > 100 and i < 90000:
            coordFlag = False
            coordStopLine = i

        elif line.find("&VELOCITY")>=0 and i > 55806:
            velFlag = True
            velLine = i
        elif line.find("&END VELOCITY")>=0 and i > 55806:
            velFlag = False
            velStopLine = i
    #	print i
     
        if coordFlag and i > coordLine:
            natoms+=1
    #        print line,i
            atom = line.split()[0]
            x    = float(line.split()[1])
            y    = float(line.split()[2])
            z    = float(line.split()[3])   
        
            coord.append([atom,x,y,z])

        if velFlag and i > velLine:
            
            x    = float(line.split()[0])
            y    = float(line.split()[1])
            z    = float(line.split()[2])	    
        
        vel.append([atom,x,y,z]) 

        nlines+=1

    print coord
    print vel
    #print inp

    print "Number of lines = %d"%(nlines)
    print "Number of atoms = %d"%(natoms)




def write_output(inp,out):
    
    fout = open(out,"w")
    fin = open(inp).readlines()

    for i,line in enumerate(fin):
        
        if line.find("SEED")>=0:
            fout.write("SEED   %d\n")

#        elif i == coordLine+1:
#            for atom in range(0,natoms,1):
#            fout.write("%s %25.16e %25.16e %25.16e\n"%(coord[atom][0],coord[atom][1],coord[atom][2],coord[atom][3]))
                   
#        elif i > coordLine+1 and i < coordStopLine:
#            pass

#        elif i == velLine+1:
#        for atom in range(0,natoms,1):
#            fout.write("  %25.16e %25.16e %25.16e\n"%(vel[atom][0],vel[atom][1],vel[atom][2]))

#        elif i > velLine+1 and i < velStopLine:
#        pass

        else:
            fout.write(line)

    fout.close()
    return

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
	    c3[i].append(CoordinationNumber(h,o,1.2))
    tmp=[]	    
    for c in c3:
        tmp.append(sum(c))
    print tmp
    tmp=heapq.nlargest(2,tmp)

    isADP=0
    if sum(c1)<0.5 and sum(c2)>3.7 and geomean(tmp)>0.7:
        isADP=1

    return [isADP,sum(c1),sum(c2),geomean(tmp)]

def basins(fpdb,bpdb):

    """READ IN ATOM COORDINATES"""
    #forward trajectory
    forward = []
    fin = open(fpdb).readlines()
    gammaP = [float(fin[5874].split()[5]),float(fin[5874].split()[6]),float(fin[5874].split()[7])]

    for i,line in enumerate(fin):
        if i<4 or i==len(fin)-1:
            pass
	else:
  	    atom=[float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]
            if distance(gammaP,atom) < 8:
		forward.append([int(line.split()[1]),line.split()[-1],
			        float(line.split()[5]),
			        float(line.split()[6]),
			        float(line.split()[7])])
    #backward trajectory
    backward = []
    fin = open(bpdb).readlines()
    gammaP = [float(fin[5874].split()[5]),float(fin[5874].split()[6]),float(fin[5874].split()[7])]
    
    for i,line in enumerate(fin):
        if i<4 or i==len(fin)-1:
            pass
        else:
            atom=[float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]
            if distance(gammaP,atom) < 8:
                backward.append([int(line.split()[1]),line.split()[-1],
                                float(line.split()[5]),
                                float(line.split()[6]),
                                float(line.split()[7])])


    """CHECK BASINS"""

    haf=0 
    hbf=0
    hab=0
    hbb=0

    #forward trajectory
    haf=checkATP(forward)
    hbf=checkADP(forward)   
    print hbf
    #backward trajectory	` 
    hab=checkATP(backward)
    hbb=checkADP(backward)
    print hbb

    
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

    return

def reseed(inp,out,seed,coord=False):

    fin  = open(inp).readlines()
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
        basins(forw,back)

    elif cmd=="generate":
	generate(finput,iteration,lsteps,ssteps)

    elif cmd=="seed":
	reseed(finput,foutput,seed,coord)
