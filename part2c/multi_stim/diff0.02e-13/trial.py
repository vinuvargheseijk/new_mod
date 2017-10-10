########################################################################
## This program is part of 'MOOSE', the
## Messaging Object Oriented Simulation Environment.
##           Copyright (C) 2014 Upinder S. Bhalla. and NCBS
## It is made available under the terms of the
## GNU Lesser General Public License version 2.1
## See the file COPYING.LIB for the full notice.
#########################################################################


import math
import numpy
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import moose
import os.path
import xml.etree.ElementTree as ET

def makeModel():
    """
    This example illustrates how to set up a oscillatory Turing pattern 
    in 1-D using reaction diffusion calculations.
    Reaction system is::

        s ---a---> a  // s goes to a, catalyzed by a.
        s ---a---> b  // s goes to b, catalyzed by a.
        a ---b---> s  // a goes to s, catalyzed by b.
        b -------> s  // b is degraded irreversibly to s.
        s ----Receptor---> a  // Receptor is activated by the ligand to cause an increase the conc of a. Here ligand-receptor interaction is not accounted for.
    in sum, **a** has a positive feedback onto itself and also forms **b**.
    **b** has a negative feedback onto **a**.
    Finally, the diffusion constant for **a** is 1/10 that of **b**.

    This chemical system is present in a 1-dimensional (cylindrical) 
    compartment. The entire reaction-diffusion system is set up 
    within the script.
    """
    # create container for model
    r0 = 1e-6	# m
    r1 = 1e-6	# m
    num = 100
    diffLength = 1e-7 # m
    len = num * diffLength	# m
    diffConst = 1e-13 # m^2/sec
    motorRate = 1e-6 # m/sec
    concA = 1 # millimolar
    dt4 = 0.5  # for the diffusion
    dt5 = 0.2   # for the reaction

    model = moose.Neutral( 'model' )
    compartment = moose.CylMesh( '/model/compartment' )
    compartment.r0 = r0
    compartment.r1 = r1
    compartment.x0 = 0
    compartment.x1 = len
    compartment.diffLength = diffLength
    
    assert( compartment.numDiffCompts == num )

    # create molecules and reactions
    a = moose.Pool( '/model/compartment/a' )
    b = moose.Pool( '/model/compartment/b' )
    s = moose.Pool( '/model/compartment/s' )
    e1 = moose.MMenz( '/model/compartment/e1' )
    e2 = moose.MMenz( '/model/compartment/e2' )
    e3 = moose.MMenz( '/model/compartment/e3' )
    e4 = moose.MMenz( '/model/compartment/e4' )
    r1 = moose.Reac( '/model/compartment/r1' )
    rec = moose.Pool( '/model/compartment/rec' )

    moose.connect( e1, 'sub', s, 'reac' )
    moose.connect( e1, 'prd', a, 'reac' )
    moose.connect( a, 'nOut', e1, 'enzDest' )
    e1.Km = 1
    e1.kcat = 1

    moose.connect( e2, 'sub', s, 'reac' )
    moose.connect( e2, 'prd', b, 'reac' )
    moose.connect( a, 'nOut', e2, 'enzDest' )
    e2.Km = 1
    e2.kcat = 0.5

    moose.connect( e3, 'sub', a, 'reac' )
    moose.connect( e3, 'prd', s, 'reac' )
    moose.connect( b, 'nOut', e3, 'enzDest' )
    e3.Km = 0.1
    e3.kcat = 1

    moose.connect( r1, 'sub', b, 'reac' )
    moose.connect( r1, 'prd', s, 'reac' )
    r1.Kf = 0.3 # 1/sec
    r1.Kb = 0. # 1/sec

    moose.connect( e4, 'sub', s, 'reac')
    moose.connect( e4, 'prd', a, 'reac')
    moose.connect( rec, 'nOut', e4, 'enzDest')
    e4.Km = 0.001
    e4.kcat = 4

    # Assign parameters
    a.diffConst = diffConst/10
    b.diffConst = diffConst
    s.diffConst = diffConst
    rec.diffConst = diffConst/50

    # Make solvers
    ksolve = moose.Ksolve( '/model/compartment/ksolve' )
    dsolve = moose.Dsolve( '/model/dsolve' )
    # Set up clocks. The dsolver to know before assigning stoich
    moose.setClock( 4, dt4 )
    moose.setClock( 5, dt5 )
    moose.useClock( 4, '/model/dsolve', 'process' )
    # Ksolve must be scheduled after dsolve.
    moose.useClock( 5, '/model/compartment/ksolve', 'process' )

    stoich = moose.Stoich( '/model/compartment/stoich' )
    stoich.compartment = compartment
    stoich.ksolve = ksolve
    stoich.dsolve = dsolve
    stoich.path = "/model/compartment/##"
    assert( dsolve.numPools == 4 )
    a.vec.concInit = [0.1]*num
#    a.vec[50].concInit *= 1.2 # slight perturbation at one end.
    a.vec[10].concInit *=1.2
    a.vec[35].concInit *=1.2
    a.vec[60].concInit *=1.2
    a.vec[85].concInit *=1.2
    print moose.showfields(a)
    b.vec.concInit = [0.1]*num
    s.vec.concInit = [1]*num
    rec.vec.concInit = [0]*num

def writeXML(storeAvec,fileName,timePer,distS,location):
    if os.path.isfile(fileName):
       tree = ET.parse(fileName)
       root = tree.getroot()
       dist = ET.SubElement(root, 'dist'+str(distS)+'time'+str(timePer)+'location'+str(location))
       print 'in xml'
       for i in range(len(storeAvec)):
         avec = ET.SubElement(dist,'avec'+str(i))
         avec.text = ''.join(str(j)+' ' for j in storeAvec[i])+'\n'
       tree.write(fileName)
    else: 
        root = ET.Element('Data')
        dist = ET.SubElement(root, 'dist'+str(distS)+'time'+str(timePer)+'location'+str(location))
        for i in range(len(storeAvec)):
          avec = ET.SubElement(dist, 'avec'+str(i))
          avec.text = ''.join(str(j)+' ' for j in storeAvec[i])+'\n'
        tree = ET.ElementTree(root)
        tree.write(fileName)

def displayPlots():
    a = moose.element( '/model/compartment/a' )
    b = moose.element( '/model/compartment/b' )
    pos = numpy.arange( 0, a.vec.conc.size, 1 )
    pylab.plot( pos, a.vec.conc, label='a' )
    pylab.plot( pos, b.vec.conc, label='b' )
    pylab.legend()
    pylab.show()

def computeTP(t,distS):
    a = moose.element( '/model/compartment/a' )
    b = moose.element( '/model/compartment/b' )
    s = moose.element( '/model/compartment/s' )
    rec = moose.element( '/model/compartment/rec' )
    w = len(a.vec.conc)
    h = w
    v2d = [[0 for x in range(w)] for y in range(h)]
    anC = []
    mm=moose.vec('/model/compartment/mesh')
    diffL = moose.element('/model/compartment').diffLength
    print a.vec[0].nInit, 'diff length', diffL
    list = [[distS,t]]
    SourceC = 12
    extraDiff = 1e-12
    for j in list:
       anC = []
       anC2 = [] 
       for i in range(len(a.vec.conc)):
            anC.append(2.0*diffL*(SourceC/(math.sqrt(4*math.pi*extraDiff*j[1])))*math.exp(-(mm[i].Coordinates[0]**2)/(4*extraDiff*j[1])))
            for k in range(len(a.vec.conc)):
              anC2.append(4.0*diffL**2*(SourceC/(4*numpy.pi*j[1]*numpy.sqrt(extraDiff*extraDiff)))*numpy.exp(-(mm[i].Coordinates[0]**2)/(4*extraDiff*j[1])-(mm[k].Coordinates[0]**2)/(4*extraDiff*j[1])))
              v2d[i][k] = 4.0*diffL**2*(SourceC/(4*numpy.pi*j[1]*numpy.sqrt(extraDiff*extraDiff)))*numpy.exp(-(mm[i].Coordinates[0]**2)/(4*extraDiff*j[1])-(mm[k].Coordinates[0]**2)/(4*extraDiff*j[1]))
              firstCell = diffL*(SourceC/(math.sqrt(4*math.pi*extraDiff*j[1])))*math.exp(-(mm[0].Coordinates[0]**2)/(4*extraDiff*j[1]))
              firstCell2 = diffL**2*(SourceC/(4*numpy.pi*j[1]*numpy.sqrt(extraDiff*extraDiff)))*numpy.exp(-(mm[0].Coordinates[0]**2)/(4*extraDiff*j[1])-(mm[0].Coordinates[0]**2)/(4*extraDiff*j[1]))
    print 'analaytical sum 1D is',sum(anC)-firstCell, 'max of anC is', max(anC)
    print 'analaytical sum 2D is',sum(anC2)-2*firstCell-firstCell2, 'max of anC2 is', max(anC2), firstCell, firstCell2
      #plt.plot(anC)
      #plt.plot(anC2)
    #raw_input()
    tempRec = []
    anP = True
    if anP == True:
        for i in range(len(a.vec.conc)):
          tempRec.append(v2d[i][distS])
        print 'sum conc at dendrite', sum(tempRec)
    return tempRec

def main(location, distS, timePer):
    runtime = 5000
    displayInterval = 20
    makeModel()
    dsolve = moose.element( '/model/dsolve' )
    moose.reinit()
    #moose.start( runtime ) # Run the model for 10 seconds.

    a = moose.element( '/model/compartment/a' )
    b = moose.element( '/model/compartment/b' )
    s = moose.element( '/model/compartment/s' )
    rec = moose.element( '/model/compartment/rec' )
    print moose.showfields(a)

#    img = mpimg.imread( 'turingPatternTut.png' )
    #imgplot = plt.imshow( img )
    #plt.show()

    #fig = plt.figure( figsize=(12,10) )
    #png = fig.add_subplot(211)
  #  imgplot = plt.imshow( img )
    #ax = fig.add_subplot(212)
    #ax.set_ylim( 0, 2 )
    #plt.ylabel( 'Conc (mM)' )
    #plt.xlabel( 'Position along cylinder (microns)' )
    #pos = numpy.arange( 0, a.vec.conc.size, 1 )
    #line1, = ax.plot( pos, a.vec.conc, label='a' )
    #line2, = ax.plot( pos, b.vec.conc, label='b' )
    #line3, = ax.plot( pos, rec.vec.conc, label='rec')
    #line4, = ax.plot( pos, s.vec.conc, label='s')
    #timeLabel = plt.text(60, 0.4, 'time = 0')
    #plt.legend()
    #fig.canvas.draw()
    #distS = 1
    #posPer = [20,45]
    #timePer = 1
    gaussian = True
    storeAvec = []
    numStimPoints = 2
    print storeAvec
    print rec.vec.conc
    stimPoint = 0
    for ti in range( displayInterval, runtime, displayInterval ):
        if ti>550 and ti<580:
          print stimPoint
          if gaussian == True:
             tempRec = computeTP(timePer,distS)
             #for pos in posPer[stimPoint]:
             rec.vec[location].conc = rec.vec[location].conc+tempRec[0]
             for i in range(location+1,len(a.vec.conc)):
                 rec.vec[i].conc=rec.vec[i].conc+tempRec[i-location]
                 if int(location*2-i)>=0:
                    rec.vec[location*2-i].conc=rec.vec[location*2-i].conc+tempRec[i-location]
             rec.vec.conc = rec.vec.conc/2
             stimPoint = stimPoint+1
          else:
             #rec.vec[10].conc = 0.07
             #rec.vec[30].conc = 0.07
             rec.vec[22].conc = 0.08
             
          print 'receptor sum is', sum(rec.vec.conc), max(rec.vec.conc)
          #plt.figure(10)
          #plt.plot(rec.vec.conc)
          #raw_input()
          
        moose.start( displayInterval )
        storeAvec.append(a.vec.conc)
        #line1.set_ydata( a.vec.conc )
        #line2.set_ydata( b.vec.conc )
        #line3.set_ydata( rec.vec.conc *10 )
        #line4.set_ydata( s.vec.conc)
        #timeLabel.set_text( "time = %d" % ti)
        #fig.canvas.draw()

    print(a.vec.conc)
    file=open('c41','w')
    ty=0
    for tin in a.vec.conc:
        tout=str(tin)
        ty=ty+1
        file.write(str(ty) + ' ' + tout + '\n')

    fileName = 'trial.xml'
    writeXML(storeAvec,fileName,timePer,distS,location)
    return storeAvec



# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
     location = [35,40,45,50,55,60]
     distS = [5,3,1]
     timePer = [1,3,5]
     for i in location:
       for j in distS:
        for k in timePer:
	   main(i,j,k)
           moose.delete('/model')
