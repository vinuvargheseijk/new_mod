#########################################################################
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
import xml.etree.ElementTree as ET
import os.path

def makeModel():
    """
    This example illustrates how to set up a oscillatory Turing pattern 
    in 1-D using reaction diffusion calculations.
    Reaction system is::

        s ---a---> a  // s goes to a, catalyzed by a.
        s ---a---> b  // s goes to b, catalyzed by a.
        a ---b---> s  // a goes to s, catalyzed by b.
        b -------> s  // b is degraded irreversibly to s.

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
    diffLength = 1e-6 # m
    len = num * diffLength	# m
    diffConst = 1e-12 # m^2/sec
    motorRate = 1e-6 # m/sec
    concA = 1 # millimolar
    dt4 = 0.002  # for the diffusion
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
    r1 = moose.Reac( '/model/compartment/r1' )
    c=moose.Pool('/model/compartment/c')
    d=moose.Pool('/model/compartment/d')
    r2=moose.Reac('/model/compartment/r2')
    moose.connect(r2,'sub',a,'reac')
    moose.connect(r2,'sub',c,'reac')
    moose.connect(r2,'prd',d,'reac')
    r2.Kf=0.0
    r2.Kb=0.0
    #moose.connect( e1, 'sub', s, 'reac' )
    #moose.connect( e1, 'prd', a, 'reac' )
    #moose.connect( a, 'nOut', e1, 'enzDest' )
    #e1.Km = 1
    #e1.kcat = 1

    #moose.connect( e2, 'sub', s, 'reac' )
    #moose.connect( e2, 'prd', b, 'reac' )
    #moose.connect( a, 'nOut', e2, 'enzDest' )
    #e2.Km = 1
    #e2.kcat = 0.5

    #moose.connect( e3, 'sub', a, 'reac' )
    #moose.connect( e3, 'prd', s, 'reac' )
    #moose.connect( b, 'nOut', e3, 'enzDest' )
    #e3.Km = 0.1
    #e3.kcat = 1

    #moose.connect( r1, 'sub', b, 'reac' )
    #moose.connect( r1, 'prd', s, 'reac' )
    #r1.Kf = 0.3 # 1/sec
    #r1.Kb = 0. # 1/sec

    # Assign parameters
    a.diffConst = diffConst
    b.diffConst = 0
    s.diffConst = 0
    c.diffConst = 0
    d.diffConst = 0
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
    assert( dsolve.numPools == 5 )
    #a.vec.concInit = [0.1]*num
#    a.vec[50].concInit *= 1.2 # slight perturbation at one end.
    #a.vec[60].concInit *= 1.2
    #a.vec[40].concInit *= 1.2
    a.vec[0].nInit = 100
    #a.vec[30].concInit *= 1.2
    #a.vec[50].concInit *= 1.2
    #b.vec[99].concInit *= 0.5
    #for i in range(0, num-1,50):
       #c.vec[i].concInit = 0.1
       #c.vec[i].concInit = 0.5
       #d.vec[i].concInit = 0.5
      # print i
    c.vec[0].nInit = 0
    d.vec[0].nInit = 0
    #c.vec[30].concInit = 0.5
    #d.vec[30].concInit = 0.5

def writeXML(timeArr,halfWidth,maxFWHH,trialNum,fileName):
    if os.path.isfile(fileName):
       tree = ET.parse(fileName)
       root = tree.getroot()
       widths = ET.SubElement(root, 'widths'+trialNum)
       widths.text = ''.join(str(i)+' ' for i in halfWidth)+'\n'
       times = ET.SubElement(root, 'times')
       times.text = ''.join(str(j)+' ' for j in timeArr)+'\n'
       xmaxFWHH = ET.SubElement(root, 'xmaxFWHH'+trialNum)
       xmaxFWHH.text = ''.join(str(k)+' ' for k in maxFWHH)+'\n'
       #tree = ET.ElementTree(root)
       tree.write('data.xml') 
    else:
       root = ET.Element('Data')
       widths = ET.SubElement(root, 'widths'+trialNum)
       widths.text = ''.join(str(i)+' ' for i in halfWidth)+'\n'
       times = ET.SubElement(root, 'times')
       times.text = ''.join(str(j)+' ' for j in timeArr)+'\n'
       xmaxFWHH = ET.SubElement(root, 'xmaxFWHH'+trialNum)
       xmaxFWHH.text = ''.join(str(k)+' ' for k in maxFWHH)+'\n'
       tree = ET.ElementTree(root)
       tree.write('data.xml') 


def displayPlots():
    a = moose.element( '/model/compartment/a' )
    b = moose.element( '/model/compartment/b' )
    pos = numpy.arange( 0, a.vec.conc.size, 1 )
    pylab.plot( pos, a.vec.conc, label='a' )
    pylab.plot( pos, b.vec.conc, label='b' )
    pylab.legend()
    pylab.show()

def main():
    global num
    num = 100
    runtime = 11
    makeModel()
    dsolve = moose.element( '/model/dsolve' )
    moose.reinit()
    #moose.start( runtime ) # Run the model for 10 seconds.

    a = moose.element( '/model/compartment/a' )
    b = moose.element( '/model/compartment/b' )
    c = moose.element( '/model/compartment/c' )
    d = moose.element( '/model/compartment/d' )
    s = moose.element( '/model/compartment/s' )
    halfWidth = []
    timeArr = []
    maxFWHH = []
    r2=moose.element('/model/compartment/r2')
    trialNum = str(c.diffConst)+str(r2.Kf)+str(r2.Kb)
    plt.ion()
    
    #fig = plt.figure( figsize=(12,10) )
    #png = fig.add_subplot(211)
  #  imgplot = plt.imshow( img )
    #ax = fig.add_subplot(212)
    #ax.set_ylim( 0, 2 )
    #plt.ylabel( 'Conc (mM)' )
    #plt.xlabel( 'Position along cylinder (microns)' )
    #pos = numpy.arange( 0, a.vec.conc.size, 1 )
    #timeLabel = plt.text(60, 0.4, 'time = 0')
    #line1, = ax.plot( pos, a.vec.conc, label='a' )
    #line2, = ax.plot( pos, b.vec.conc, label='c' )
    #line3, = ax.plot( pos, b.vec.conc, label='d' )
    #plt.legend()
    #fig.canvas.draw()
    anC = []
    mm=moose.vec('/model/compartment/mesh')
    diffL = moose.element('/model/compartment').diffLength
    print a.vec[0].nInit
    list = [[1e-12,1], [2e-12,2], [4e-12,4], [6e-12,6], [8e-12,8], [10e-12,10]]
    for j in list:
      anC = [] 
      for i in range(len(a.vec.conc)):
        #anC.append((10/numpy.sqrt(4*numpy.pi*a.diffConst*1))*numpy.exp(-(i**2)/4*a.diffConst*1))
        #print mm[i].Coordinates[0]
        anC.append(diffL*(100/(math.sqrt(4*math.pi*a.diffConst*4)))*math.exp(-(mm[i].Coordinates[0]**2)/(4*a.diffConst*4)))
        #anC.append(4*diffL**2*(3/(4*numpy.pi*j[1]*numpy.sqrt(a.diffConst*a.diffConst)))*numpy.exp(-(mm[i].Coordinates[0]**2)/(4*a.diffConst*j[1])-j[0]/(4*a.diffConst*j[1])))
        # anC.append(4*diffL**2*(3/(4*numpy.pi*4*numpy.sqrt(4*numpy.pi*4*a.diffConst*a.diffConst*a.diffConst)))*numpy.exp(-(mm[i].Coordinates[0]**2)/(4*a.diffConst*4)-(mm[i].Coordinates[0]**2)/(4*a.diffConst*4)-j/(4*a.diffConst*4)))
        #print (3/(numpy.sqrt(4*numpy.pi*a.diffConst*1)))*numpy.exp(-(mm[i].Coordinates[0]**2)/(4*a.diffConst*1))
      print sum(anC)
      plt.plot(anC)
    raw_input()
    time = 0    
    for t in range( 0, runtime):
        moose.start(0.5 )
        time = time + 0.5
        if int(time == 4):
          print sum(a.vec.n)
          temp = a.vec.n
          plt.figure(2)
          plt.plot(temp)
          raw_input()
        aList = a.vec.conc
        maxa = max(aList)
        #indMax = numpy.where(aList==maxa)
        #indMax = num/2
        indMax = 0
        halfMax = maxa/2
        indHalfMax = (numpy.abs(aList[0:num]-halfMax)).argmin()
        #print 'maximum value is', maxa, 'postion is', numpy.where(aList==maxa)
        print 'maximum value is', maxa
        print 'half height is at', indHalfMax
        print 'half maximum is', aList[indHalfMax]
        #halfWidth.append(int(indMax[0])-indHalfMax)
        #halfWidth.append(indMax-indHalfMax) 
        halfWidth.append(indHalfMax-indMax)
        timeArr.append(t)
        #line1.set_ydata( a.vec.conc )
        #line2.set_ydata( b.vec.conc )
        #line2.set_ydata( c.vec.conc )
        #line3.set_ydata( d.vec.conc )
        #timeLabel.set_text( "time = %d" % t )
        #fig.canvas.draw()
        print 'minimum and maximum of a', min(a.vec.conc), max(a.vec.conc)

    maxFWHH.append(max(halfWidth))       
    print 'max FWHH is', max(halfWidth)
    fileName = 'data.xml'
    writeXML(timeArr,halfWidth,maxFWHH,trialNum,fileName)
    #treeName = ET.parse('data.xml')
    
    #plotXML(treeName)

# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
	main()
