import moose
import numpy as np
import rdesigneur as rd
import matplotlib.pyplot as plt


params={
        #'diffusionL':1,
        'diffusionL':1e-07,
        #'diffA':1000,
       'diffA':0, 
       #'diffB':0,
       'diffB':0,
       'dendDia':10e-06,
        #'dendDia':10,
        #'dendL':100,
       'dendL':500e-06,
       }

def makePassiveSoma(name,length,diamater):
	elecid=moose.Neuron('/library/'+name)
        dend=moose.Compartment(elecid.path+'/soma')
        dend.diameter=params['dendDia']
        dend.length=params['dendL']
        dend.x=params['dendL']
        return elecid


def makeChemProto(name='hydra'):
        chem=moose.Neutral('/library/'+name)
        compt=moose.CubeMesh('/library/'+name + '/' + name)
        A=moose.Pool(compt.path+'/A')
        B=moose.Pool(compt.path+'/B')
        #A.diffConst=params['diffA']
        #B.diffConst=params['diffB']
        Adot = moose.Function( A.path + '/Adot' )
        Bdot = moose.Function( B.path + '/Bdot' )
        #Adot.expr="0.001*((x0^2/x1)+1)-1*x0+0.5*x0"
        #Bdot.expr="0*x0+0.001*x0^2-1*x1+1*x1"
        Adot.expr="0.1*((x0^2/x1)+1)-1*x0+0.5*x0"
        Bdot.expr="0*x0+0.1*x0^2-1*x1+1*x1"
 
        
        print "$$$$> ", Adot, Bdot
        print Adot.expr, Bdot.expr
        print moose.showmsg(Adot)
        Adot.x.num = 2 #2
        Bdot.x.num = 2 #2
        #A.nInit=10
        #B.nInit=5
        A.nInit=100
        B.nInit=100
        moose.connect( A, 'nOut', Adot.x[0], 'input' )
        moose.connect( B, 'nOut', Adot.x[1], 'input' )
        moose.connect( Adot, 'valueOut', A, 'increment' )

        moose.connect( A, 'nOut', Bdot.x[0], 'input' )
        moose.connect( B, 'nOut', Bdot.x[1], 'input' )
        moose.connect( Bdot, 'valueOut', B, 'increment' )
        return compt

def main():
    library=moose.Neutral('/library')
    makePassiveSoma( 'cell', params['dendL'], params['dendDia'] )
    makeChemProto()
    rdes=rd.rdesigneur(
        turnOffElec=True,
        chemPlotDt = 0.1, 
        diffusionLength=params['diffusionL'],
        cellProto=[['cell','soma']],
        chemProto=[['hydra','hydra']],
        chemDistrib=[['hydra', 'soma', 'install', '1']],
        plotList=[
            ['soma', '1', 'dend/A', 'n', '# of A'],
            ['soma', '1', 'dend/B', 'n', '# of B']
        ],
    #    moogList = [['soma', '1', 'dend/A', 'n', 'num of A (number)']]
        )
    moose.le( '/library' )
    moose.le( '/library/hydra' )
    moose.showfield( '/library/soma/soma' )
    rdes.buildModel()
    #moose.element('/model/chem/dend/B').vec[50].nInit=15.5
    
    A = moose.vec('/model/chem/dend/A')
    B = moose.vec('/model/chem/dend/B') 
#    A.concInit=1
#    B.concInit=10
    moose.element('/model/chem/dend/A').vec[500].nInit=110
    #moose.element('/model/chem/dend/A').vec[490].nInit=101
    #moose.element('/model/chem/dend/B').vec[25].nInit=0
    #A.nInit=1
    #B.nInit=15
    for i in range(0,499):
        #moose.element('/model/chem/dend/A').vec[i].diffConst=9.45e-12
        moose.element('/model/chem/dend/A').vec[i].diffConst=1e-14
        moose.element('/model/chem/dend/B').vec[i].diffConst=1e-14
        #moose.element('/model/chem/dend/B').vec[i].diffConst=0.27e-09
        
    for i in range(500,999):
        #moose.element('/model/chem/dend/A').vec[i].diffConst=9.45e-12
        moose.element('/model/chem/dend/A').vec[i].diffConst=1e-14
        moose.element('/model/chem/dend/B').vec[i].diffConst=1e-14
        #moose.element('/model/chem/dend/B').vec[i].diffConst=0.27e-09
        
    for i in range(0,200):
        moose.element('/model/chem/dend/A').vec[i].diffConst=1e-14
        
    for i in range(700,999):
        moose.element('/model/chem/dend/A').vec[i].diffConst=1e-14
    #avec = moose.vec( '/model/chem/dend/A' ).n
    #bvec = moose.vec( '/model/chem/dend/B' ).n
        
    print "diffconst is ",moose.element('/model/chem/dend/A').vec[50].diffConst
    #plt.ion()
    #fig=plt.figure(figsize=(12,10))
    #ax=fig.add_subplot(111)
    #ax.set_ylim( 0, 5)
    #ax1=fig.add_subplot(212)
    #plt.xlabel('position')
    #plt.ylabel('# of A')
    #timeLabel = plt.text(1, 1, 'time=0.0')
    #pos=np.arange(0,A.n.size,1)
    #line1,=ax.plot(pos,A.n,label='A')
    #line2,=ax.plot(pos,B.n,label='B')
    #fig.canvas.draw()
    #displayInterval=2
    #runTime=1000
    moose.reinit()
    moose.start(2000)
    avec = moose.vec( '/model/chem/dend/A' ).n
    bvec = moose.vec( '/model/chem/dend/B' ).n
    plt.plot(avec)
    
    #moose.start(10)    
    #plt.plot(avec)
    #plt.show()
    
    #print "Before for loop"
    #for t in range(displayInterval,runTime,displayInterval):
        #print t, "entering loop"
        #moose.start(displayInterval)
        #timeLabel.set_text( "Time = %f" % t)
        #line1.set_ydata(A.n)
        #line2.set_ydata(B.n)
        #plt.plot(avec)
        #plt.show()
        #fig.canvas.draw()

        #plt.plot(avec)

    #plt.ioff()'''
        
    #for i in range(98,99):
        #moose.element('/model/chem/dend/A').vec[i].diffConst=1e-06
        #moose.element('/model/chem/dend/B').vec[i].diffConst=0
    

    #Adot = moose.element('/model/chem/dend/A/Adot')
    #Bdot = moose.element('/model/chem/dend/B/Bdot')
   
    #print "\n\n\t\t Before Run", Adot, Bdot, A.nInit, B.nInit, A.n, B.n, A.concInit, B.concInit
    
    #moose.reinit()
    #moose.start(2000)
    
    #print "\n\n\t\t After Run", A.nInit, B.nInit, A.n, B.n, A.concInit, B.concInit
    
   # rdes.display()
   # rdes.displayMoogli( 1, 400, 0.001 )
    #avec = moose.vec( '/model/chem/dend/A' ).n
    #bvec = moose.vec( '/model/chem/dend/B' ).n
    #tab = moose.element( '/model/graphs/plot0')
    #dt = tab.dt
    #tvec=[]
    #tvec.append( tab.vector )
   # xplot = []
   # xplot.append( avec )
  ##  t = np.arange( 0, len( tvec[0] ), 1.0 ) * dt
    #plt.plot(bvec)
    
    #print bvec, len(bvec), bvec
    
    return bvec,avec

if __name__ == '__main__':
    bvec,avec = main()
     
 

                   
