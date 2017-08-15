import moose
import numpy as np
import rdesigneur as rd
import matplotlib.pyplot as plt


params={
        'diffusionL':1e-06,
        'diffA':1e-06,
        'diffB':1e-06,
        'dendL':100e-06,
        'dendDia':10e-06,
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
        A.diffConst=params['diffA']
        B.diffConst=params['diffB']
        Adot = moose.Function( A.path + '/Adot' )
        Bdot = moose.Function( B.path + '/Bdot' )
        Adot.expr="-0.28*x0*x1"
        Bdot.expr="0.28*x0*x1-0.03*x0"
        
        Adot.x.num = 2 #2
        Bdot.x.num = 2
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
            ['soma', '1', 'dend/B', 'n', '# of B'],
        ],
    #    moogList = [['soma', '1', 'dend/A', 'n', 'num of A (number)']]
        )
    moose.le( '/library' )
    moose.le( '/library/hydra' )
    moose.showfield( '/library/soma/soma' )
    rdes.buildModel()
    #moose.element('/model/chem/dend/B').vec[50].nInit=30.5
    moose.element('/model/chem/dend/A').vec[50].nInit=10000
    
    A = moose.vec('/model/chem/dend/A')
    B = moose.vec('/model/chem/dend/B')
#    A.concInit=1
#    B.concInit=10
    #A.nInit=1
    #B.nInit=30
    A.nInit=10000
    B.nInit=10000
    
    

    #Adot = moose.element('/model/chem/dend/A/Adot')
    #Bdot = moose.element('/model/chem/dend/B/Bdot')
   
    #print "\n\n\t\t Before Run", Adot, Bdot, A.nInit, B.nInit, A.n, B.n, A.concInit, B.concInit
    
    moose.reinit()
    moose.start(500)
    
    #print "\n\n\t\t After Run", A.nInit, B.nInit, A.n, B.n, A.concInit, B.concInit
    
   # rdes.display()
   # rdes.displayMoogli( 1, 400, 0.001 )
    avec = moose.vec( '/model/chem/dend/A' ).n
    bvec = moose.vec( '/model/chem/dend/B' ).n
    #tab = moose.element( '/model/graphs/plot0')
    #dt = tab.dt
    #tvec=[]
    #tvec.append( tab.vector )
   # xplot = []
   # xplot.append( avec )
  ##  t = np.arange( 0, len( tvec[0] ), 1.0 ) * dt
    plt.plot(bvec)
    
    #print bvec, len(bvec), bvec
    
    return avec

if __name__ == '__main__':
    main()
     
 

                   
