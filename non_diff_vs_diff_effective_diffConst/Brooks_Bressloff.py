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
       #'dendL':100e-06,
       }
       
numDendSegments=100
comptLenBuff=14.4e-06
comptLen=0.216e-06
#comptDia=1e-06
RM=1.0
RA=10.0
CM=0.001

      
def makeDendProto():
    dend=moose.Neuron('/library/dend')
    prev=rd.buildCompt(dend,'soma',RM=RM,RA=RA,CM=CM,dia=0.3e-06,x=0,dx=comptLen)
    x=comptLenBuff
    y=0.0
    comptDia=0.3e-06

    for i in range(numDendSegments):
      dx=comptLen
      dy=0
      #comptDia +=1.7e-08
      compt=rd.buildCompt(dend,'dend'+str(i),RM=RM,RA=RA,CM=CM,x=x,y=y,dx=dx,dy=dy,dia=comptDia)
      moose.connect(prev,'axial',compt,'raxial')
      prev=compt
      x+=dx
      y+=dy
      
    #compt=rd.buildCompt(dend,'dendL',RM=RM,RA=RA,CM=CM,x=x,y=y,dx=comptLenBuff,dy=dy,dia=comptDia)
    #moose.connect(prev,'axial',compt,'raxial')

    return dend


def makeChemProto(name='hydra'):
        chem=moose.Neutral('/library/'+name)
        compt=moose.CubeMesh('/library/'+name + '/' + name)
        A=moose.Pool(compt.path+'/A')
        #B=moose.Pool(compt.path+'/B')
        C=moose.Pool(compt.path+'/C')
        r1=moose.Reac(compt.path+'/r1')
        #e1=moose.MMenz(compt.path+'/e1')
        moose.connect(r1,'sub',A,'reac')
        moose.connect(r1,'prd',C,'reac')
        #moose.connect(C,'nOut',e1,'enzDest')
        r1.Kf=1
        r1.Kb=0.1
        A.nInit=0
        #B.nInit=0
        C.nInit=0
        return compt

def main():
    library=moose.Neutral('/library')
    #makePassiveSoma( 'cell', params['dendL'], params['dendDia'] )
    makeDendProto()
    makeChemProto()
    rdes=rd.rdesigneur(
        turnOffElec=True,
        chemPlotDt = 0.1, 
        diffusionLength=params['diffusionL'],
        #cellProto=[['cell','soma']],
        cellProto=[['elec','dend']],
        chemProto=[['hydra','hydra']],
        chemDistrib=[['hydra', '#soma#,#dend#', 'install', '1']],
        plotList=[
            ['soma', '1', 'dend/A', 'n', '# of A']
            #['soma', '1', 'dend/B', 'n', '# of B']
        ],
    #    moogList = [['soma', '1', 'dend/A', 'n', 'num of A (number)']]
        )
    moose.le( '/library' )
    moose.le( '/library/hydra' )
    #moose.showfield( '/library/soma/soma' )
    rdes.buildModel()
    #moose.element('/model/chem/dend/B').vec[50].nInit=15.5
    
    A = moose.element('/model/chem/dend/A')
    #B = moose.element('/model/chem/dend/B')
    C = moose.element('/model/chem/dend/C')
    A.diffConst=1e-12
    #B.diffConst=0
    C.diffConst=0
 

    #moose.element('/model/chem/dend/A').vec[200].nInit=1.2
    #moose.element('/model/chem/dend/A').vec[100].nInit=1.2
    avec=moose.vec('/model/chem/dend/A').n
    savec=avec.size
    randpera=np.random.uniform(1,2,savec)
    randperc=np.random.uniform(1,2,savec)
    #print randper, randper.size
    for j in range(0,savec-1,10):
        moose.element('/model/chem/dend/C').vec[j].nInit=1
    
    for i in range(0,savec-1,1):
        moose.element('/model/chem/dend/A').vec[i].nInit=randpera[i]
        #moose.element('/model/chem/dend/C').vec[i].nInit=randperc[i]
        #print moose.element('/model/chem/dend/A').vec[i].nInit

    moose.reinit()
    #moose.start(2)
    avec=moose.vec('/model/chem/dend/A').n
    storeAvec=[]
    #storeBvec=[]
    storeCvec=[]
    for t in range(0,4000,100):
        moose.start(100)
        avec=moose.vec('/model/chem/dend/A').n
        #bvec=moose.vec('/model/chem/dend/B').n
        cvec=moose.vec('/model/chem/dend/C').n
        storeAvec.append(avec)
        #storeBvec.append(bvec)
        storeCvec.append(cvec)
        
        #plt.plot(avec)
        
    #avec = moose.vec( '/model/chem/dend/A' ).n
    #bvec = moose.vec( '/model/chem/dend/B' ).n
    #cvec = moose.vec( '/model/chem/dend/C' ).n
    #plt.plot(bvec) 
    

        

    
    return storeAvec,storeCvec

if __name__ == '__main__':
    storeAvec,storeCvec = main()
     
 

                   
