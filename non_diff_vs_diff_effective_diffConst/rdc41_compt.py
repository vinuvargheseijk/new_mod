import moose
import numpy as np
import rdesigneur as rd
import pylab
import matplotlib.pyplot as plt


numDendSegments=100
comptLenBuff=14.4e-06
comptLen=0.36e-06
#comptDia=1e-06
RM=1.0
RA=10.0
CM=0.001

def makeChemyOscillator( name = 'osc', parent = '/library' ):
    model = moose.Neutral( parent + '/' + name )
    compt = moose.CubeMesh( model.path + '/kinetics' )
    """
    This function sets up a simple oscillatory chemical system within
    the script. The reaction system is::

        s ---a---> a  // s goes to a, catalyzed by a.
        s ---a---> b  // s goes to b, catalyzed by a.
        a ---b---> s  // a goes to s, catalyzed by b.
        b -------> s  // b is degraded irreversibly to s.

    in sum, **a** has a positive feedback onto itself and also forms **b**.
    **b** has a negative feedback onto **a**.
    Finally, the diffusion constant for **a** is 1/10 that of **b**.
    """
    # create container for model
    diffConst = 1e-13 # m^2/sec
    motorRate = 1e-6 # m/sec
    concA = 1 # millimolar

    # create molecules and reactions
    a = moose.Pool( compt.path + '/a' )

    c=moose.Pool(compt.path+'/c')
    d=moose.BufPool(compt.path+'/d')
    r2=moose.Reac(compt.path+'/r2')
    moose.connect(r2,'sub',a,'reac')
    moose.connect(r2,'sub',c,'reac')
    moose.connect(r2,'prd',d,'reac')
    r2.Kf=1
    r2.Kb=10


    # Assign parameters
    a.diffConst = diffConst/10
#    b.diffConst = diffConst
#    s.diffConst = diffConst
    c.diffConst = diffConst
    d.diffConst = 0
    return compt

def makeDendProto():
    dend=moose.Neuron('/library/dend')
    prev=rd.buildCompt(dend,'soma',RM=RM,RA=RA,CM=CM,dia=0.3e-06,x=0,dx=comptLenBuff)
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
      
    compt=rd.buildCompt(dend,'dendL',RM=RM,RA=RA,CM=CM,x=x,y=y,dx=comptLenBuff,dy=dy,dia=comptDia)
    moose.connect(prev,'axial',compt,'raxial')

    return dend

moose.Neutral('/library')
makeDendProto()
makeChemyOscillator()

rdes=rd.rdesigneur(
        #chanProto=[['make_HH_Na()','Na'],['make_HH_K()','K']],
        diffusionLength=1e-07,
        turnOffElec=True,
        cellProto=[['elec','dend']],
        chemProto=[['dend','osc']],
        chemDistrib=[['osc','#soma#,#dend#','install','1']],
        #chanDistrib=[
               #['Na','#','Gbar','12000*(dia<1.5e-6)'],
               #['K','#','Gbar','3600*(dia<1.5e-6)']],
        #stimList=[['soma','1','.','inject','(t>0.01&&t<0.2)*1e-10']],
        plotList=[
            ['#', '1', 'dend/a', 'conc', 'conc of A'],
            ['#', '1', 'dend/c', 'conc', 'conc of B']
                 ],
        #moogList=[['#','1','.','Vm','Vm']]
                  )

rdes.buildModel()
av=moose.vec('/model/chem/dend/a')
cv=moose.vec('/model/chem/dend/c')
dv=moose.vec('/model/chem/dend/d')
a=moose.element('/model/chem/dend/a')
d=moose.element('/model/chem/dend/d')
c=moose.element('/model/chem/dend/c')
avec=moose.vec('/model/chem/dend/a').conc
savec=avec.size
print savec
for i in range(0, savec-1,100):
       #c.vec[i].concInit = 0.1
       c.vec[i].concInit = 0.5
       d.vec[i].concInit = 0.5
      # print i
c.vec[60].concInit = 0.5
d.vec[60].concInit = 0.5

moose.reinit()
moose.start(100)
#rdes.displayMoogli(0.00005, 0.05, 0.0)
avec=av.conc
plt.plot(avec)
