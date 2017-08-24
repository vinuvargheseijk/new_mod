import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import numpy as np
 
def plotXML(tree,trials):
    xaxis = tree.find('times')
    plt.ion()
    fig = plt.figure( figsize=(6,6) )
    xValues = [float(j) for j in xaxis.text.split()]
    cList = ['blue','green','red', 'black','yellow','cyan']
    cCount = 0
    yValues1=[]
    diffMax = []
    for i in trials:
      yaxis = tree.find('widths'+i)
      yValues = [float(j) for j in yaxis.text.split()]
      if cCount == 0:
        diffV = yValues
      else:
        tempmaxlist = list(np.array(diffV)-np.array(yValues))
        diffMax.append(max(tempmaxlist))

      yaxis1 = tree.find('xmaxFWHH'+i)
      temp = [float(j) for j in yaxis1.text.split()]
      yValues1.append(temp)
      print 'length of yvalues', yValues
      #xticks = [str(int(x)) for x in xValues]
      #plt.xticks(range(len(xValues)),xticks)
      #plt.yticks(range(len(yValues)),[str(int(y)) for y in yValues])
      plt.xlabel('time')
      plt.ylabel('FWHH')
      plt.plot(xValues,yValues, color = cList[cCount])
      cCount = cCount + 1
    #tree = ET.parse('data.xml')
    #root = tree.getroot()
    #new = ET.SubElement(root,'')
    #new.text = 'sdsd'
    #tree.write('data.xml')
    print diffMax
    plt.figure(2)
    plt.plot(yValues1)
    plt.figure(3)
    plt.plot(diffMax)
    print( "Hit 'enter' to exit" )
    raw_input()


def main():
    xmlfilename = 'data.xml'
    tree = ET.parse(xmlfilename)
    #diffList = ['0.05.01.0','0.05.00.5','0.05.00.2','0.05.00.1','0.05.00.0']
    diffList = ['1e-12bind0.1unbind1e-06dissoc0.1assoc0.0','1e-12bind0.1unbind1e-06dissoc0.0assoc0.0','1e-12bind0.1unbind1e-06dissoc0.01assoc0.0']
    plotXML(tree,diffList)


if __name__ == '__main__':
    main()
