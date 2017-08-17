import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
 
def plotXML(tree,trials):
    xaxis = tree.find('times')
    plt.ion()
    fig = plt.figure( figsize=(12,10) )
    xValues = [float(j) for j in xaxis.text.split()]
    for i in trials:
      yaxis = tree.find('widths'+i)
      yValues = [float(j) for j in yaxis.text.split()]
      print 'length of yvalues', yValues
      #xticks = [str(int(x)) for x in xValues]
      #plt.xticks(range(len(xValues)),xticks)
      #plt.yticks(range(len(yValues)),[str(int(y)) for y in yValues])
      plt.plot(xValues,yValues)
    #tree = ET.parse('data.xml')
    #root = tree.getroot()
    #new = ET.SubElement(root,'')
    #new.text = 'sdsd'
    #tree.write('data.xml')
    print( "Hit 'enter' to exit" )
    raw_input()


def main():
    xmlfilename = 'data.xml'
    tree = ET.parse(xmlfilename)
    diffList = ['0.01.01.0','1e-121.01.0']
    plotXML(tree,diffList)


if __name__ == '__main__':
    main()
