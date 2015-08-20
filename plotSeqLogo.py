import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['pdf.fonttype'] = 42

import matplotlib.pyplot as plt
from matplotlib import transforms
import matplotlib.patheffects

import operator

# from http://www.python-forum.de/viewtopic.php?f=1&t=30856
class Scale(matplotlib.patheffects._Base):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy
       
    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine=affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)
 
def plotLogo(fig, ax, heights, charXDist=10):
    """ draw a sequence logo onto axis. heights is a list of dictionaries with letter -> float 
    >>> freqs = [{"A":0.5, "C":0.3}, {"T":0.3, "G":0.7}]
    >>> fig.set_size_inches(3,0.5*len(freqs))
    >>> plotLogo(fig, ax, freqs)
    """
    #ax.set_xlim(0,len(heights))
    ax.set_ylim(0,1)
    plt.axis('off')
    t = plt.gca().transData
    xPos = 0.0
    charToCol = {"C":"b", "A":"r", "T":"k", "G":"y"}
    for yDict in heights:
        charFreqs = yDict.items()
        charFreqs.sort(key=operator.itemgetter(1))
        for char, freq in charFreqs:
            col = charToCol.get(char.upper(), "k")
            text = plt.text(0, 0, char, transform=t, fontsize=50, color=col, family="monospace")
            text.set_path_effects([Scale(1,freq)])
            fig.canvas.draw()
            ex = text.get_window_extent(fig.canvas.renderer)
            t = transforms.offset_copy(text._transform, y=ex.height*freq, units='dots')

        xPos += ex.width+charXDist
        t = plt.gca().transData
        t = transforms.offset_copy(t, x=xPos, units='dots')
            #text2 = plt.text(0, 0, "C", transform=t)
            #text2.set_path_effects([Scale(1,0.5)])

def main():
    #fig = plt.gcf()
    fig = plt.figure()
    ax = plt.subplot()
    freqs = [{"A":0.5, "C":0.3}, {"T":0.3, "G":0.7}]
    fig.set_size_inches(3,0.5*len(freqs))
    plotLogo(fig, ax, freqs)
    #ax = plt.gcf()
    #print ex
    #text.draw(fig.canvas)
    plt.savefig("temp.png")
    print "wrote to temp.png"
    #print ex.height()



main()
