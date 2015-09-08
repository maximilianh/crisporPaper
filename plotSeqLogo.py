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
    ax.set_ylim(-1,1)
    ax.axis('off')
    #t = plt.gca().transData
    t = ax.transData
    xPos = 0.0
    charToCol = {"C":"b", "A":"r", "T":"k", "G":"y"}
    for yDict in heights:
        if len(yDict)==0:
            yDict["A"]=0.0

        charFreqs = yDict.items()
        charFreqs.sort(key=operator.itemgetter(1))
        #yZeroT = transforms.offset_copy(t, units="dots")
        lastFreq = None
        # need to draw something, otherwise we don't know the width to advance
        for char, freq in charFreqs:
            # jump back to y=0 when we switch from negative to positive freqs
            #print char, freq
            if lastFreq == None or (lastFreq < 0 and freq > 0):
                #t = transforms.offset_copy(text._transform, y=0, units='dots')
                #t = yZeroT
                t = transforms.offset_copy(ax.transData, x=xPos, units='dots')
            lastFreq = freq

            alpha = 1.0
            if freq < 0:
                alpha = 0.5
            col = charToCol.get(char.upper(), "k")

            text = ax.text(0, 0, char, transform=t, fontsize=50, color=col, family="monospace", alpha=alpha)
            text.set_path_effects([Scale(1,freq)])
            fig.canvas.draw()
            ex = text.get_window_extent(fig.canvas.renderer)
            t = transforms.offset_copy(text._transform, y=ex.height*freq, units='dots')

        xPos += ex.width+charXDist
        #t = transforms.offset_copy(ax.transData, x=xPos, units='dots')
        #text2 = plt.text(0, 0, "C", transform=t)
        #text2.set_path_effects([Scale(1,0.5)])

def main():
    #fig = plt.gcf()
    fig = plt.figure()
    ax = plt.subplot()
    freqs = [{"A":-0.5, "C":-0.3, "G":-0.3, "C":0.2}, {"T":-0.3, "G":-0.7, "C":-0.2}]
    fig.set_size_inches(3,0.5*len(freqs))
    plotLogo(fig, ax, freqs)
    #ax = plt.gcf()
    #print ex
    #text.draw(fig.canvas)
    plt.savefig("temp.png")
    print "wrote to temp.png"
    #print ex.height()



if __name__=="__main__":
    main()
