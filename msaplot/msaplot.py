import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms
from matplotlib.transforms import Affine2D
import math
from matplotlib.patheffects import RendererBase
import matplotlib.patheffects as PathEffects
import seaborn as sns

# Calculate the font size that match the unit length on x-axis or y-axis
def CalculateFontsize(xlim, ylim, ax, fig, rows, cols, unit_scale=1):
    # Get axis limits
    axXlim = ax.get_xlim()
    axYlim = ax.get_ylim()

    # Get figure dimensions in pixels
    fig_width, fig_height = fig.get_size_inches() * fig.dpi

    # Get axis dimensions in pixels
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    ax_width, ax_height = bbox.width * fig.dpi, bbox.height * fig.dpi

    # Calculate font size proportional to axis units
    fontsize_x = unit_scale * ax_width / (axXlim[1] - axXlim[0]) / cols * (xlim[1] - xlim[0])
    fontsize_y = 0.8 * unit_scale * ax_height / (axYlim[1] - axYlim[0]) / rows * (ylim[1] - ylim[0])

    # Use the minimum of the two to keep the font size consistent
    fontsize = min(fontsize_x, fontsize_y)

    return fontsize


def GetStartEnd(msa, start, end):
    length = len(msa[0])
    if (start == None):
        start = 0
    elif (start < 0):
        start = length + start 
    if (end == None):
        end = length - 1 
    elif (end < 0):
        end = length + end 
    return start, end


def GetColorMap(msa, color_order = None, palette = None):
    color_map = {}
    if (color_order != None):
        for i,c in enumerate(color_order):
            color_map[c] = sns.color_palette(palette)[i]
            
    # The case some of the alphabet color is not specified in the order     
    for a in msa:
        for c in a:
            if c not in color_map:
                size = len(color_map)
                color_map[c] = sns.color_palette(palette)[size]
    color_map['-'] = [1,1,1] # white
    color_map['.'] = [1,1,1] # white
    return color_map


def DrawMSA(msa, seq_names = None, start = None, end = None,
            axlim = None, color_map = None, palette=None, ax=None, fig=None):
    # Get the canvas attributes.
    ax = ax or plt.gca()
    
    fig = fig or ax.get_figure()
    renderer = fig.canvas.get_renderer()
    
    height = len(msa)
    length = len(msa[0])
    
    # start, end: draw the [start,end] (both inclusive) region of the MSA
    start, end = GetStartEnd(msa, start, end)
    
    if (axlim == None):
        fontsize = CalculateFontsize(ax.get_xlim(), ax.get_ylim(), ax, fig, height, end - start + 1)
    else:
        fontsize = CalculateFontsize(axlim[0], axlim[1], ax, fig, height, end - start + 1)
    
    color_map = GetColorMap(msa, color_order=None, palette=palette)
    
    lengthUnit = 1 / (end - start + 1)
    heightUnit = 1 / height 
    if (axlim != None):
        lengthUnit = (axlim[0][1] - axlim[0][0]) / (end - start + 1)
        heightUnit = (axlim[1][1] - axlim[1][0]) / height
    
    for i, a in enumerate(msa):
        for j,c in enumerate(a[start:end+1]):
            text = ax.text(x=(j + 0.5)*lengthUnit, y=(i+0.5) * heightUnit, s=c, color="black",
                va="center_baseline", ha="center", fontsize=fontsize, 
                transform=ax.transAxes if axlim == None else ax.transData)
            linewidth = min(2,fontsize/50)
            text.set_path_effects([PathEffects.withStroke(linewidth=linewidth, 
                                                         foreground='w')])
            ax.add_patch( patches.Rectangle(xy=(j * lengthUnit, i * heightUnit),
                                           width = lengthUnit, height=heightUnit,
                                          facecolor=color_map[c], linewidth=linewidth, edgecolor="white",
                                          transform=ax.transAxes if axlim == None else ax.transData))
    if (axlim == None):
        ax.set_xlim(-0.5, end - start + 1 - 0.5)
        ax.set_ylim(0-0.5, height-0.5)
        
        # Set the x ticks adaptively at point easy to count
        ticks = []
        tickLabels = []
        
        candidateSteps = [1,5,10,20,50,100,500,1000,5000,10000]
        step = 1
        for i,s in enumerate(candidateSteps):
            if (s * 5 > end - start + 1):
                if (i > 0):
                    step = candidateSteps[i - 1]
                break
        
        tickStart = (int)(start / step) * step + step
        if (tickStart != start):
            ticks.append(0)
            tickLabels.append(start+1)
        for i in range(start, end + 1, step):
            ticks.append(i - start)
            tickLabels.append(i+1)
        if (tickLabels[-1] != end - 1):
            ticks.append(end - start)
            tickLabels.append(end+1)
        ax.set_xticks(ticks, tickLabels)
        
        # Set the y
        ticks = []
        tickLabels = []
        for i in range(height):
            ticks.append(i)
            tickLabels.append(i if seq_names == None else seq_names[i])
        ax.set_yticks(ticks, tickLabels)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)    
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    return ax, color_map

def GetConsensus(msa, start = None, end = None):
    start, end = GetStartEnd(msa, start, end)
    consensus = []
    consensusComposition = []
    for j in range(start, end + 1):
        composition = {}
        if (j >= len(msa[0])):
            consensus.append('')
            consensusComposition.append({"":0})
            continue
            
        for i,a in enumerate(msa):
            if (a[j] not in composition):
                composition[a[j]] = 0
            composition[a[j]] += 1
        result = ""
        maxCnt = 0
        for c in composition:
            if (composition[c] > maxCnt):
                maxCnt = composition[c]
                result = c
        consensus.append(result)
        consensusComposition.append(composition)
    return consensus, consensusComposition

def SetConsensusAxTicks(consensus, consensusComposition, ax):
    ticks = []
    tickLabels = []
    for i in range(len(consensus)):
        ticks.append(i)
        label = consensus[i]
        for c in consensusComposition[i]:
            if (c == label):
                continue
            if (consensusComposition[i][c] == consensusComposition[i][label]):
                label = 'X'
                break
        tickLabels.append(label)
    ax.set_xticks(ticks, tickLabels)

def DrawConsensusHisto(msa, color = [0, 0, 1], color_map = None, start = None, end = None, ax=None):
    ax = ax or plt.gca()
    consensus, consensusComposition = GetConsensus(msa, start, end)
    binHeight = []
    colors = []
    for i in range(len(consensus)):
        if (sum(consensusComposition[i].values()) == 0):
            binHeight.append(0)
        else:
            binHeight.append( consensusComposition[i][consensus[i]] / 
                         sum(consensusComposition[i].values()) )
        colors.append([1-binHeight[i] * (1-color[0]), 
                       1-binHeight[i] * (1-color[1]), 
                       1-binHeight[i] * (1-color[2])]) # color-base in blue
    ax.bar(x=list(range(0, len(consensus))), height=binHeight, 
           color = colors,
           width=0.95)
    ax.set(ylim=(0, 1))
    ax.set(xlim=(-0.5, len(consensus)-0.5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_yaxis().set_ticks([])
    SetConsensusAxTicks(consensus, consensusComposition, ax)


# Code from pyseqlogo for scale font to one direction
class Scale(RendererBase):
    """Scale alphabets using affine transformation"""

    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = Affine2D().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def CalculateEntropy(count):
    s = sum(count)
    if (s == 0):
        return 0
    return sum([-c/s * math.log2(c/s) for c in count])

def DrawSeqLogo(msa, color_map, alphabet_size = None, start = None, end = None, ax=None):
    ax = ax or plt.gca()
    fig = ax.get_figure()
    start, end = GetStartEnd(msa, start, end)
    
    alphabet_size = len(color_map) - 2 # 2 for "-" and "."
    consensus, consensusComposition = GetConsensus(msa, start, end)

    # Definition from https://en.wikipedia.org/wiki/Sequence_logo
    # Search for appropriate height. 
    r = []
    adjuste = 1 / math.log(2) * (alphabet_size - 1) / (2 * len(msa[0]))
    for i,c in enumerate(consensus):
        entropy = CalculateEntropy(list(consensusComposition[i].values()))
        r.append(math.log2(alphabet_size) - (entropy + adjuste))
    
    ax.set_xlim(-0.5, end-start+1-0.5)
    ax.set_ylim(0, max(r))
    
    fontsize = CalculateFontsize(ax.get_xlim(), ax.get_ylim(), 
                                 ax, fig, 1, end - start + 1)
    
    lengthUnit = 1
    for i,c in enumerate(consensusComposition):
        prevy = 0
        totalCount = sum(list(c.values()))
        if (totalCount == 0):
            continue
        for j,item in enumerate(sorted(c.items(), key=lambda x:x[1], reverse=True)):
            k = item[0]
            v = item[1]
            text = ax.text(x=i * lengthUnit, y=prevy, s=k, fontsize=fontsize,
               va="baseline", ha="center", color = (color_map[k] if (k not in ['.','-']) else [0,0,0]),
               transform=ax.transData)
            
            height = v / totalCount * r[i]
            tbox = text.get_window_extent(text._renderer).transformed(ax.transData.inverted())
            scale = height / (tbox.y1 - prevy) 
            #print(i, j, height, scale, tbox)
            prevy = prevy + height 
            text.set_path_effects([Scale(1.0, scale)])
    ax.axis('off')
    SetConsensusAxTicks(consensus, consensusComposition, ax)

# Add the annotation for the sequence alignment
#  annotations: a list of numbers like [['a',0,3]]: msa[0..3] (both inclusive) is annotated as the name 'a'
def DrawAnnotation(msa, annotations, color_map=None,start = None, end = None, ax=None):
    ax = ax or plt.gca()
    fig = ax.get_figure()
    start, end = GetStartEnd(msa, start, end)
     
    ax.set_xlim(start - 0.5, end + 0.5)
    ax.set_ylim(0, 1)
    fontsize = CalculateFontsize(ax.get_xlim(), ax.get_ylim(), ax, fig, 1, end - start + 1)
    
    for a in annotations:
        text = ax.text(x=(a[1]+a[2])/2, y=0.5, s=a[0], fontsize=fontsize,
                       va="center", ha="center", color="black", clip_on=True)
        tbox = text.get_window_extent(text._renderer).transformed(ax.transData.inverted())
        # Draw the bracket
        ax.plot([a[1], a[1]], [0, 1], color="black", clip_on=True)
        ax.plot([a[2], a[2]], [0, 1], color="black", clip_on=True)
        ax.plot([a[1], tbox.x0], [0.5, 0.5], color="black", clip_on=True)
        ax.plot([tbox.x1, a[2]], [0.5, 0.5], color="black", clip_on=True)

    ax.axis('off')

# Draw multipanel MSA
def DrawComplexMSA(msa, panels=[], seq_names = None, panel_height_ratios=None, panel_params=[], 
                   color_map=None, start=None, end=None, wrap=None, figsize=None):
    color_map = color_map or GetColorMap(msa)
    start,end = GetStartEnd(msa, start,end)
    wrap = wrap or (end - start + 1)
    chunks = math.ceil((end - start + 1) / wrap)
    
    height_ratios = None
    if (panel_height_ratios is None):
        panel_height_ratios = []
        for p in panels:
            if (p == DrawMSA):
                panel_height_ratios.append(len(msa))
            elif (p == DrawAnnotation):
                panel_height_ratios.append(0.5)
            else:
                panel_height_ratios.append(1)
    height_ratios = panel_height_ratios * chunks
    fig,axes = plt.subplots(len(panels) * chunks, 1, constrained_layout=True,
                            figsize=figsize, height_ratios = height_ratios)
    
    axidx = 0
    for i in range(start, end + 1, wrap):
        for j,func in enumerate(panels):
            extraParam = panel_params[j].copy()
            if func is DrawMSA:
                extraParam['seq_names'] = seq_names
            func(msa, color_map = color_map, start=i, end=i+wrap-1,ax=axes[axidx], **extraParam)
            axidx += 1
    return axes
