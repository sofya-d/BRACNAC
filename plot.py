# This script contains functions for drawing plot

import math
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

def getXticks(amplToExon,amplToIntron,ampliconsOrder):
    # Contains numbers of amplicon
    xticks1=[]
    # Contains numbers of exons/introns
    xticks2=[]
    # Exon or intron number for previous amplicon number
    prevExonIntronNum=None
    # Previous amplicon number
    prevAmplNum=None
    i=0
    for amplNum in ampliconsOrder:
        i+=1
        # Number of exon or intron for current amplicon number
        exonIntronNum=None
        if amplNum in amplToExon.keys():
            exonIntronNum='ex'+str(amplToExon[amplNum])
        elif amplNum in amplToIntron.keys():
            exonIntronNum='in'+str(amplToIntron[amplNum])
        else:
            print('ERROR! The following amplicon was not related to any exon or intron:')
            print(amplNum)
            exit(1)
        # If it is the 1st amplicon number
        if prevExonIntronNum==None:
            prevExonIntronNum=exonIntronNum
            prevAmplNum=amplNum
            continue
        # If the next intron or exon started
        if exonIntronNum!=prevExonIntronNum or i==len(ampliconsOrder):
            if i==len(ampliconsOrder):
                prevExonIntronNum=exonIntronNum
                prevAmplNum=amplNum
            xticks1.append(prevAmplNum)
            prevExonIntronNumWithoutExIn=int(prevExonIntronNum.replace('ex','').replace('in',''))
            if 'ex' in prevExonIntronNum:
                prefix='ex'
            elif 'in' in prevExonIntronNum:
                prefix='in'
            else:
                print('ERROR! Unknown type of exon or intron number:')
                print(prevExonIntronNum)
                exit(1)
            if 200>prevExonIntronNumWithoutExIn>100:
                xticks2.append(prefix+str(prevExonIntronNumWithoutExIn-100))
            elif prevExonIntronNumWithoutExIn>200:
                xticks2.append(prefix+str(prevExonIntronNumWithoutExIn-200))
            elif prevExonIntronNumWithoutExIn==0:
                xticks2.append('intron')
            elif prevExonIntronNumWithoutExIn==100:
                xticks2.append('up')
            elif prevExonIntronNumWithoutExIn==200:
                xticks2.append('up')
            else:
                print('ERROR! Unknown xticks2')
                print(prevExonIntronNumWithoutExIn)
                exit(1)
        prevExonIntronNum=exonIntronNum
        prevAmplNum=amplNum
    return(xticks1,xticks2)

# Draws plot with normalized values for each amplicon
def showPlot(vals,text,colors,yperc,fileName,xticks1,xticks2):
    fig,ax1 = plt.subplots( nrows=1, ncols=1 )
    fig.set_size_inches(30,10)
    titleText=plt.title(text,fontsize=20,fontweight='bold')
    titlePos=titleText.get_position()
    titleText.set_position((titlePos[0],titlePos[1]+0.1))
    plt.ylabel('Normalized coverage',size=26,fontweight='bold')
    plt.axis([-2,len(vals),0,4.1])
    x=list(range(len(vals)))
    plt.fill_between(x,yperc[0],yperc[1],facecolor='#C2C2C2',edgecolor='#C2C2C2') #,step='pre'
    plt.scatter(x,vals,c=colors,lw=0,s=120)
    plt.plot(x,vals,'-b')
    plt.plot(x,[2]*len(vals),'k',x,[1]*len(vals),'k-',x,[3]*len(vals),'k-',x,[1.5]*len(vals),'k--',x,[2.5]*len(vals),'k--')
    newXticks1=[]
    newXticks1.append(-0.5)
    for i,xtick in enumerate(xticks1):
        newXticks1.append(xtick+0.5)
    plt.vlines(newXticks1,0,4.1,linestyle='dashed')
    newXticks1=[]
    newXticks2=[]
    for i,xtick in enumerate(sorted(xticks1)):
        if i==0:
            newXticks1.append(0)
        else:
            newXticks1.append((xtick+sorted(xticks1)[i-1]+2)/2-0.5)
        newXticks2.append(xticks2[xticks1.index(xtick)])
        if (i>0 and xticks2[xticks1.index(xtick)]!='up' and
            xticks2[xticks1.index(sorted(xticks1)[i-1])]!='up'):
            curExInNum=int(xticks2[xticks1.index(xtick)].replace('ex','').replace('in',''))
            prevExInNum=int(xticks2[xticks1.index(sorted(xticks1)[i-1])].replace('ex','').replace('in',''))
            # If amplicons of BRCA2 gene started
            if prevExInNum>curExInNum:
                lines=plt.vlines([sorted(xticks1)[i-1]+0.5],0,4.1,
                                 linestyle='solid',
                                 colors=['k'])
                lines.set_linewidth(3)
                plt.text(sorted(xticks1)[i-2]/2,
                         -0.7,'BRCA1 exons and amplicons',ha='center',
                         weight='bold',size=24)
                plt.text((sorted(xticks1)[-1]+sorted(xticks1)[i-2])/2,
                         -0.7,'BRCA2 exons and amplicons',ha='center',
                         weight='bold',size=24)
    ax2=ax1.twiny()
    ax1.set_xticks(newXticks1[::2])
    ax1.set_xticklabels(newXticks2[::2],
                        rotation=60)
    ax1.set_xlim(newXticks1[0]-0.5,
                 xticks1[-1]+0.5)
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.set_xticks(newXticks1[1::2])
    ax2.set_xticklabels(newXticks2[1::2],
                        rotation=60)
    ax1.tick_params(axis='x',labelsize=18,length=5,
                    width=2,direction='out')
    ax2.tick_params(axis='x',labelsize=18,length=5,
                    width=2,direction='out')
    ax2.set_xlim(ax1.get_xlim())
    plt.tick_params(axis='y',labelsize=18,right='off')
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(fileName,
                dpi=300)
    plt.close()
