import argparse
import numpy as np
import statistics as stat
from operator import itemgetter
import math
import xlsxwriter as xls
from scipy import stats
import os,time,sys
from PyQt5.QtCore import QObject,pyqtSignal
# Our modules
from exon_coordinates import getExonCoordinates
from patient_table import getSampleInfo
from amplicons_file import getAmpliconCoordinates,matchAmpliconsToExonsAndIntrons
from plot import getXticks,showPlot
from input_file import readInputFile
from output_file import createOutputFile
from cnv_calling import detect_CNVs

class BRACNAC(QObject):
    
    finished=pyqtSignal()
    done_percent=pyqtSignal(int)
    
    def __init__(self,inFile,ampliconFile,patFile,outFile,
                 refVersion,notClust,hardMinScore,
                 hardMaxPvalue,minScore,maxPvalue,minCov,
                 permNum,exonCoveredWhole,exonCoveredPart,exonNonCovered,
                 delTh1,delTh2,duplTh1,duplTh2):
        # Intialization, checking parameters
        super(BRACNAC,self).__init__()
        self.inFile=inFile
        self.ampliconFile=ampliconFile
        self.patFile=patFile
        if outFile=='':
            outFile=inFile[:inFile.rfind('.')]+'.CNVs.xls'
        self.outFile=outFile
        self.refVersion=refVersion
        self.notClust=notClust
        self.hardMinScore=hardMinScore
        self.hardMaxPvalue=hardMaxPvalue
        self.minScore=minScore
        self.maxPvalue=maxPvalue
        self.minCov=minCov
        self.permNum=permNum
        self.exonCoveredWhole=exonCoveredWhole
        self.exonCoveredPart=exonCoveredPart
        self.exonNonCovered=exonNonCovered
        self.delTh1=delTh1
        self.delTh2=delTh2
        self.duplTh1=duplTh1
        self.duplTh2=duplTh2
        
    def run(self):
        if self.refVersion not in ['hg19','hg38']:
            print('ERROR! Reference genome can be only "hg19" or "hg38"')
            exit(1)
        # Getting CDS coordinates
        brca1exons,brca2exons,posToExon1,posToExon2=getExonCoordinates(version=self.refVersion)

        # Getting patient table
        if self.patFile:
            sampleIds=getSampleInfo(self.patFile)
        else:
            sampleIds=None
        # Getting amplicon coordinates
        ampliconInfo=getAmpliconCoordinates(self.ampliconFile,posToExon1,posToExon2)
        amplToExon,exonToAmpls,amplToIntron=ampliconInfo[:3]
        intronToAmpls,ampliconsOrder,multToAmpls,amplToMults=ampliconInfo[3:7]
        amplNumUserSortedToAmplName,amplNumCoordSortedToAmplName=ampliconInfo[7:9]
        amplNameToAmplNumSortedByCoord,amplNameToAmplNum=ampliconInfo[9:11]
        brca1AmpliconNum=ampliconInfo[11]

        # Make x axis ticks
        xticks1,xticks2=getXticks(amplToExon,
                                  amplToIntron,
                                  ampliconsOrder)
        # Read input
        inputData=readInputFile(self.inFile,
                                self.outFile,
                                multToAmpls,
                                amplToMults,
                                amplNumUserSortedToAmplName,
                                amplNumCoordSortedToAmplName,
                                amplNameToAmplNumSortedByCoord,
                                amplNameToAmplNum,
                                self.minCov,
                                self.notClust)
        data3,normCovPatients,medianCovs=inputData[:3]
        lowCovAmplicons,maxi,maxj,colNames,sampleNames=inputData[3:8]
        dirName,outFilePart,minYs,maxYs,allValVars=inputData[8:]

        # Create file for output
        wb,ws,f0,f1=createOutputFile(self.outFile,
                                     colNames,
                                     amplNumUserSortedToAmplName,
                                     amplNameToAmplNumSortedByCoord,
                                     amplToExon,
                                     amplToIntron)
        # start rowNum from which we begin to write results
        rowNum=2
        # We perform calling CNV in two steps:
        # 1) We search for obvious CNVs
        # 2) We search all CNVs
        for step in range(2):
            # Contains list of samples for which we call CNVs
            permSamples=normCovPatients[:]
            # If this is the 2nd step
            if step==1:
                # If this is the 2ns step, we remove samples with obvious CNVs,
                # because they bias normalization
                for sample in obviousCnvSamples:
                    permSamples.remove(sample)
                print('\nSearching all CNVs...')
            # If it is the 1st step
            else:
                # And begin to save samples with obvious CNVs
                obviousCnvSamples=set()
                print('Searching CNVs with the highest P-value...')
            # Go through all samples
            for i in range(maxi)[:]:
                ads=[]
                ais=[]
                # We analyze only samples with enough median coverage
                if i in normCovPatients:
                    # Initially we set color for all dots blue
                    ampliconsColor=['b']*len(data3[i])
                    # Dots for amplicons with low median coverage
                    # are set to black.
                    # amplicon numbers in lowCovAmplicons starts from 0
                    for lca in lowCovAmplicons:
                        ampliconsColor[lca]='k'
                    detectCnvOut=detect_CNVs(data3,i,
                                             permSamples,
                                             sampleIds,sampleNames,
                                             obviousCnvSamples,
                                             ampliconsColor,
                                             lowCovAmplicons,
                                             medianCovs,allValVars,
                                             amplToExon,
                                             exonToAmpls,
                                             amplToIntron,
                                             brca1AmpliconNum,
                                             step,self)

                    if step==1:
                        plotTitle,newCols,ampliconsColor=detectCnvOut
                    else:
                        obviousCnvSamples=detectCnvOut
                elif step==1:
                    if sampleIds!=None:
                        newCols=[sampleNames[i]]+list(map(str,data3[i,]))+[sampleIds[sampleNames[i].replace('patient_','')],str(medianCovs[i]),allValVars[i],'Low_Coverage','-','-','-','-','-']
                    else:
                        newCols=[sampleNames[i]]+list(map(str,data3[i,]))+['',str(medianCovs[i]),allValVars[i],'Low_Coverage','-','-','-','-','-']
                if step==1:
                    for j,newCol in enumerate(newCols):
                        f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center'})
                        if 1<=j<=len(data3[i,]):
                            if 2.4<float(newCol)<=2.6 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#ff00ff'})
                            elif 2.6<float(newCol)<=2.8 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#ff00cc'})
                            elif 2.8<float(newCol)<=3.0 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#ff0099'})
                            elif 3.0<float(newCol)<=4.0 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#ff0066'})
                            elif 4.0<float(newCol) and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#ff0033'})
                            elif 1.6<float(newCol)<=1.8 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#00ffff'})
                            elif 1.4<float(newCol)<=1.6 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#00ddff'})
                            elif 1.2<float(newCol)<=1.4 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#00bbff'})
                            elif 1.0<float(newCol)<=1.2 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#0099ff'})
                            elif 0.5<float(newCol)<=1.0 and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#0077ff'})
                            elif 0.5>float(newCol) and i in normCovPatients:
                                f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,'align':'center','bg_color':'#0055ff'})
                            try:
                                if np.isnan(float(newCol)):
                                    ws.write(rowNum,j,newCol,f2)
                                else:
                                    ws.write(rowNum,j,round(float(newCol),3),f2)
                            except TypeError:
                                print(newCol,type(newCol))
                                exit(0)
                        elif j==0:
                            ws.write(rowNum,j,newCol,f1)
                        else:
                            ws.write(rowNum,j,newCol,f2)
                    rowNum+=1
                    if i in normCovPatients:
                        # Plot all values of this patient
                        if sampleIds!=None:
                            showPlot(data3[i,:],plotTitle,
                                     ampliconsColor,(minYs,maxYs),
                                     ''.join([dirName,outFilePart,'_plots/',
                                              outFilePart,'_',sampleNames[i],'_',
                                              sampleIds[sampleNames[i].replace('patient_','')],
                                              '.png']),
                                     xticks1,xticks2)
                        else:
                            showPlot(data3[i,:],plotTitle,
                                     ampliconsColor,(minYs,maxYs),
                                     ''.join([dirName,outFilePart,'_plots/',
                                              outFilePart,'_',sampleNames[i],
                                              '.png']),
                                     xticks1,xticks2)

                self.showPercWork(i+1,maxi)
        wb.close()
        self.finished.emit()
        
    # Shows percent of work done
    def showPercWork(self,done,allWork):
        percDoneWork=math.floor((done/allWork)*100)
        self.done_percent.emit(percDoneWork)
        sys.stdout.write("\r"+str(percDoneWork)+"%")
        sys.stdout.flush()


if __name__=='__main__':
    # Section of input arguments
    par=argparse.ArgumentParser(description='This script normalize coverage data by multiplexes')
    par.add_argument('--input-file','-in',
                     dest='inFile',type=str,
                     help='file with amplicon coverages with following columns: Patient#, Patient_ID, Barcodes, Median_Coverage, Number_<30, amplicon#1...'
                     '(Patient_ID, Barcodes, Median_Coverage, Number_<30 are optional columns)',required=True)
    par.add_argument('--amplicon-file','-af',
                     dest='ampliconFile',type=str,
                     help='TSV-file with coordinates for each amplicon and numbers of multiplexes',required=True)
    par.add_argument('--patient-table','-pt',
                     dest='patFile',type=str,
                     help='TSV-file with patient IDs and names. Column names are required',required=False)
    par.add_argument('--output-file','-out',
                     dest='outFile',type=str,
                     help='output file',required=True)
    par.add_argument('--reference-genome-version','-ref',
                     dest='refVersion',type=str,
                     help='version of the reference genome (hg19 or hg38). Default: hg19',
                     default='hg19',required=False)
    par.add_argument('--not-clusterize','-notclust',
                     dest='notClust',action='store_true',
                     help='use this parameter if you do not want to clusterize samples during normalization',
                     required=False)
    par.add_argument('--hard-score-treshold','-hard_score',
                     dest='hardMinScore',type=float,
                     help='minimal hard score for large rearrangement detection (default: 9.9)',
                     default=9.9,required=False)
    par.add_argument('--hard-pvalue-treshold','-hard_pvalue',
                     dest='hardMaxPvalue',type=float,
                     help='maximal p-value for hard filtering mutations (default: 0.01)',
                     default=0.01,required=False)
    par.add_argument('--score-treshold','-score',
                     dest='minScore',type=float,
                     help='minimal score for large rearrangement detection (default: 2)',
                     default=2,required=False)
    par.add_argument('--pvalue-treshold','-pvalue',
                     dest='maxPvalue',type=float,
                     help='maximal p-value for filtering mutations (default: 0.05)',
                     default=0.05,required=False)
    par.add_argument('--coverage-treshold','-cov',
                     dest='minCov',type=int,
                     help='minimal median coverage for patients (default: 100)',
                     default=100,required=False)
    par.add_argument('--permutation-number','-perm',
                     dest='permNum',type=int,
                     help='number of permutations for calculating p-value (default: 1000)',
                     default=1000,required=False)
    par.add_argument('--exon-covered-whole-score','-whole',
                     dest='exonCoveredWhole',type=int,
                     help='value for calculating score when whole exon is covered (default: 1)',
                     default=1,required=False)
    par.add_argument('--exon-covered-part-score','-part',
                     dest='exonCoveredPart',type=float,
                     help='value for calculating score when part of exon is covered (default: 0.5)',
                     default=0.5,required=False)
    par.add_argument('--exon-non-covered-score','-non',
                     dest='exonNonCovered',type=int,
                     help='value for calculating score when exon is not covered (default: 0)',
                     default=0,required=False)
    par.add_argument('--deletion-threshold1','-del1',
                     dest='delTh1',type=float,
                     help='normalized value of coverage for considering '
                          'an amplicon as likely deleted (default: 1.3)',
                     default=1.3,required=False)
    par.add_argument('--deletion-threshold2','-del2',
                     dest='delTh2',type=float,
                     help='normalized value of coverage for considering '
                          'an amplicon as probably deleted (default: 1.4)',
                     default=1.8,required=False)
    par.add_argument('--duplication-threshold1','-dupl1',
                     dest='duplTh1',type=float,
                     help='normalized value of coverage for considering '
                          'an amplicon as likely duplicated (default: 2.8)',
                     default=2.7,required=False)
    par.add_argument('--duplication-threshold2','-dupl2',
                     dest='duplTh2',type=float,
                     help='normalized value of coverage for considering '
                          'an amplicon as probably duplicated (default: 2.7)',
                     default=2.4,required=False)
    args=par.parse_args()
    bracnac=BRACNAC(args.inFile,args.ampliconFile,args.patFile,args.outFile,
                    args.refVersion,args.notClust,args.hardMinScore,args.hardMaxPvalue,
                    args.minScore,args.maxPvalue,args.minCov,args.permNum,
                    args.exonCoveredWhole,args.exonCoveredPart,args.exonNonCovered,
                    args.delTh1,args.delTh2,args.duplTh1,args.duplTh2) 
    bracnac.run()
