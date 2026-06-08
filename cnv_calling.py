# This script detects CNVs in normalized values of
# coverage for one sample

import math
import numpy as np
import statistics as stat
import itertools
import pandas as pd

def write_cnv_table(sample: str, cnv_lst: list, score_lst: list, pvalue_lst: list, CNVTable_dirpath: str):
    """
    Create a CNV table with CNVs from the output plot title.
    """
    # Round p-values
    pvalue_lst = [round(float(pvalue), 3) for pvalue in pvalue_lst]
    # Divide each CNV to CNV type, gene, and exon lists to create a data frame later
    cnv_lst_lst = [cnv.split(',') for cnv in cnv_lst]
    cnv_type_lst = []
    gene_lst = []
    exon_lst = []
    for cnv_lst in cnv_lst_lst:
        cnv_types = []
        genes = []
        exons = []
        for cnv in cnv_lst:
            cnv_type, gene, exon = cnv.split('_')
            cnv_types.append(cnv_type)
            genes.append(gene)
            exons.append(exon)
        same_cnv_type = len(set(cnv_types)) == 1
        same_gene = len(set(genes)) == 1
        if same_cnv_type and same_gene:
            cnv_type_lst.append(cnv_types[0])
            gene_lst.append(genes[0])
            exon_lst.append('_'.join(exons))
        else:
            raise ValueError(f'Multiple non-consecutive exon CNV {cnv_lst} has either different CNV type, or gene, or both for different exons in sample {sample}.')
    # Create a data frame
    cnv_df = pd.DataFrame({'gene': gene_lst, 'cnv_type': cnv_type_lst, 'exon': exon_lst, 'score': score_lst, 'pvalue': pvalue_lst})
    cnv_df = cnv_df.sort_values('gene')
    CNVTable_filepath = f'{CNVTable_dirpath}/{sample}_CNV_BRACNAC.tsv'
    cnv_df.to_csv(CNVTable_filepath, sep='\t', index=False)

def detect_CNVs(data3,
                sampleNum,
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
                step,args,
                CNVTable_dirpath):
    if step==1:
        minScore=args.minScore
        maxPvalue=args.maxPvalue
    else:
        minScore=args.hardMinScore
        maxPvalue=args.hardMaxPvalue
    # Extract potential InDels
    # Contains all normalized coverage values of this sample
    vals=np.array(data3[sampleNum])
    # Contains amplicon numbers of BRCA1
    # for which normalized coverage is less than 2
    dels1=[i for i,v in enumerate(vals[:brca1AmpliconNum]<2) if v]
    # Contains amplicon numbers of BRCA2
    # for which normalized coverage is less than 2
    dels2=[i+brca1AmpliconNum for i,v in enumerate(vals[brca1AmpliconNum:]<2) if v]
    # Contains amplicon numbers of BRCA1
    # for which normalized coverage is more than 2
    inss1=[i for i,v in enumerate(vals[:brca1AmpliconNum]>2) if v]
    # Contains amplicon numbers of BRCA2
    # for which normalized coverage is more than 2
    inss2=[i+brca1AmpliconNum for i,v in enumerate(vals[brca1AmpliconNum:]>2) if v]
    # Convert lists of amplicon numbers to ranges of numbers
    delRanges1=[]
    delRanges2=[]
    insRanges1=[]
    insRanges2=[]
    for deletion in dels1:
        delRanges1.append((deletion,deletion))
    for deletion in dels2:
        delRanges2.append((deletion,deletion))
    for insertion in inss1:
        insRanges1.append((insertion,insertion))
    for insertion in inss2:
        insRanges2.append((insertion,insertion))
    ##################################
    ##### DELS for BRCA1 #############
    ##################################
    newRanges=delRanges1[:]
    # Join amplicon ranges
    # We go through all amplicons with coverage less than 2
    for k,delRange1 in enumerate(delRanges1):
        # and compare with all other amplicons with coverage less than 2
        for delRange2 in delRanges1[k+1:]:
            # If start of the 2nd amplicon range is less than or
            # equal to the end of the 1st amplicon range +2
            if (delRange2[0]<=delRange1[1]+2 and
                # and 
                (int((delRange2[0]+delRange1[1])/2) in lowCovAmplicons or
                 data3[sampleNum,int((delRange2[0]+delRange1[1])/2)]<=2)):
                newRanges.append((delRange1[0],delRange2[1]))
                delRange1=(delRange1[0],delRange2[1])
            else:
                break
    delRanges1=[]
    for delRange in newRanges:
        # Prepare list of normally covered amplicons
        normCovAmplicons=[]
        for i in range(delRange[0],delRange[1]+1):
            if i not in lowCovAmplicons:
                normCovAmplicons.append(i)
        rangeVals=data3[sampleNum,normCovAmplicons]
        if (checkCoveredExons(delRange,
                              amplToExon,
                              exonToAmpls,
                              amplToIntron,
                              args)[0]==args.exonCoveredWhole and
            len(rangeVals[rangeVals<=args.delTh1])>=1 and
            len(rangeVals[rangeVals>args.delTh2])<=math.ceil((delRange[1]-delRange[0]+1)/5) and
            len(rangeVals[rangeVals<=args.delTh2])>=(delRange[1]-delRange[0]+1)/2):
            delRanges1.append(delRange)
##    if sampleNum==71:
##        print(delRanges1,sampleIds[sampleNames[sampleNum].replace('patient_','')])
    ##################################
    ##### DELS for BRCA2 #############
    ##################################
    newRanges=sorted(delRanges2[:])
    # Join amplicon ranges
    for k,delRange1 in enumerate(delRanges2):
        for delRange2 in delRanges2[k+1:]:
            if (delRange2[0]<=delRange1[1]+2 and
                (int((delRange2[0]+delRange1[1])/2) in lowCovAmplicons or
                 data3[sampleNum,int((delRange2[0]+delRange1[1])/2)]<=2)):
                newRanges.append((delRange1[0],delRange2[1]))
                delRange1=(delRange1[0],delRange2[1])
            else:
                break
    delRanges2=[]
    for delRange in newRanges:
        # Prepare list of normally covered amplicons
        normCovAmplicons=[]
        for i in range(delRange[0],delRange[1]+1):
            if i not in lowCovAmplicons:
                normCovAmplicons.append(i)
        rangeVals=data3[sampleNum,normCovAmplicons]
        if (checkCoveredExons(delRange,
                              amplToExon,
                              exonToAmpls,
                              amplToIntron,
                              args)[0]==args.exonCoveredWhole and
            len(rangeVals[rangeVals<=args.delTh1])>=1 and
            len(rangeVals[rangeVals>args.delTh2])<=math.ceil((delRange[1]-delRange[0]+1)/5) and
            len(rangeVals[rangeVals<=args.delTh2])>=(delRange[1]-delRange[0]+1)/2):
            delRanges2.append(delRange)
    ##################################
    ##### DUPLICATIONS for BRCA1 #####
    ##################################
    newRanges=insRanges1[:]
    # Join amplicon ranges
    for k,insRange1 in enumerate(insRanges1):
        for insRange2 in insRanges1[k+1:]:
            if (insRange2[0]<=insRange1[1]+2 and
                (int((insRange2[0]+insRange1[1])/2) in lowCovAmplicons or
                 data3[sampleNum,int((insRange2[0]+insRange1[1])/2)]>=2)):
                newRanges.append((insRange1[0],insRange2[1]))
                insRange1=(insRange1[0],insRange2[1])
            else:
                break
    insRanges1=[]
##    if sampleNum==53:
##        print(newRanges,sampleIds[sampleNames[sampleNum].replace('patient_','')])
    for insRange in newRanges:
        # Prepare list of normally covered amplicons
        normCovAmplicons=[]
        for i in range(insRange[0],insRange[1]+1):
            if i not in lowCovAmplicons:
                normCovAmplicons.append(i)
        rangeVals=data3[sampleNum,normCovAmplicons]
        if (checkCoveredExons(insRange,
                              amplToExon,
                              exonToAmpls,
                              amplToIntron,
                              args)[0]==args.exonCoveredWhole and
            len(rangeVals[rangeVals>=args.duplTh1])>=1 and
            len(rangeVals[rangeVals<args.duplTh2])<=math.ceil((insRange[1]-insRange[0]+1)/5) and
            len(rangeVals[rangeVals>=args.duplTh2])>=(insRange[1]-insRange[0]+1)/2):
            insRanges1.append(insRange)
##    if sampleNum==53:
##        print(insRanges1,sampleIds[sampleNames[sampleNum].replace('patient_','')])
    ##################################
    ##### DUPLICATIONS for BRCA2 #####
    ##################################
    newRanges=insRanges2[:]
    # Join amplicon ranges
    for k,insRange1 in enumerate(insRanges2):
        for insRange2 in insRanges2[k+1:]:
            if (insRange2[0]<=insRange1[1]+2 and
                (int((insRange2[0]+insRange1[1])/2) in lowCovAmplicons or
                 data3[sampleNum,int((insRange2[0]+insRange1[1])/2)]>=2)):
                newRanges.append((insRange1[0],insRange2[1]))
                insRange1=(insRange1[0],insRange2[1])
            else:
                break
    insRanges2=[]
    for insRange in newRanges:
        # Prepare list of normally covered amplicons
        normCovAmplicons=[]
        for i in range(insRange[0],insRange[1]+1):
            if i not in lowCovAmplicons:
                normCovAmplicons.append(i)
        rangeVals=data3[sampleNum,normCovAmplicons]
        if (checkCoveredExons(insRange,
                              amplToExon,
                              exonToAmpls,
                              amplToIntron,
                              args)[0]==args.exonCoveredWhole and
            len(rangeVals[rangeVals>=args.duplTh1])>=1 and
            len(rangeVals[rangeVals<args.duplTh2])<=math.ceil((insRange[1]-insRange[0]+1)/5) and
            len(rangeVals[rangeVals>=args.duplTh2])>=(insRange[1]-insRange[0]+1)/2):
            insRanges2.append(insRange)
    newDelRanges=delRanges1+delRanges2
    newInsRanges=insRanges1+insRanges2

    ##############################################################
    #### Get affetced exons and calculate p-values and scores ####
    ##############################################################
    ads=[]
    adsValues=[]
    adsScores=[]
    bestAds=[]
    ais=[]
    aisValues=[]
    aisScores=[]
    bestAis=[]
    if len(newDelRanges+newInsRanges)>0:
        for delRange in newDelRanges:
            delValue,deletedExons=checkCoveredExons(delRange,
                                                    amplToExon,
                                                    exonToAmpls,
                                                    amplToIntron,
                                                    args)
            if len(deletedExons)==0:
                continue
            adv=data3[sampleNum,delRange[0]:delRange[1]+1]
            # Calculate score
            score=delValue*(sum(4-adv[adv<=args.delTh2]))
            if adv[0]>args.delTh2:
                score-=4
            if adv[-1]>args.delTh2:
                score-=4
            amplNums=list(range(delRange[0],delRange[1]+1))
            amplNumsData=data3[sampleNum][amplNums]
            moreThanThreshNum=len(amplNumsData[amplNumsData>args.delTh1])
            score/=args.delTh1**(moreThanThreshNum/len(amplNums))
            if 102 not in deletedExons:
                if len(amplNums)==1 and len(deletedExons)>0:
                    score/=args.delTh1**2
                score/=args.delTh1**(1/len(deletedExons))
            dists=[]
            for k,amplNum in enumerate(amplNums):
                if k==1:
                    continue
                dists.append(abs(data3[sampleNum][amplNum]-data3[sampleNum][amplNums[k-1]]))
            score/=args.delTh1**stat.median(dists)
            if score<minScore:
                continue
            deletedAmpls=list(np.array(list(range(delRange[0],delRange[1]+1)))+1)
            permNum1=args.permNum
            pValue=calc_pValue(np.array(data3[sampleNum]),
                               np.array(data3[:,amplNums]),
                               amplNums,sampleNum,permNum1,
                               args.permNum-permNum1,
                               permSamples[:],0,
                               score,args)
            if pValue>maxPvalue:
                continue
            elif step==0:
                obviousCnvSamples.add(sampleNum)
            delExs=list(to_ranges(deletedExons))
            deletedExons=[]
            for delEx in delExs:
                # Convert IDs of exons to their numbers: BRCA1 has IDs from 101; BRCA2 from 201
                if delEx[0]>=100 and delEx[0]<200:
                    if delEx[0]==100:
                        deletedExons.append('del_BRCA1_up-ex'+str(delEx[1]-100))
                    else:
                        deletedExons.append('del_BRCA1_ex'+str(delEx[0]-100)+'-ex'+str(delEx[1]-100))
                elif delEx[0]>200:
                    deletedExons.append('del_BRCA2_ex'+str(delEx[0]-200)+'-ex'+str(delEx[1]-200))
            if ','.join(deletedExons) not in ads:
                ads.append(','.join(deletedExons))
                adsScores.append(str(round(score,1)))
                adsValues.append(str(pValue))
                bestAds.append(delEx)
            elif (pValue<float(adsValues[ads.index(','.join(deletedExons))]) or
                  (pValue==float(adsValues[ads.index(','.join(deletedExons))]) and
                   round(score,1)>float(adsScores[ads.index(','.join(deletedExons))]))):
                adsScores[ads.index(','.join(deletedExons))]=str(round(score,1))
                adsValues[ads.index(','.join(deletedExons))]=str(pValue)
            # Change color of deleted amplicons
            for ampl in amplNums:
                ampliconsColor[ampl]='r'
        for insRange in newInsRanges:
            insValue,insertedExons=checkCoveredExons(insRange,
                                                     amplToExon,
                                                     exonToAmpls,
                                                     amplToIntron,
                                                     args)
            if insertedExons==[]:
                continue
            aiv=data3[sampleNum,insRange[0]:insRange[1]+1]
            # Calculate score
            score=insValue*(sum(aiv[aiv>=args.duplTh2]))
            if aiv[0]<args.duplTh2:
                score-=4
            if aiv[-1]<args.duplTh2:
                score-=4
            amplNums=list(range(insRange[0],insRange[1]+1))
            amplNumsData=data3[sampleNum][amplNums]
            lessThanThreshNum=len(amplNumsData[amplNumsData<args.duplTh1])
            score/=args.delTh1**(lessThanThreshNum/len(amplNums))
            if (len(insertedExons)>0 and 101 not in insertedExons):
                if len(amplNums)==1 and len(insertedExons)>0:
                    score/=args.delTh1**2
                score/=args.delTh1**(1/len(insertedExons))
            dists=[]
            for k,amplNum in enumerate(amplNums):
                if k==1:
                    continue
                dists.append(abs(data3[sampleNum][amplNum]-data3[sampleNum][amplNums[k-1]]))
            score/=args.delTh1**stat.median(dists)
            if score<minScore:
                continue
            insertedAmpls=list(np.array(list(range(insRange[0],insRange[1]+1)))+1)
            permNum1=args.permNum
            pValue=calc_pValue(np.array(data3[sampleNum]),
                               np.array(data3[:,amplNums]),
                               amplNums,sampleNum,permNum1,
                               args.permNum-permNum1,
                               permSamples[:],1,score,args)
            if pValue>maxPvalue:
                continue
            elif step==0:
                obviousCnvSamples.add(sampleNum)
            insExs=list(to_ranges(insertedExons))
            insertedExons=[]
            for insEx in insExs:
                # Convert IDs of exons to their numbers: BRCA1 has IDs from 101; BRCA2 from 201
                if insEx[0]>=100 and insEx[0]<200:
                    if insEx[0]==100:
                        insertedExons.append('dupl_BRCA1_up-ex'+str(insEx[1]-100))
                    else:
                        insertedExons.append('dupl_BRCA1_ex'+str(insEx[0]-100)+'-ex'+str(insEx[1]-100))
                elif insEx[0]>200:
                    insertedExons.append('dupl_BRCA2_ex'+str(insEx[0]-200)+'-ex'+str(insEx[1]-200))
            if ','.join(insertedExons) not in ais:
                ais.append(','.join(insertedExons))
                aisScores.append(str(round(score,1)))
                aisValues.append(str(pValue))
                bestAis.append(insEx)
            elif (pValue<float(aisValues[ais.index(','.join(insertedExons))]) or
                  (pValue==float(aisValues[ais.index(','.join(insertedExons))]) and
                   round(score,1)>float(aisScores[ais.index(','.join(insertedExons))]))):
                aisScores[ais.index(','.join(insertedExons))]=str(round(score,1))
                aisValues[ais.index(','.join(insertedExons))]=str(pValue)
            # Change color of amplified amplicons
            for ampl in amplNums:
                ampliconsColor[ampl]='r'
        if step==1:
            if len(ads)>0:
                newAds=[]
                newAdsValues=[]
                newAdsScores=[]
                consideredExons=set()
                for delCNV in sorted(ads,
                                     key=lambda key:(float(adsValues[ads.index(key)]),
                                                     -float(adsScores[ads.index(key)]))):
                    if len(consideredExons.intersection(set(range(bestAds[ads.index(delCNV)][0],
                                                                  bestAds[ads.index(delCNV)][1]+1))))>0:
                        continue
                    newAds.append(delCNV)
                    newAdsValues.append(adsValues[ads.index(delCNV)])
                    newAdsScores.append(adsScores[ads.index(delCNV)])
                    consideredExons.update(set(range(bestAds[ads.index(delCNV)][0],
                                                     bestAds[ads.index(delCNV)][1]+1)))
                if sampleIds!=None:
                    newCols=[sampleIds[sampleNames[sampleNum].replace('patient_','')],
                             str(medianCovs[sampleNum]),allValVars[sampleNum],
                             ';'.join(newAds),';'.join(newAdsScores),';'.join(newAdsValues)]
                else:
                    newCols=['',str(medianCovs[sampleNum]),allValVars[sampleNum],
                             ';'.join(newAds),';'.join(newAdsScores),';'.join(newAdsValues)]
            else:
                if sampleIds!=None:
                    newCols=[sampleIds[sampleNames[sampleNum].replace('patient_','')],str(medianCovs[sampleNum]),allValVars[sampleNum],'-','-','-']
                else:
                    newCols=['',str(medianCovs[sampleNum]),allValVars[sampleNum],'-','-','-']
            if len(ais)>0:
                newAis=[]
                newAisValues=[]
                newAisScores=[]
                consideredExons=set()
                for insCNV in sorted(ais,
                                     key=lambda key:(float(aisValues[ais.index(key)]),
                                                     -float(aisScores[ais.index(key)]))):
                    if len(consideredExons.intersection(set(range(bestAis[ais.index(insCNV)][0],
                                                                  bestAis[ais.index(insCNV)][1]+1))))>0:
                        continue
                    newAis.append(insCNV)
                    newAisValues.append(aisValues[ais.index(insCNV)])
                    newAisScores.append(aisScores[ais.index(insCNV)])
                    consideredExons.update(set(range(bestAis[ais.index(insCNV)][0],
                                                     bestAis[ais.index(insCNV)][1]+1)))
                newCols.extend([';'.join(newAis),';'.join(newAisScores),';'.join(newAisValues)])
            else:
                newCols.extend(['-','-','-'])
            newCols=[sampleNames[sampleNum]]+list(map(str,data3[sampleNum,]))+newCols
    elif step==1:
        if sampleIds!=None:
            newCols=[sampleNames[sampleNum]]+list(map(str,data3[sampleNum,]))+[sampleIds[sampleNames[sampleNum].replace('patient_','')],str(medianCovs[sampleNum]),allValVars[sampleNum],'-','-','-','-','-','-']
        else:
            newCols=[sampleNames[sampleNum]]+list(map(str,data3[sampleNum,]))+['',str(medianCovs[sampleNum]),allValVars[sampleNum],'-','-','-','-','-','-']
    if step==1:
        if sampleIds!=None:
            text=sampleNames[sampleNum]+' ('+sampleIds[sampleNames[sampleNum].replace('patient_','')]+')\n'
        else:
            text=sampleNames[sampleNum]+'\n'
        
        # Create variables to write CNV table data to
        sample = sampleNames[sampleNum]
        cnv_lst = []
        score_lst = []
        pvalue_lst = []
        # Create variables to write CNV table data to
        
        if len(ads)>0:
            consideredExons=set()
            for delCNV in sorted(ads,
                                 key=lambda key:(float(adsValues[ads.index(key)]),
                                                 -float(adsScores[ads.index(key)]))):
                if len(consideredExons.intersection(set(range(bestAds[ads.index(delCNV)][0],
                                                              bestAds[ads.index(delCNV)][1]+1))))>0:
                    continue
                text+='Deletion: '+delCNV+' (Score='+adsScores[ads.index(delCNV)]+', p-value='+adsValues[ads.index(delCNV)]+')\n'
                consideredExons.update(set(range(bestAds[ads.index(delCNV)][0],
                                                 bestAds[ads.index(delCNV)][1]+1)))

        # Uppend CNV data to lists
        score = adsScores[ads.index(delCNV)]
        pvalue = adsValues[ads.index(delCNV)]
        cnv_lst.append(delCNV)
        score_lst.append(score)
        pvalue_lst.append(pvalue)
        # Uppend CNV data to lists
        
        else:
            text+='Deletions: No\n'
        if len(ais)>0:
            consideredExons=set()
            for insCNV in sorted(ais,
                                 key=lambda key:(float(aisValues[ais.index(key)]),
                                                 -float(aisScores[ais.index(key)]))):
                if len(consideredExons.intersection(set(range(bestAis[ais.index(insCNV)][0],
                                                              bestAis[ais.index(insCNV)][1]+1))))>0:
                    continue
                text+='Insertion: '+insCNV+' (Score='+aisScores[ais.index(insCNV)]+', p-value='+aisValues[ais.index(insCNV)]+')\n'
                consideredExons.update(set(range(bestAis[ais.index(insCNV)][0],
                                                 bestAis[ais.index(insCNV)][1]+1)))

        # Uppend CNV data to lists
        score = aisScores[ais.index(insCNV)]
        pvalue = aisValues[ais.index(insCNV)]
        cnv_lst.append(insCNV)
        score_lst.append(score)
        pvalue_lst.append(pvalue)
        # Uppend CNV data to lists
        
        else:
            text+='Insertions: No'
        
        write_cnv_table(sample, cnv_lst, score_lst, pvalue_lst, CNVTable_dirpath)
        
        return(text,newCols,ampliconsColor)
    else:
        return(obviousCnvSamples)

def checkCoveredExons(ampls,
                      amplToExon,
                      exonToAmpls,
                      amplToIntron,
                      args):
    # Check for overlapping with exons
    coveredExons={}
    coveredIntrons={}
    for ampl in range(ampls[0],ampls[1]+1):
        if int(ampl) in amplToExon.keys():
            if amplToExon[int(ampl)] in coveredExons:
                coveredExons[amplToExon[int(ampl)]]+=1
            else:
                coveredExons[amplToExon[int(ampl)]]=1
        elif int(ampl) in amplToIntron.keys():
            if amplToIntron[int(ampl)] in coveredIntrons:
                coveredIntrons[amplToIntron[int(ampl)]]+=1
            else:
                coveredIntrons[amplToIntron[int(ampl)]]=1
    # Go through all somehow covered exons (with all covered amplicons or only part)
    affectedExons=[]
    # Stores for each somehow covered exon amplicon number affected
    valuePre=[]
    if 100 in coveredIntrons.keys():
        affectedExons.append(100)
    for key,item in sorted(coveredExons.items(),
                           key=lambda item:item[0]):
        if item==len(exonToAmpls[key]):
            affectedExons.append(key)
            valuePre.append(1)
        else:
            valuePre.append(0)
    if len(valuePre)==valuePre.count(1):
        delValue=args.exonCoveredWhole
    elif len(valuePre)==valuePre.count(0):
        delValue=args.exonNonCovered
    else:
        delValue=args.exonCoveredPart
    return(delValue,affectedExons)

def calc_pValue(patRow,
                amplCols,
                amplNums,
                sampleNum,
                permNum1,
                permNum2,
                normCovPatients,
                inDel,
                minScore,
                args):
    # Remove from patient values that correspond to the potential CNV
    patRow=np.delete(patRow,amplNums,None)
    # Save amplCols shape before removing values (to calculate p-value later)
    initAmplColSize=amplCols.shape[0]
    # Remove from the region set values that correspond to the potential CNV
    amplCols=np.delete(amplCols,sampleNum,0)
    # Correct indices because our current array doesn't include the current sample row
    # so we need to decrease indices for samples that are after the current sample
    for i,normCovPatient in enumerate(normCovPatients):
        if normCovPatient>=sampleNum:
            normCovPatients[i]-=1
    # Prepare array only for samples that have enough coverage
    amplCols=amplCols[normCovPatients,:]
    if amplCols.shape[0]==0 or amplCols.shape[1]==0:
        return(1)
    # Stores number of permutations with scores more than minimal acceptable score
    score1=1
    # Stores number of performed permutations
    permNum=0
    for i in range(permNum1): 
        np.random.shuffle(patRow)
        # Score for the current permutation
        score=0
        k=0 # ???
        # Number of regions studied that have value more
        # than threshold for insertion or deletion
        moreThanThreshNum=0
        # Go through amplicon indices that are included into the potential CNV
        for j in range(min(len(patRow),
                           len(amplNums))):
            # If the current value is more than threshold and type of CNV is deletion
            if patRow[k+j]<=args.delTh2 and inDel==0:
                score+=4-patRow[k+j]
                if patRow[k+j]>args.delTh1:
                    moreThanThreshNum+=1            
            # If the current value is more than the threshold and type of CNV is insertion
            elif patRow[k+j]>=args.duplTh2 and inDel==1:
                score+=patRow[k+j]
                if patRow[k+j]<args.duplTh1:
                    moreThanThreshNum+=1
        # Divide the obtained score on the values depending on the length of CNV
        score/=args.delTh1**(moreThanThreshNum/len(amplNums))
        if score>=minScore:
            score1+=1
    score2=0
    for i in range(amplCols.shape[0]):
        score=0
        moreThanThreshNum=0
        for j in range(min(len(patRow),
                           amplCols.shape[1])):
            if amplCols[i,j]<=args.delTh2 and inDel==0:
                score+=4-amplCols[i,j]
                if amplCols[i,j]>args.delTh1:
                    moreThanThreshNum+=1
            elif amplCols[i,j]>=args.duplTh2 and inDel==1:
                score+=amplCols[i,j]
                if amplCols[i,j]<args.duplTh1:
                    moreThanThreshNum+=1
        score/=args.delTh1**(moreThanThreshNum/len(amplNums))
        if score>=minScore:
            score2+=1
    try:
        pValue=min(1,score1/permNum1)+score2/initAmplColSize
    except ZeroDivisionError:
        print('ERROR (3)! Some value in divider is zero:')
        print('permNum1:',permNum1)
        print('amplCols:',amplCols)
        print('amplCols shape:',amplCols.shape)
        exit(3)
    return(pValue)

# Converts some list to list of ranges
def to_ranges(iterable):
    iterable = sorted(set(iterable))
    keyfunc = lambda t: t[1] - t[0]
    for key, group in itertools.groupby(enumerate(iterable), keyfunc):
        group = list(group)
        yield group[0][1],group[-1][1]
