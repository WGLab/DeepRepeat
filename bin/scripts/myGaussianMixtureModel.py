
import os;
import sys;
import string;
import math;

import logging

import numpy as np
from sklearn.mixture import GaussianMixture as GM #GMM

from . import myHeader

def getLargerOverSmallRatio(p1, p2, commonoptions):
	mp = [int(p1),int(p2)]; #

	if math.fabs(mp[0]-mp[1])>=5:
		mratio = 3**(-mp[1]/float(mp[0]))
	else:
		if abs(mp[0]-mp[1])==4: mratio=0.50;
		elif abs(mp[0]-mp[1])==3: mratio=0.65;
		elif abs(mp[0]-mp[1])==2: mratio=0.80;
		else: mratio=0.95

	if mratio<0.001: mratio=0.001;

	if (not commonoptions==None) and ('pcr' in commonoptions) and commonoptions['pcr']==1:
		mratio = 0.8
		mratio = 0.5

	return mratio

def getRatioInfo(x_index2a, x_index2b, a_sum, b_sum, commonoptions):
	cur_ratio = round(b_sum/float(a_sum),3)
	cur_ratio_threshold = getLargerOverSmallRatio(x_index2a, x_index2b, commonoptions)

	return [cur_ratio, cur_ratio_threshold]

def getWindowForCounts(cur_point):
	cur_window = int(cur_point/200.0+0.75) + 1
	return cur_window

def getNeighbors_reads_fixed(lendict, cur_point, cur_window):
	cur_sum_read = 0; total_k = cur_window*2+1; has_k = 0;
	for curk in range(cur_point-cur_window, cur_point+cur_window+1):
		if curk in lendict:
			cur_sum_read += lendict[curk]
			has_k += 1;
	return [cur_sum_read, has_k, total_k]

def getNeighbors_reads(lendict, cur_point, minreads):
	cur_window = getWindowForCounts(cur_point) #int(cur_point/200.0+0.75)
	cur_sum_read = 0; total_k = cur_window*2+1; has_k = 0;
	for curk in range(cur_point-cur_window, cur_point+cur_window+1):
		if curk in lendict:
			cur_sum_read += lendict[curk]
			has_k += 1;
	#print cur_point, [cur_sum_read, has_k, total_k]
	
	#if cur_window==0:
	#	for i in [-1, 1]:
	#		if lendict.has_key(cur_point+i): 
	#			print cur_point+i, lendict[cur_point+i], cur_sum_read
	#			if lendict[cur_point+i]>lendict[cur_point]: 
	#				if lendict[cur_point+i]<cur_sum_read:
	#					print 'Error!!!! ', lendict[cur_point+i], cur_sum_read, cur_point, i
	#				cur_sum_read = lendict[cur_point+i]
	
	return [cur_sum_read, has_k, total_k]

def checkSmallSupport1(msup, lendict, minreads):
	cursum = getNeighbors_reads(lendict, msup, minreads)
	if cursum[0]/2<minreads:
		return False;
	else: return True;

def checkSmallSupport(peak2, lendict, minreads):
	newpeak2 = []
	for p in peak2:
		if checkSmallSupport1(p, lendict, minreads):
			newpeak2.append(p)
	return newpeak2

def getNeighbors2(x_index2a, x_index2b, lendict, minreads):
	neighbor3 = [0,0, 1,1,1,1]; point3 = [x_index2a, x_index2b]
	for pi in range(len(point3)):
		neighbor3[pi], neighbor3[pi+2], neighbor3[pi+4] = getNeighbors_reads(lendict, point3[pi], minreads)

	return neighbor3

def getClose(lendict, cur_point, minreads):
   cur_window = getWindowForCounts(cur_point) #int(cur_point/200.0+0.75)
   cur_max_red = 0;
   for curk in range(cur_point-cur_window, cur_point+cur_window+1):
      if curk in lendict:
         if lendict[curk]>cur_max_red: cur_max_red = lendict[curk]
   return cur_max_red

def selectFromTwoX(x_index2a, x_index2b, lendict, minreads, commonoptions):
	if x_index2a>x_index2b:
		x_index2a, x_index2b = x_index2b, x_index2a

	neighbor3 = getNeighbors2(x_index2a, x_index2b, lendict, minreads)

	msmall1 = False;
	if checkSmallSupport1(x_index2a, lendict, minreads): msmall1 = True;
	msmall2 = False;
	if checkSmallSupport1(x_index2b, lendict, minreads): msmall2 = True;

	if msmall1 and msmall2: #getRatioInfo(x_index2a, x_index2b, a_sum, b_sum):
		ratioInfo = getRatioInfo(x_index2a, x_index2b, neighbor3[0], neighbor3[1], commonoptions)
		risingle = getRatioInfo(x_index2a, x_index2b, getClose(lendict, x_index2a, minreads), getClose(lendict, x_index2b, minreads), commonoptions)
		#print ("pcrtest", x_index2a, x_index2b, ratioInfo, neighbor3)
		#print ("pcrtest", risingle, getClose(lendict, x_index2a, minreads), getClose(lendict, x_index2b, minreads))
		#if ratioInfo[0]<ratioInfo[1] and risingle[0]<risingle[1]:
		#if ratioInfo[0]<ratioInfo[1]:
		if ((not commonoptions==None) and ('pcr' in commonoptions) and commonoptions['pcr']==1):
			if ratioInfo[0]<ratioInfo[1] or risingle[0]<risingle[1]:
		  		x_index2 = x_index2a
			else: x_index2 = x_index2b
		else:
			if ratioInfo[0]<ratioInfo[1] and risingle[0]<risingle[1]:
				x_index2 = x_index2a
			else: x_index2 = x_index2b
	else:
		if msmall1: x_index2 = x_index2a
		elif msmall2: x_index2 = x_index2b
		else:
			if x_index2a<x_index2b:  x_index2 = x_index2b
			else: x_index2 = x_index2a

	return x_index2


def getNewRepeatForIllumina(lendict, MinSup, commonoptions, oldp2):
	if not ( len(oldp2)==0 or (len(oldp2)==1 and (not oldp2[0]==0)) or (len(oldp2)>1 and oldp2[0]==oldp2[1] and (not oldp2[0]==0)) ): return oldp2

	newrepeatsdictkeys = sorted(list(lendict.keys()))
	max1 = [0,0]; max2 = [0,0]
	allreads = 0;
	for nrkey in newrepeatsdictkeys:
		allreads = allreads + lendict[nrkey]
		if lendict[nrkey]>=MinSup:
			if lendict[nrkey]>max1[1]:
				max2 = [max1[0], max1[1]]
				max1 = [nrkey, lendict[nrkey]]
			else:
				if lendict[nrkey]>max2[1]:
					max2 = [nrkey, lendict[nrkey]]
	if max1[1]<MinSup:
		if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_WARNING: print('Warning!!! less than '+str(MinSup), max1, max2)
	if commonOptions['outlog'] <= myHeader.M_INFO: print(max1, max2, max2[1]/float(max1[1]), (max2[1]+max1[1])/float(allreads))
	if max2[1]==0:
		newrepeats = [max1[0], max1[0]]
	else:
		if max2[1]/float(max1[1])>0.6 or (max2[1]+max1[1])/float(allreads)>0.8:
			newrepeats = [max1[0], max2[0]]
		else: newrepeats = [max1[0], max1[0]]

	newrepeats = sorted(newrepeats)

	if len(oldp2)>0 and (oldp2[0] not in newrepeats):
		if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_WARNING: print('Warning!!! not in ', oldp2, newrepeats)

	return newrepeats

def myGMM(lendict, MinSup=2, mparameters=None, commonoptions=None):
	if len(lendict)<1: return [[0,0], []]
	elif len(lendict)==1: return [[lendict[lendict.keys()[0]], lendict[lendict.keys()[0]]], []];

	isillumina = False;
	if (not commonoptions==None) and 'SeqTech' in commonoptions and commonoptions['SeqTech']=="Illumina":
		isillumina = True;

	truecounts = None;
	if (not commonoptions==None) and 'truecounts' in commonoptions:
		truecounts = commonoptions['truecounts']

	minreads = int(MinSup);
	if minreads<2: minreads = 2

	ldkeys = sorted(list(lendict.keys()))
	allocr = 'allocr:';
	for ldk in ldkeys:
		allocr += ('%d:%d, ' % (ldk, lendict[ldk]))
	logging.info(allocr)
	if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_INFO: print (allocr, minreads, MinSup)
	#print (allocr)
	minrepcount = 5;
	for ldk in ldkeys:
		if ldk<minrepcount: del lendict[ldk]
	ldkeys = sorted(list(lendict.keys()))

	ldkeys = sorted(list(lendict.keys()))
	while len(ldkeys)>1:
		lastk = ldkeys[-1]; secondlk=ldkeys[-2];
		curw = getWindowForCounts(secondlk)
		if lendict[lastk]<minreads and lastk-secondlk>curw*10:
			del lendict[lastk]
		else: break;
		ldkeys = sorted(list(lendict.keys()))
	while len(ldkeys)>2:
		firstk = ldkeys[0]; secondk = ldkeys[1];
		curw = getWindowForCounts(secondk)
		if lendict[firstk]<minreads and secondk-firstk>curw*10:
			del lendict[firstk]
		else: break;
		ldkeys = sorted(list(lendict.keys()))
	#print(ldkeys)
	#for special case;
	if len(ldkeys)<1: return [[0], allocr]
	elif len(ldkeys)<2: peak2 = [ldkeys[0]];
	elif len(ldkeys)==2 or len(ldkeys)==3 or isillumina:
		maxk = ldkeys[0]
		if maxk==0: maxk = ldkeys[1]
		for ldk in ldkeys:
			if lendict[ldk]>=lendict[maxk] and ldk>0: maxk = ldk
		peak2 = [maxk]
	if len(ldkeys)<4 or isillumina:
		peak2 = checkSmallSupport(peak2, lendict, minreads)
		if len(peak2)==0: peak2 = [0]

		if isillumina:
			peak2 = getNewRepeatForIllumina(lendict, minreads, commonoptions, peak2)
		peak2 = sorted(peak2)
		return [peak2, allocr[:-1]];

	total_point = 0;
	mkeys = sorted(list(lendict.keys())); 
	for mk in mkeys:
		curnum = lendict[mk] - minreads;
		if curnum<1 or mk<minrepcount: 
			continue;
		total_point += curnum
	lowcov = False;
	if total_point<50:
		total_point = 0;
		lowcov = True;
		for mk in mkeys:
			curnum = lendict[mk]
			if curnum<1 or mk<minrepcount:
				continue;
			total_point += curnum

	X = np.zeros((total_point,1))
	xi = 0;
	for mk in mkeys:
		if not lowcov:
			curnum = lendict[mk] - minreads;
		else:
			curnum = lendict[mk]
		#if curnum<1 or mk<minrepcount:
		if lendict[mk]<2 or mk<minrepcount: 
			continue;
		for j in range(curnum):
			X[xi][0] = mk; xi += 1;
	#f len(X)<200: print total_point, lowcov, len(X), X

	default_n_components = [4] #[4,3,2,5,6,4];
	#for nc in range(2, 7):
	for nc in range(3, 7):
		if nc>=total_point: break;
		if nc==4: continue
		default_n_components.append(nc)
	default_n_components.append(4)
	gmm_times = 20;
	small_covars_threhold = 0.01
	if len(X)<3: return [[0,0], allocr[:-1]] 
	for cur_n_component_ind in range(len(default_n_components)):
		cur_n_component = default_n_components[cur_n_component_ind]
		
		atedge = True;
		for run_time in range(gmm_times):
			N = np.arange(1, (cur_n_component+1))
			models = [None for i in range(len(N))]

			for i in range(len(N)):
				models[i] = GM(N[i]).fit(X) 

			# compute the AIC and the BIC
			AIC = [m.aic(X) for m in models]
			calbic = False;
			if calbic:
				BIC = [m.bic(X) for m in models]

			mbest = models[np.argmin(AIC)]

			if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_DEBUG:
				print ('aic', np.argmin(AIC)) 
				for aic in range(len(AIC)):
					if aic==np.argmin(AIC):
						print ('<%.3f>' % (AIC[aic]))
					else:	print (' %.3f ' % (AIC[aic]))
				print ('')
				if calbic:
					print ('bic', np.argmin(BIC))
					for bic in range(len(BIC)):
						if bic==np.argmin(BIC):
							print ('<%.3f>' % (BIC[bic]))
						else: print (' %.3f ' % (BIC[bic]))
					print ('')
			
			if 0<np.argmin(AIC)<len(AIC)-1: 
				atedge = False;
				break;
			
		has_too_small_std = False;
		for i in range(len(mbest.means_)):
			if mbest.covariances_[i,0][0]<small_covars_threhold:
				has_too_small_std = True;
		#if (not has_too_small_std) and (not atedge): break;
		if (not atedge): break;
		elif cur_n_component_ind==len(default_n_components)-1:
			if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_WARNING:
				print ('Warning!!!! could not find optimized model')
				logging.info('Warning!!!! could not find optimized model')
			
	#print mbest.covariances_
	mean_covars = []; 
	for i in range(len(mbest.means_)):
		curk = int(mbest.means_[i,0]+0.5) #75)
		if curk in lendict:
			if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_DEBUG: print ('>>', i, ('%9.3fm' % (mbest.means_[i,0])), ('%6d' % (lendict[curk])), ('%20.9fst' % (mbest.covariances_[i,0][0])))
			mean_covars.append([curk, mbest.means_[i,0], lendict[curk], mbest.covariances_[i,0][0]])
		else:
			closedif = sys.maxsize; closekey = -1;
			for mk in mkeys:
				if closedif>abs(mk-curk):
					closedif=abs(mk-curk)
					closekey = mk
			if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_DEBUG: print ('>>', i, ('%9.3fm' % (mbest.means_[i,0])), ('%6d' % (lendict[closekey])), ('%20.9fst' % (mbest.covariances_[i,0][0])), closekey)
			mean_covars.append([closekey, mbest.means_[i,0], lendict[closekey], mbest.covariances_[i,0][0]])

	fixed_boundarywidth = 50; close_ratio_threhold = 0.8
	if (not commonoptions==None) and ('pcr' in commonoptions) and commonoptions['pcr']==1:
		close_ratio_threhold = 0.4
	
	remove_larger_covar_smaller_means = []
	for i in range(len(mean_covars)):
		mean_covars[i].append(getNeighbors_reads(lendict, mean_covars[i][0], minreads)[0])
	for i in range(len(mean_covars)):
		should_remove = False;
		for j in range(len(mean_covars)):
			if i==j: continue;
			if mean_covars[j][3]<small_covars_threhold: continue;

			if mean_covars[i][3]<=mean_covars[j][3] and mean_covars[i][1]<mean_covars[j][1]: pass
			elif mean_covars[i][3]>mean_covars[j][3] and mean_covars[i][1]<mean_covars[j][1]:
				meandif = (mean_covars[j][1]-mean_covars[i][1])*3
				#if meandif>10: meandif = 10
				if meandif>20: meandif = 20
				if mean_covars[i][3]>mean_covars[j][3]*meandif: should_remove = True
				else:
					cur_window_i = getWindowForCounts(mean_covars[i][0]) #int(mean_covars[i][0]/200.0+0.75)
					cur_window_j = getWindowForCounts(mean_covars[j][0]) #int(mean_covars[j][0]/200.0+0.75)
					if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_INFO: print (should_remove, mean_covars[i][0], mean_covars[j][0], cur_window_i, cur_window_j, mean_covars[i][0], mean_covars[j][0], mean_covars[i][4], mean_covars[j][4], mean_covars[j][4]*close_ratio_threhold, abs(cur_window_i-cur_window_j)==1, mean_covars[i][4]<mean_covars[j][4]*close_ratio_threhold)
					
					if abs(cur_window_i-cur_window_j)==1 and abs(mean_covars[i][0]-mean_covars[j][0])<fixed_boundarywidth:
						newneighbors = getNeighbors_reads_fixed(lendict, mean_covars[i][0], cur_window_j)
						if newneighbors[0]<mean_covars[j][4]*close_ratio_threhold: should_remove = True;
					elif mean_covars[i][4]<mean_covars[j][4]*close_ratio_threhold:
						should_remove = True;
					if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_INFO: print (should_remove)
		if not should_remove:
			remove_larger_covar_smaller_means.append(mean_covars[i])

	furtherremove = []
	for i in range(len(remove_larger_covar_smaller_means)):
		should_remove = False;
		for j in range(len(remove_larger_covar_smaller_means)):
			if i==j: continue;
			if abs(remove_larger_covar_smaller_means[i][1]-remove_larger_covar_smaller_means[j][1])<1.5 and (remove_larger_covar_smaller_means[i][3]<small_covars_threhold or remove_larger_covar_smaller_means[j][3]<small_covars_threhold):
				if remove_larger_covar_smaller_means[i][3]<small_covars_threhold and (not remove_larger_covar_smaller_means[j][3]<small_covars_threhold):
					should_remove = True;
				elif (not remove_larger_covar_smaller_means[i][3]<small_covars_threhold) and (remove_larger_covar_smaller_means[j][3]<small_covars_threhold): pass
				else:
					cur_window_i = getWindowForCounts(remove_larger_covar_smaller_means[i][0])
					cur_window_j = getWindowForCounts(remove_larger_covar_smaller_means[j][0])
					if abs(cur_window_i-cur_window_j)==1 and abs(remove_larger_covar_smaller_means[i][0]-remove_larger_covar_smaller_means[j][0])<fixed_boundarywidth:
						if cur_window_i<cur_window_j:
							newneighbors = getNeighbors_reads_fixed(lendict, remove_larger_covar_smaller_means[i][0], cur_window_j)
							if newneighbors[0]<remove_larger_covar_smaller_means[j][4]: should_remove = True;
						else:
							newneighbors = getNeighbors_reads_fixed(lendict, remove_larger_covar_smaller_means[j][0], cur_window_i)
							if newneighbors[0]>remove_larger_covar_smaller_means[i][4]: should_remove = True;
					elif remove_larger_covar_smaller_means[i][4]<remove_larger_covar_smaller_means[j][4]: should_remove = True;
		if not should_remove:
			furtherremove.append(remove_larger_covar_smaller_means[i]);

	remove_larger_covar_smaller_means, furtherremove = furtherremove, remove_larger_covar_smaller_means

	if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_DEBUG:
		print ('keep')
		for i in range(len(remove_larger_covar_smaller_means)):
			print ('>>', i, ('%9.3fm' % (remove_larger_covar_smaller_means[i][1])), ('%6d' % (remove_larger_covar_smaller_means[i][2])), ('%20.9fst' % (remove_larger_covar_smaller_means[i][3])), mean_covars[i][4])

	peak2 = []; max_p1 = 0; max_p1_reads = 0; max_p1_1_r = 0;
	for i in range(len(remove_larger_covar_smaller_means)):
		cur_window_i = getWindowForCounts(remove_larger_covar_smaller_means[i][0]) #int(remove_larger_covar_smaller_means[i][0]/200.0+0.75)
		cur_window_max = getWindowForCounts(max_p1) #int(max_p1/200.0+0.75)
		if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_INFO: print (max_p1, max_p1_reads, cur_window_i, cur_window_max, remove_larger_covar_smaller_means[i][0],max_p1, '<', abs(cur_window_i-cur_window_max)==1, '>', remove_larger_covar_smaller_means[i][4]>max_p1_reads, remove_larger_covar_smaller_means[i])
		if abs(cur_window_i-cur_window_max)==1 and abs(remove_larger_covar_smaller_means[i][0]-max_p1)<fixed_boundarywidth:
			newneighbors = getNeighbors_reads_fixed(lendict, remove_larger_covar_smaller_means[i][0], cur_window_max) 
			if newneighbors[0] > max_p1_reads:
				max_p1 = remove_larger_covar_smaller_means[i][0]
				max_p1_reads = remove_larger_covar_smaller_means[i][4]
				max_p1_1_r = remove_larger_covar_smaller_means[i][2] 
		elif remove_larger_covar_smaller_means[i][4]>max_p1_reads:
			max_p1 = remove_larger_covar_smaller_means[i][0]
			max_p1_reads = remove_larger_covar_smaller_means[i][4]
			max_p1_1_r = remove_larger_covar_smaller_means[i][2]

	peak2.append(max_p1)
	secondpeak = {}
	for i in range(len(remove_larger_covar_smaller_means)):
		if remove_larger_covar_smaller_means[i][0]==max_p1: continue;
		if selectFromTwoX(max_p1, remove_larger_covar_smaller_means[i][0], lendict, minreads, commonoptions)==remove_larger_covar_smaller_means[i][0]:
			secondpeak[remove_larger_covar_smaller_means[i][0]] = []
		else:
			#if remove_larger_covar_smaller_means[i][0]<max_p1 and remove_larger_covar_smaller_means[i][4]>max_p1_reads*close_ratio_threhold:
			#	secondpeak[remove_larger_covar_smaller_means[i][0]] = []
			if remove_larger_covar_smaller_means[i][4]>max_p1_reads*close_ratio_threhold or remove_larger_covar_smaller_means[i][2]>max_p1_1_r*close_ratio_threhold:
				secondpeak[remove_larger_covar_smaller_means[i][0]] = []

	if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_INFO: print ('peak2', peak2, secondpeak)
	secondpeakkeys = secondpeak.keys();
	for spk in secondpeakkeys:
		for spkj in secondpeakkeys:	
			if spkj in [max_p1, spk]: continue;
			if selectFromTwoX(spk, spkj, lendict, minreads, commonoptions)==spk:
				if spk not in secondpeak[spkj]: secondpeak[spkj].append(spk)
			else:
				if spkj not in secondpeak[spk]: secondpeak[spk].append(spkj)

	for spk in secondpeakkeys:
		if  len(secondpeak[spk])==0: 
			peak2.append(spk)

	if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_INFO: print ('peak2', peak2, secondpeak	)
	if len(secondpeakkeys)>0 and len(peak2)<2:
		max_p2 = 0; max_p2_reads = 0;
		for i in range(len(remove_larger_covar_smaller_means)):
			if remove_larger_covar_smaller_means[i][0] in secondpeak:
				if remove_larger_covar_smaller_means[i][4]>max_p2_reads: 
					max_p2_reads = remove_larger_covar_smaller_means[i][4]
					max_p2 = remove_larger_covar_smaller_means[i][0]

		peak2.append(max_p2)
	if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_DEBUG: print ('peak2', peak2)
	peak2 = checkSmallSupport(peak2, lendict, minreads)
	if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_DEBUG: print ('peak2', peak2)
	
	if len(peak2)==2:
		if peak2[0]==0 and (not peak2[1]==0): peak2[0] = peak2[1]
		if peak2[1]==0 and (not peak2[0]==0): peak2[1] = peak2[0]
	if len(peak2)==1 or (len(peak2)>1 and abs(peak2[0]-peak2[1])>2):
		curlen = len(peak2)
		if curlen>2: peak2 = 2;
		for i in range(curlen):
			cur_window_max = getWindowForCounts(peak2[i])
			for npdif in range(1, cur_window_max+1):
				newpeak = peak2[i]+npdif
				if newpeak in lendict and lendict[newpeak]>1.5*lendict[peak2[i]] and lendict[newpeak]>5*minreads:
					peak2[i] = newpeak
				newpeak = peak2[i]-npdif
				if newpeak in lendict and lendict[newpeak]>1.5*lendict[peak2[i]] and lendict[newpeak]>5*minreads:
					peak2[i] = newpeak
			#for newpeak in range(peak2[i]-cur_window_max, peak2[i]+cur_window_max+1):
			#	if lendict.has_key(newpeak) and lendict[newpeak]>2*lendict[peak2[i]]:
			#		peak2[i] = newpeak
		if len(peak2)==1:
			peak2.append(peak2[0]);
	if len(peak2)==0: peak2 = [0,0]

	peak2 = sorted(peak2)

	if not truecounts==None:
		if abs(truecounts[0]-peak2[0])>5 or abs(truecounts[1]-peak2[1])>5:
			if (not commonoptions==None) and ('outlog' in commonoptions) and commonoptions['outlog'] <= myHeader.M_INFO: print ('Big dif, please check', peak2, truecounts)

	return [peak2, allocr[:-1]]

def get2Peaks(lengd, MinSup, commonoptions=None):
	lendict = {};
	for l in lengd:
		l = int(l+0.5)
		if l not in lendict: lendict[l] = 0;
		lendict[l] += 1;

	return myGMM(lendict, MinSup, commonoptions=commonoptions)


if __name__=='__main__':
   a = 0;
