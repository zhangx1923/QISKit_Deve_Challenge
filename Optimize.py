# -*- coding: utf-8 -*-
from methods import *
MAX_TIMES = 4
class Optimize:
	#the init function
	def __init__(self,dagcircuit,coupling_map):
		#the number of the qubits
		self.qNum = dagcircuit.width()
		self.cmap = coupling_map

		#interQRecord is used to record the interQ of two qubits 
		#which aren't connected with each other directly.
		#if self.__getInterQ and self.__computeCost, we will record value in it
		#the format of it is a two-dimen list, the first dimen is the [cq,tq]
		#and the second dimen is the interq of cq and tq. They have the same index
		self.interQRecord = [[],[]]

		#partMap is used to record the part id map of two qubit
		#the first dimen stands for the tow qubit:[fixed,swapped]
		#the second dimen stands for the map relation
		#restore value in self.__getPartMap
		self.partMap = [[],[]]

	########################the private method##############################

	#get the allowed cnot
	def __getAllow(self):
		allowCnot = []
		onewayCnot = []
		dcmap = copy.deepcopy(self.cmap)
		for k in self.cmap.keys():
			for i in range(0,len(self.cmap[k])):
				if int(self.cmap[k][i]) in dcmap.keys():
					if int(k) not in dcmap[int(self.cmap[k][i])]:
						dcmap[int(self.cmap[k][i])].append(int(k))
					else:
						pass
				else:
					dcmap[int(self.cmap[k][i])] = [int(k)]
				tmp = [int(k),int(self.cmap[k][i])]
				tmpReser = [int(self.cmap[k][i]),int(k)]
				if tmp not in allowCnot:
					allowCnot.append(tmp)
					onewayCnot.append(tmp)
				if tmpReser not in allowCnot:
					allowCnot.append(tmpReser)
		self.allowSets = allowCnot
		self.onewayCnot = onewayCnot
		self.dcmap = dcmap
		#return allowCnot

	#get the cnot connectivity of this QP and count the number of each cnot
	def __getCNOT(self,qasm = None):
		cnotSet = []

	    #record the whole record of the cnot
		cnotTimes = []
		if qasm == None:
			qasm = self.qasm
		for q in qasm.split(";\n"):
			m = re.match(r'cx q\[(\d+)\],q\[(\d+)\]',q)
			if m != None:
				#the id of the control qubit 
				cqID = m.group(1)
				#the id of the target qubit
				tqID = m.group(2)
				tmp = [int(cqID),int(tqID)]
				tmpReser = [int(tqID),int(cqID)]
				#if the record has been stored in the set, then pass it
				if tmp in cnotSet or tmpReser in cnotSet:
					pass
				else:
					cnotSet.append(tmp)
				cnotTimes.append(tmp)

		res = [cnotSet,cnotTimes]
		return res

	#get the illegal cnot
	def __getIll(self,cnotSets = None,allowSets = None):
		if cnotSets == None:
			cnotSets = self.cnotSets
		if allowSets == None:
			allowSets = self.allowSets
		bakSets = copy.deepcopy(cnotSets)
		for cnot in cnotSets:
			if cnot in allowSets:
				bakSets.remove(cnot)
		return bakSets

	#adjust the overall id by dynamic programming 
	def __adjustAllMap1(self,bc,cnotSet,amap,costs,mapRecord):
		q1 = bc[0]
		q2 = bc[1]
		canUse = {}
		canUse[q1] = self.dcmap[q2]
		canUse[q2] = self.dcmap[q1]		
		index = cnotSet.index(bc)
		alternativeMap = []
		for k in bc:
			for k in [q1,q2]:
				for q in canUse[k]:
					tmpMap = {k:q,q:k}
					map_bool = True
					for ind in range(0,index+1):
						cs = cnotSet[ind]
						if cs[0] in tmpMap.keys():
							csq1 = tmpMap[cs[0]]
						else:
							csq1 = cs[0]
						if cs[1] in tmpMap.keys():
							csq2 = tmpMap[cs[1]]
						else:
							csq2 = cs[1]
						if [csq1,csq2] not in self.allowSets:
							map_bool = False
							break
					if map_bool and tmpMap not in alternativeMap:
						alternativeMap.append(tmpMap)

		if alternativeMap == []:
			#we have no choice
			#compute the cost
			if amap in mapRecord:
				return True
			cost = self.__evaluateCost2(cnotSet)
			if costs != [] and cost >= min(costs):
				return True
			costs.append(cost)
			mapRecord.append(amap)
			return True
			
		for am in alternativeMap:
			tmpBC = []
			tmpCS = copy.deepcopy(cnotSet)
			tmpAM = copy.deepcopy(amap)

			#change the cnotSet and set the badCnot
			for ind in range(0,len(tmpCS)):
				if tmpCS[ind][0] in am.keys():
					c1 = am[tmpCS[ind][0]]
				else:
					c1 = tmpCS[ind][0]
				if tmpCS[ind][1] in am.keys():
					c2 = am[tmpCS[ind][1]]
				else:
					c2 = tmpCS[ind][1]
				if [c1,c2] not in self.allowSets:
					#bab cnot
					tmpBC.append([c1,c2])
				tmpCS[ind] = [c1,c2]	
			#compute the new amap
			deleteKey = []
			for key in am.keys():
				if key not in tmpAM.keys():
					tmpAM[key] = am[key]
				else:
					for k1 in tmpAM.keys():
						if tmpAM[k1] == key:
							if k1 == am[key]:
								deleteKey.append(k1)
							else:
								tmpAM[k1] = am[key]
							break
			for dk in deleteKey:
				del tmpAM[dk]
			#iterate
			if tmpBC == []:
				mapRecord.append(tmpAM)
				costs.append(0)
				return True
			else:
				self.__adjustAllMap1(tmpBC[0],tmpCS,tmpAM,costs,mapRecord)	

	#evaluate the cost of the cnot set in the qasm code
	#and the function will be called by __adjustAllMap1 and ccIterate
	#cnotSets
	def __evaluateCost(self,cnotSet):
		#argument:
		#cnotSet: the cnot sets
		cost = 0
		for csind in range(len(cnotSet)):
			cs = cnotSet[csind]
			#if cs in self.allowSets and cs not in self.onewayCnot:
				#we need to reverse the cnot
				#cost += 4
			if cs not in self.allowSets:
				#bad cnot
				cost += (len(self.__getInterQ(cs[0],cs[1]))*34)
		return cost

	#cnotTimes
	def __evaluateCost2(self,cnotSet):
		cost = 0
		if len(cnotSet) == 0:
			return 0
		lastCnot = copy.deepcopy(cnotSet[0])
		badTimes = -1
		badcnotTimes = 0
		badCost = []
		badDepth = []
		for csind in range(len(cnotSet)):
			cs = cnotSet[csind]
			if cs in self.allowSets and cs not in self.onewayCnot:
				badTimes += 1
				badCost.append(4)
				badDepth.append(len(cnotSet) - csind)
			elif cs not in self.allowSets:
				badcnotTimes += 1		
				if (lastCnot[0] == cs[0] and lastCnot[1] == cs[1]) or (lastCnot[1] == cs[0] and lastCnot[0] == cs[1]):
					#don't repeat computing the cost
					pass
				else:
					badTimes += 1
					iqs = self.__getInterQ(cs[0],cs[1])
					badCost.append(34 * len(iqs))
					badDepth.append(len(cnotSet) - csind)
			else:
				pass
			lastCnot = copy.deepcopy(cs)
		first = True
		for j in range(badTimes):
			tmpcost = 0
			if first:
				tmpcost = badCost[j]
				if tmpcost > 4:
					#bad CNOT
					first = False
			else:
				if badCost[j] == 4:
					continue
				tmpcost = badCost[j] * (badDepth[j]  / len(cnotSet)) * ((badTimes-j)/badTimes)
			cost += tmpcost
		return cost

	#get the inter qubit list of control-q and target-q
	def __getInterQ(self,q1,q2):
	    #argument:
	    #q1: the id of the control qubit
	    #q2: the id of the target qubit
		#print(self.interQRecord)
		if [q1,q2] in self.interQRecord[0]:
			index = self.interQRecord[0].index([q1,q2])
			return copy.deepcopy(self.interQRecord[1][index])
	    #the interQ is a one-dimen list which records the inter qubit, for example: [0]
	    #which means that the first bad cnot can connect with each other by qubit '0'
	    #or [1,2]
	    #it means that bad cnot can connect with each other by 'cq--2--1--tq'
		interQ = []

		#the parameter actually stands for a tree and the root node is the control qubit
		#interQs_cq[n] stands for the nth layer of the tree
		#[[1],[2,3][1,3,-1,4,5]] means the root node is 1 and the second layer is 2 and 3
		#-1 is a separator, we can ignore it and regard 1,3 as 2's children and 4,5 as 3's children
		interQs_cq = [[q1]]

		#"q_num" stands for the number of the inter qubit we should use to connect q1 and q2
		#and it must be no bigger than (total qubits - 2) 
		q_num = 0
		try:
			while q_num <= self.qNum-2 and interQ == []:
			    #the layer of the current inter qubit
				layer = interQs_cq[q_num]
				newLayer = []
				for j in range(0,len(layer)):
					q = layer[j]
			        #"-1" stands for the separator
					if q == -1:
						continue
					if j != 0:
						#add the separator
						newLayer.append(-1)						
					if q2 in self.dcmap[q]:
						#has found!
						aimQ = q
						for i in range(q_num,0,-1):
							interQ.append(aimQ)
							layer = interQs_cq[i]
					        #computer the pre-layer
							if i > 1:
								index = layer[0:layer.index(aimQ)+1].count(-1)
								#the index is not the actual index
								#because we should pass the -1
								tmpList = interQs_cq[i-1]
								tmp1 = 0
								tmp2 = -1
								for ii in range(0,len(tmpList)):
									if tmp2 == index:
										break
									if tmpList[ii] == -1:
										tmp1 += 1
									else:
										tmp2 += 1
								index += tmp1
								aimQ = tmpList[index]
						raise GetOutLoop()
					else:
						for element in self.dcmap[q]:
							if element in newLayer:
								continue
							newLayer.append(element)

				interQs_cq.append(newLayer)
				q_num += 1
		except GetOutLoop:
			pass
		self.__addRecord(q1,q2,interQ)
		return interQ 

	#add the interQ record to self.interQRecord
	def __addRecord(self,q1,q2,interQ):
		self.interQRecord[0].append([q1,q2])
		self.interQRecord[0].append([q2,q1])
		tmp = copy.deepcopy(interQ)
		self.interQRecord[1].append(tmp)
		tmp1 = copy.deepcopy(interQ)
		tmp1.reverse()
		self.interQRecord[1].append(tmp1)		

	#get the part map
	def __getPartMap(self,q1,q2,interQ):
		#starttime = datetime.datetime.now()
		#q1 is fixed qubit
		#q2 is swapped qubit
		partMap = {}
		cnot = [q1,q2]
		if cnot in self.partMap[0]:
			index = self.partMap[0].index([q1,q2])
			return self.partMap[1][index]
		for iq in interQ:
			iq_bool = False
			tq_bool = False
			for k in partMap.keys():
				swapElement = -1
				if k == iq:
					iq_bool = True
					swapElement = q2
				elif k == q2:
					tq_bool = True
					swapElement = iq
				else:
					continue
				value = partMap[k]
				key = -1
				while value != k:
					key = value
					value = partMap[value]
				partMap[key] = swapElement						
			if iq_bool == False:
				partMap[iq] = q2
			if tq_bool == False:
				partMap[q2] = iq
			q2 = iq	
		
		#restore the partMap
		self.partMap[0].append(cnot)
		self.partMap[1].append(partMap)
		#endtime = datetime.datetime.now()
		#print((endtime - starttime).microseconds)
		#return the value
		return partMap	

	##########################the greedy method################################
	#compute the cost of adjusting the bad CNOT
	def __computeCost1(self,q1,q2,interQ,qasm):
		#Argument:
		#q1: the id of the control qubit
		#q2: the id of the target qubit
		#interQ: the inter qubit list of q1 and q2
		#qasm: the remained qasm list, and the type is a list. we don't need to split the code again
		#Return:[a,b,c]
		#a stands for the id of the qubit which is to be swapped
		#b stands for the remained qasm code
		#c stands for the cnots set in the remained qasm code
		costs = []
		exeRecord = []
		#the currnet depth
		depth = 0

		cost = 0
		cRecord = [[],[]]
		self.__ccIterate(q1,q2,interQ,qasm,0,depth,costs,exeRecord,cost,cRecord)
		minIndex = costs.index(min(costs))
		return exeRecord[minIndex]

		#use random to decide
		# costSum = 0
		# times = []
		# for i in costs:
		# 	costSum += i
		# 	times.append(0)
		# costtmp = [i / costSum for i in costs]
		# costN = copy.deepcopy(costs)
		# for j in range(len(costs)):
		# 	tmp = 0
		# 	for k in range(j+1):
		# 		tmp += costtmp[k]
		# 	costN[j] = tmp
		# import random
		# for i in range(len(costs)*5):
		# 	a = random.random()
		# 	for l in range(len(costs)):
		# 		if a <= costN[l]:
		# 			times[l] += 1
		# 			break
		# maxIndex = times.index(max(times))
		# return exeRecord[maxIndex]

	def __ccIterate(self,q1,q2,interQ,qasm,lines,depth,costs,exeRecord,cost,cr):
		######################## greedy ##########################
		#compute the part map of the two case:
		#1.swap the q1 and the interQ
		#2.swap the q2 and the interQ

		cost += 34 * len(interQ)

		for i in range(0,2):
			tmpInterQ = copy.deepcopy(interQ)
			tmpQASM = copy.deepcopy(qasm)
			#tmpBadCnot = []
			if i == 0:
				tmpInterQ.reverse()
			cnot = [q1,q2]
			tmpPartMap = {}
			#set the value of the partMaps			
			for iq in tmpInterQ:
				iq_bool = False
				tq_bool = False
				for k in tmpPartMap.keys():
					swapElement = -1
					if k == iq:
						iq_bool = True
						swapElement = cnot[i] 
					elif k == cnot[i]:
						tq_bool = True
						swapElement = iq
					else:
						continue
					value = tmpPartMap[k]
					key = -1
					while value != k:
						key = value
						value = tmpPartMap[value]
					tmpPartMap[key] = swapElement						
				if iq_bool == False:
					tmpPartMap[iq] = cnot[i]
				if tq_bool == False:
					tmpPartMap[cnot[i]] = iq
				cnot[i] = iq	

			#modify the qasm code and set the contSets
			badLine = -1
			badCnot = []
			TmpCnotSet = []
			measureMap = {}
			for j in range(lines,len(tmpQASM)):
				q = tmpQASM[j]
				qids = re.findall(r'\[(\d+)\]',q)
				elementStr = re.split(r'\[\d+\]',q)
				if 'measure' in q:
					#the measurement statement
					measureMap[int(qids[0])] = int(qids[1])
					continue
				elif q == "":
					#the blank statement
					continue
				else:
					#the common statement
					tmpQASMcode = elementStr[0]
					for t in range(0,len(qids)):
						try:
							qid = str(tmpPartMap[int(qids[t])])
						except KeyError:
							qid = str(qids[t])
						tmpQASMcode += "["
						tmpQASMcode += qid
						tmpQASMcode += "]"
						tmpQASMcode += elementStr[t+1]
					tmpQASM[j] = tmpQASMcode

				#judge whether the current statement is the cnot statement
				m = re.match(r'cx q\[(\d+)\],q\[(\d+)\]',tmpQASM[j])	
				if m != None:
					TmpCnotSet.append([int(m.group(1)),int(m.group(2))])
					if badLine == -1:
						if [int(m.group(1)),int(m.group(2))] in self.allowSets:
							if [int(m.group(1)),int(m.group(2))] in self.onewayCnot:
								pass
							else:
								cost += 4
						else:
							#find the bad cnot
							badLine = j
							badCnot = [int(m.group(1)),int(m.group(2))]
			#add the measure statement
			mmBak = copy.deepcopy(measureMap)
			for key in tmpPartMap.keys():
				value = tmpPartMap[key]
				measureMap[value] = mmBak[key]

			#the final statement is a blank
			tmp = 1
			total = len(tmpQASM)
			for mk in measureMap.keys():
				tmpQASM[total-1-tmp] = "measure q["+ str(mk) +"] -> c["+ str(measureMap[mk]) +"]"
				tmp += 1	
			#restore value to cr and executeRecord		
			tmpcr = copy.deepcopy(cr)
			tmpcr[1].append(tmpQASM)
			#tmpcr[2].append(tmpBadCnot)
			tmpCost = cost
			if i == 0:
				tmpCost += 4
				tmpcr[0].append(q1)
			else:
				tmpcr[0].append(q2)
			if badLine == -1:
				#no bad cnot
				costs.append(tmpCost)
				exeRecord.append(tmpcr)
				return True
			elif depth == MAX_TIMES-1:
				tmpCost += self.__evaluateCost2(TmpCnotSet)
				costs.append(tmpCost)
				exeRecord.append(tmpcr)
				return True			
			else:
				tmpInterQ = self.__getInterQ(badCnot[0],badCnot[1])
				self.__ccIterate(badCnot[0],badCnot[1],tmpInterQ,tmpQASM,badLine,depth + 1,costs,exeRecord,tmpCost,tmpcr)

	#use swap gates to make bad cnot feasible 
	def __adjustPartMap1(self,interQList = None):
	    #Argument:
		#interQList: the inter qubit list
		newQasm = ""
		qasmList = self.qasm.split(";\n")
		lines = len(qasmList)
		badCnot = self.badCnot

		for ind in range(0,lines):
			q = qasmList[ind]
			m = re.match(r'cx q\[(\d+)\],q\[(\d+)\]',q)
			if m != None:
				#cnot
				cq = int(m.group(1))
				tq = int(m.group(2))
				#print(badCnot)
				if [cq,tq] in badCnot or [tq,cq] in badCnot:
					#the current cnot is the bad cnot
					pass

				else:
					#the cnot statement is allowed in current layout
					newQasm += q
					newQasm += ";\n"
					continue
				
				interQ = self.__getInterQ(cq,tq)
				#judge whether we should use cq tp swap or tq to swap
				res = self.__computeCost1(cq,tq,interQ,qasmList[ind:lines])
				#change the following qasm code
				
				qasmList[ind:lines] = res[1][0]
				swapQID = res[0][0]
				fixedQID = cq if swapQID == tq else tq

				#get the new bad cnots set
				badCnot = []
				for qq in qasmList[ind:lines]:
					m = re.match(r'cx q\[(\d+)\],q\[(\d+)\]',qq)
					if m != None and [int(m.group(1)),int(m.group(2))] not in self.allowSets:
						#badCnot
						badCnot.append([int(m.group(1)),int(m.group(2))])
				if fixedQID != cq:
					interQ.reverse()

				for iq in interQ:
		            #the code depends on one hypothesis that there is no isolated qubit
					newQasm += self.__getSwapCode(iq,swapQID)
					swapQID = iq

				#add a new cnot statement 
				if fixedQID != cq:
					newQasm += "u2(0.0,3.14159265359) q[" + str(fixedQID) + "];\n"
					newQasm += "u2(0.0,3.14159265359) q[" + str(interQ[-1]) + "];\n"
				newQasm += "cx q[" + str(fixedQID) + "],q[" + str(interQ[-1]) + "];\n"
				if fixedQID != cq:
					newQasm += "u2(0.0,3.14159265359) q[" + str(fixedQID) + "];\n"
					newQasm += "u2(0.0,3.14159265359) q[" + str(interQ[-1]) + "];\n"					
				#newQasm += newcnot
		        ############act the SWAP gate again to restore the state of the qubit######### 
				# interQ.reverse()
				# if res == None:
				# 	for iq in range(0,len(interQ)-1):
				# 		newQasm += self.__getSwapCode(interQ[iq],interQ[iq+1])
				# newQasm += self.__getSwapCode(interQ[-1],bak_q)
		        ############################end of this method################################
			elif q == "":
				pass

			else:   
				newQasm += q
				newQasm += ";\n"     

		#choose the appropriate direction of the swap gate
		#so that we can use less H gates to adjust the direction of the elemental cnot gates
		#print(newQasm)
		return newQasm
	##########################the greedy method################################

	def __judgeReverse(self,qid1,qid2):
	    #Argument:
	    #qid1 & qid2: judge whether "cx q[qid1],q[qid2]" can be executed by add H to reverse the direction
	    
		if (qid1 not in self.cmap.keys()) or (qid1 in self.cmap.keys() and qid2 not in self.cmap[qid1]): 
			if qid2 in self.cmap.keys() and qid1 in self.cmap[qid2]:
				return True
		return False	   

	#get the code of the swap gate     
	def __getSwapCode(self,qid1,qid2):
	    #Argument:
	    #qid1 & qid2: the two qubits need to swap

		if qid1 in self.cmap.keys() and qid2 in self.cmap[qid1]:
			pass
		else:
			tmp = qid1
			qid1 = qid2
			qid2 = tmp 

		swap = "cx q[" + str(qid1) + "],q[" + str(qid2) + "];\n"
		#cnot qid2,qid1
		swap += "u2(0.0,3.14159265359) q[" + str(qid1) + "];\n"
		swap += "u2(0.0,3.14159265359) q[" + str(qid2) + "];\n"
		swap += "cx q[" + str(qid1) + "],q[" + str(qid2) + "];\n"
		swap += "u2(0.0,3.14159265359) q[" + str(qid1) + "];\n"
		swap += "u2(0.0,3.14159265359) q[" + str(qid2) + "];\n"
		swap += "cx q[" + str(qid1) + "],q[" + str(qid2) + "];\n"
		return swap

	#change the qubit's id according to the map
	#the type of qasm is string
	#the funciton will be called by self.generateQASM only 
	def __changeID(self,qasm,maps):
		#Argument:
		#qasm: QASM code:string
		#maps: a dict, and the key is the original id, the value is new id.
		newQASM = ""

		#if the maps is none, which means that the id of qubits has not been modified
		if maps == {}:
			return qasm

		measureMap = {}
		for q in qasm.split(";\n"):
			qids = re.findall(r'\[(\d+)\]',q)
			elementStr = re.split(r'\[\d+\]',q)
			tmpQASM = elementStr[0]
			if q == "":
				#the final statement
				continue
			elif 'measure' in q:
				#the statement is "measurement"
				measureMap[int(qids[0])] = int(qids[1])
				continue
			else:
				for t in range(0,len(qids)):
					try:
						qid = str(maps[int(qids[t])])
					except KeyError:
						qid = str(qids[t])
					tmpQASM += "["
					tmpQASM += qid
					tmpQASM += "]"
					tmpQASM += elementStr[t+1]
			newQASM += tmpQASM
			newQASM += ";\n"
		#add the measure statement
		mmBak = copy.deepcopy(measureMap)
		for k in maps.keys():
			value = maps[k]
			measureMap[value] = mmBak[k]

		for mk in measureMap.keys():
			newQASM += "measure q["+ str(mk) +"] -> c["+ str(measureMap[mk]) +"];\n"

		return newQASM

	#the type of qasm is list
	#the function will be called by self.__adjustPartMap1 and __adjustPartMap2
	def __changeID2(self,qasm,maps):
		measureMap = {}
		badCnot = []
		for j in range(0,len(qasm)):
			q = qasm[j]
			qids = re.findall(r'\[(\d+)\]',q)
			elementStr = re.split(r'\[\d+\]',q)
			if 'measure' in q:
				#the measurement statement
				measureMap[int(qids[0])] = int(qids[1])
				continue
			elif q == "":
				#the blank statement
				continue
			else:
				#the common statement
				tmpQASM = elementStr[0]
				for t in range(0,len(qids)):
					try:
						qid = str(maps[int(qids[t])])
					except KeyError:
						qid = str(qids[t])
					tmpQASM += "["
					tmpQASM += qid
					tmpQASM += "]"
					tmpQASM += elementStr[t+1]
				qasm[j] = tmpQASM

			#judge whether the current statement is the cnot statement
			m = re.match(r'cx q\[(\d+)\],q\[(\d+)\]',qasm[j])	
			if m != None:
				bad = [int(m.group(1)),int(m.group(2))]
				badReverse = [int(m.group(2)),int(m.group(1))]
				if bad not in self.allowSets:
					if bad in badCnot or badReverse in badCnot:
						pass
					else: 
						badCnot.append(bad)

		#add the measure statement
		mmBak = copy.deepcopy(measureMap)
		for key in maps.keys():
			value = maps[key]
			measureMap[value] = mmBak[key]

		#the final statement is a blank
		tmp = 1
		total = len(qasm)
		for mk in measureMap.keys():
			qasm[total-1-tmp] = "measure q["+ str(mk) +"] -> c["+ str(measureMap[mk]) +"]"
			tmp += 1

		return [badCnot,qasm]

	def __computeAllCost(self):
		cost = 0
		for cs in self.cnotTimes:
			if cs not in self.allowSets:
				#badCnot
				cost += len(self.__getInterQ(cs[0],cs[1]))
		return cost
	#################end of private method##################################

	###########################the public method############################

	#set the value of self.value and call other functions 
	def setQASM(self,qasm):
		self.qasm = qasm
		#set the value of self.allowSets and self.dcmap
		#self.dcmap stands for the full coupling map
		self.__getAllow()
		#set the value of self.cnotSets and self.cnotTimes
		cnotRecord = self.__getCNOT()
		self.cnotSets = cnotRecord[0]
		self.cnotTimes = cnotRecord[1]
    	#the parameter records the bad cnot, which means that the cnot can't be supplied by the layout
    	#but we can make the bad cnot possible by adjusting the id or using SWAP gates
		self.badCnot = self.__getIll()
		# print(self.cnotTimes)
		# print("\n")
		#print(self.cnotSets)
		# print("\n")
		#print(self.cmap)
		#print("!!!!!!!!!!!!!!!!!!")
		# print("\n")
		# print(self.badCnot)

	#adjust the id of qubits so that it can satisfies the requirement of the device    
	def adjustID(self):
		#idMap = self.__adjustAllMap()
		costs = []
		mapRecord = []
		css = copy.deepcopy(self.cnotSets)
		self.__adjustAllMap1(self.badCnot[0],css,{},costs,mapRecord)
		mincost = min(costs)
		idMap = mapRecord[costs.index(mincost)]
		if idMap != {}:
			ss = self.__changeID(self.qasm,idMap)
			self.setQASM(ss)
		#make the impossible cnot gate possible with the help of SWAP
		res = self.__adjustPartMap1()
		return res

	#reduce gates if there are two single gates act on one qubit continuously
	def reduceSameGate(self,qasm):
		newQASM = ""
		measureQASM = ""
		qQASM = []
		for i in range(0,self.qNum):
			tmp = [[]]
			qQASM.append(tmp)
		cnotList = []
		ql = qasm.split("\n")
		for i in range(0,len(ql)):
			q = ql[i]
			m1 = re.match(r'u(\d)\((.+)\) q\[(\d+)\];',q)
			m2 = re.match(r'cx q\[(\d+)\],q\[(\d+)\];',q)
			if m1 != None:
				#single qubit gate statement
				qid = m1.group(3)
				#change u1 and u2 to u3
				qQASM[int(qid)][-1].append(q)
			elif m2 != None:
				q1 = int(m2.group(1))
				q2 = int(m2.group(2))
				cnot = [q1,q2]
				cnotList.append(cnot)
				qQASM[q1].append([])
				qQASM[q2].append([])
			elif q == "":
				pass
			elif 'measure' in q or "barrier" in q:
				measureQASM += q
				measureQASM += '\n'
			else:
				newQASM += q
				newQASM += "\n"

		#add the remanent statement
		cnotList.append([m for m in range(0,self.qNum)])
		for c in cnotList:
			for qq in c:
				statement = ""
				num = len(qQASM[qq][0])
				if num == 0:
					pass
				elif num == 1:
					statement  = qQASM[qq][0][0]
					statement += "\n"
				else:
					#find the statement can be merged
					statement = mergeSingleGate(qQASM[qq][0])
				newQASM += statement
				del qQASM[qq][0]
			if len(c) == 2:
				newQASM += "cx q[" + str(c[0]) + "],q[" + str(c[1]) + "];\n"
		newQASM += measureQASM
		return newQASM

	#generate the qasm code accorrding to the qasm list
	#and in the process of generating, adjust the direction of special cnot
	def generateQASM(self,qasm,idMap=None):
	    #argument:
	    #qasm: the qasm code
	    #idMap: a dict, and the key is the original id, the value is new id.

		newQASM = ""

		if idMap == None:
			idMap = {}
			for i in range(0,self.qNum):
				idMap[i] = i

		newQASM = self.__changeID(qasm,idMap)
		result = ""
		insertGateMap = False

		for q in newQASM.split("\n"):
			#generate the code
			m = re.match(r'cx q\[(\d+)\],q\[(\d+)\];',q)
			if m != None and self.__judgeReverse(int(m.group(1)),int(m.group(2))):
				id1 = str(m.group(2))
				id2 = str(m.group(1))
				#add H gate to reverse the direction of the cnot gate
				tmpStr = "u2(0.0,3.14159265359) q[" + id1 + "];\n"
				tmpStr += "u2(0.0,3.14159265359) q[" + id2 + "];\n"
				tmpStr += "cx q[" + id1 + "],q[" + id2 + "];\n"
				tmpStr += "u2(0.0,3.14159265359) q[" + id1 + "];\n"
				tmpStr += "u2(0.0,3.14159265359) q[" + id2 + "];\n"                
				result += tmpStr
			elif q == "}" and insertGateMap == False:
				result += q
				result += "\n"
				result += "gate u2(phi,lambda) q\n{\n  U((pi/2),phi,lambda) q;\n}\n"
				insertGateMap = True
			elif q == "":
				pass
			else:
				result += q
				result += "\n"

		return result

	#######################end of public method#############################



#the exception class is used to getting out of a multi-loop
class GetOutLoop(Exception):
    pass


	# #choose the feasible solution by backtracing
	# def __backTrace(self,depth,solution,ss):
	# 	#Argument:
	# 	#depth: the current depth of the solution space tree
	# 	#N: the max depth of the solution space
	# 	#solution: the current solution, which is a dict
	# 	#cnots: the cnot connectivity in the QP
	# 	#allowcnots: the cnot connectivity of the current platform
	# 	#ss: the solution space tree, which is actually a dict
		
	# 	#print(solution)
	# 	#clean the redundant data
	# 	for j in range(depth,self.qNum):
	# 	    solution[j] = -1

	# 	#if we have found the solution
	# 	if depth >= self.qNum:
	# 		#construct the map
	# 		maps = {}
	# 		for i in range(0,self.qNum):
	# 			maps[i] = int(solution[i])
	# 		if self.__checkConstraint(maps):
	# 			return maps
	# 		else:
	# 			return None
	# 	else:
	# 		for i in range(0,len(ss[depth])):
	# 			if ss[depth][i] in solution[0:depth]:
	# 				continue
	# 			solution[depth] = ss[depth][i]
	# 		    #print(solution)
	# 		    #pruning 
	# 			if self.__pruning(solution,depth):
	# 				self.__backTrace(depth + 1,solution,ss)

	# #the pruning function
	# def __pruning(self,solution,depth):
	#     #Argument:
	#     #solution: the solution in this case
	#     #depth: the current depth

	# 	for c in self.cnotSets:
	# 		if c[0] <= depth and c[1] <= depth:
	# 			tmp = [solution[c[0]],solution[c[1]]]
	# 			if tmp not in self.allowSets:
	# 				solution[depth] = -1
	# 				return False            
	# 	return True

	# def __checkConstraint(self,maps):
	#     #Argument:
	#     #cnots: the cnot connectivity appeared in this QP
	#     #acnots: the allowed cnot connectivity
	#     #maps: the adjusting map

	# 	for c in self.cnotSets:
	# 		c[0] = maps[int(c[0])]
	# 		c[1] = maps[int(c[1])]
	# 		if c not in self.allowSets:
	# 			return False
	# 	return True
	# #########################the dynamic programming method####################
	# def __computeCost2(self,q1,q2,cnots,cost,costRecord,solution,sRecord):
	# 	#argument:
	# 	#q1 & q2 : the control qubit and the target qubit
	# 	#cnots: the cnot list and the first element is a bad cnot
	# 	#cost: the cost of the current alternative
	# 	#costRecord: a list recorded all the possible cost
	# 	#solution: the swapped id of the current alternative
	# 	#sRecord: a list recorded all the possible solution
	# 	######################## dynamic #########################
	# 	partMaps = [{},{}]
	# 	cnot = [q1,q2]

	# 	interQ = self.__getInterQ(q1,q2)
	# 	#compute the part map of the two case:
	# 	#1.swap the q1 and the interQ
	# 	#2.swap the q2 and the interQ
	# 	item = copy.deepcopy(interQ)
	# 	item1 = copy.deepcopy(interQ)
	# 	item1.reverse()
	# 	interQList = [item1,item]
	# 	cost += 34 * len(interQ)

	# 	#pruning function to reduce the iteration times
	# 	if costRecord != [] and cost >= min(costRecord):
	# 		return False
		
	# 	#set the value of the partMaps
	# 	for i in range(0,2):
	# 		fixed = 0 if i == 1 else 1
	# 		partMap = self.__getPartMap(cnot[fixed],cnot[i],interQList[i])	
	# 		partMaps[i] = copy.deepcopy(partMap)		

	# 	#change the qasm code and call the function if necessary
	# 	cnotSets = [copy.deepcopy(cnots),copy.deepcopy(cnots)]

	# 	for k in range(1,-1,-1):
	# 		tmpcost = cost
	# 		tmpSolution = copy.deepcopy(solution)
	# 		if k == 0:
	# 			target = q1
	# 			tmpcost += 4
	# 		else:
	# 			target = q2
	# 		tmpSolution.append(target)
	# 		badLines = -1
	# 		badq1 = -1
	# 		badq2 = -1
	# 		for j in range(0,len(cnotSets[k])):
	# 			q = cnotSets[k][j]
	# 			#cnot statement
	# 			qid = [-1,-1]
	# 			#print(q)
	# 			for t in range(0,2):
	# 				try:
	# 					qid[t] = str(partMaps[k][q[t]])
	# 				except KeyError:
	# 					qid[t] = str(q[t])	
	# 			#print(qid)					
	# 			cnotSets[k][j] = [int(qid[0]),int(qid[1])]
	# 			if badLines == -1 and cnotSets[k][j] not in self.allowSets:
	# 				#bad cnot
	# 				badLines = j
	# 				badq1 = int(qid[0])
	# 				badq2 = int(qid[1])
			
	# 		if badLines == -1:
	# 			costRecord.append(tmpcost)
	# 			sRecord.append(tmpSolution)
	# 			#there is no bad cnot 
	# 		else:
	# 			self.__computeCost2(badq1,badq2,cnotSets[k][badLines:len(cnotSets[k])],tmpcost,costRecord,tmpSolution,sRecord)
	# 	######################## dynamic #########################
	# #use swap gates to make bad cnot feasible 
	# def __adjustPartMap2(self,interQList = None):
	#     #Argument:
	# 	#interQList: the inter qubit list

	# 	newQasm = ""
	# 	qasmList = self.qasm.split(";\n")
	# 	lines = len(qasmList)
	# 	badCnot = self.badCnot
	# 	badIndex = self.cnotTimes.index(badCnot[0])

	# 	costs = []
	# 	solutions = []
	# 	self.__computeCost2(badCnot[0][0],badCnot[0][1],self.cnotTimes[badIndex:len(self.cnotTimes)],0,costs,[],solutions)
	# 	minCostIndex = costs.index(min(costs))
	# 	solution = solutions[minCostIndex]
	# 	s_i = -1
	# 	for ind in range(0,lines):
	# 		q = qasmList[ind]
	# 		m = re.match(r'cx q\[(\d+)\],q\[(\d+)\]',q)
	# 		if m != None:
	# 			#cnot
	# 			cq = int(m.group(1))
	# 			tq = int(m.group(2))
	# 			if [cq,tq] in badCnot or [tq,cq] in badCnot:
	# 				#the current cnot is the bad cnot
	# 				interQ = self.__getInterQ(cq,tq)
	# 				s_i += 1

	# 			else:
	# 				#the cnot statement is allowed in current layout
	# 				newQasm += q
	# 				newQasm += ";\n"
	# 				continue

	# 			swapQID = solution[s_i]

	# 			#################make the impossible cnot possible#############
	# 			fixedQID = cq if swapQID == tq else tq

	# 			#get the part map
	# 			partMap = self.__getPartMap(fixedQID,swapQID,interQ)

	# 			bak_q = swapQID
	# 			if fixedQID != cq:
	# 				interQ.reverse()
				
	# 			for iq in interQ:
	# 	            #the code depends on one hypothesis that there is no isolated qubit
	# 				newQasm += self.__getSwapCode(iq,swapQID)
	# 				swapQID = iq
	# 			#add a new cnot statement 
	# 			if fixedQID != cq:
	# 				newQasm += "u2(0.0,3.14159265359) q[" + str(fixedQID) + "];\n"
	# 				newQasm += "u2(0.0,3.14159265359) q[" + str(interQ[-1]) + "];\n"
	# 			newQasm += "cx q[" + str(fixedQID) + "],q[" + str(interQ[-1]) + "];\n"
	# 			if fixedQID != cq:
	# 				newQasm += "u2(0.0,3.14159265359) q[" + str(fixedQID) + "];\n"
	# 				newQasm += "u2(0.0,3.14159265359) q[" + str(interQ[-1]) + "];\n"					
				
	# 			#change the following qasm code and get the bad cnot sets
	# 			res = self.__changeID2(qasmList[ind:lines],partMap)
	# 			badCnot = res[0]
	# 			qasmList[ind:lines] = res[1]

	# 	        ############act the SWAP gate again to restore the state of the qubit######### 
	# 	        # for iq in interQ:
	# 	        #     newQasm += getSwapCode(iq,bak_q,cmap)
	# 	        #     bak_q = iq
	# 	        ############################end of this method################################
	# 		elif q == "":
	# 			pass

	# 		else:   
	# 			newQasm += q
	# 			newQasm += ";\n"     

	# 	#choose the appropriate direction of the swap gate
	# 	#so that we can use less H gates to adjust the direction of the elemental cnot gates
	# 	return newQasm
	# #########################the dynamic programming method####################


	# adjust the id to eliminate the bad cnots by back tracing
	# def __adjustBT(self):
	#     #the solution of adjusting
	# 	solution = [-1] * self.qNum

	#     #construct the solution space
	# 	solution_space = {}

	#     #record the root node
	# 	rootNode = []
	# 	for index in range(0,self.qNum):
	# 		rootNode.append(index)
	# 	solution_space[0] = rootNode

	#     #record the node whose depth is no less than 1
	# 	for ids in range(0,self.qNum):
	# 		tmp = []    
	# 		for cnot in self.allowSets:
	# 			cq = int(cnot[0])
	# 			tq = int(cnot[1])
	# 			if ids == cq and tq not in tmp:
	#                 #record the id of the qubit
	# 				tmp.append(tq)
	# 			elif ids == tq and cq not in tmp:
	# 				tmp.append(cq)
	# 			else:
	# 				pass
	# 		solution_space[ids+1] = tmp
	# 	#the solution_space is a 2-dimens list, and ss[0] means the qubit's id
	# 	#which can connect to '0' qubit.
	# 	idMap = self.__backTrace(0,solution,solution_space)
	# 	return idMap