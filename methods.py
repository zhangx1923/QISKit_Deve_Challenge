# -*- coding: utf-8 -*-
import copy
import re
import sympy
import math

def mergeSingleGate(qasm:list):
	#Argument:
	#qasm: a list records the single gate qasm code which act on one qubit

	#Return:
	#statement: the merged single gate qasm code 
	qid = ''
	num = len(qasm)
	if num == 0:
		#no code
		return ""
	elif num == 1:
		#just one statement, don't need to merge
		statement = qasm[0]
		statement += "\n"
		return statement
	else:
		parameters = []
		#more than one statement, we should merge them
		countactID = -1
		for q in range(0,len(qasm)):
			if q == countactID:
				continue
			ss = qasm[q]
			if q != len(qasm)-1 and ss == qasm[q+1]:
				countactID = q+1
				continue
			m = re.match(r'u(\d)\((.+)\) q\[(\d+)\];',ss)
			qid = m.group(3)
			parameter = []
			if m.group(1) == '1':
				#u1
				parameter.append(float(m.group(2)))
			elif m.group(1) == '2':
				#u2
				plist = m.group(2).split(",")
				parameter.append(float(plist[0]))
				parameter.append(float(plist[1]))
			else:
				#u3
				plist = m.group(2).split(",")
				parameter.append(float(plist[0]))
				parameter.append(float(plist[1]))
				parameter.append(float(plist[2]))
			parameters.append(parameter)
		############################# begin merging ############################
		if parameters == []:
			return ""
		elif len(parameters) == 1:
			return qasm[0]
		else:
			pass
		#merge the two single gate
		#print(parameters)
		tmp = mergeU(parameters[0],parameters[1])
		#print(tmp)
		#merge the remanent single gate
		for i in range(2,len(parameters)):
			t = mergeU(tmp[-1],parameters[i])
			if len(t) == 1:
				tmp[-1] = t[0]
			else:
				tmp[-1] = t[0]
				tmp.append(t[1])
		############################# end merging ##############################
		statement = ""
		for ts in tmp: 
			if len(ts) == 1:
				#generate a u1 gate
				s = "u1(" + str(ts[0]) + ") q[" + qid + "];\n"
			elif len(ts) == 2:
				#generate a u2 gate
				s = "u2(" + str(ts[0]) + "," + str(ts[1]) + ") q[" + qid + "];\n"
			else:		
				#generate a u3 gate
				s = "u3(" + str(ts[0]) + "," + str(ts[1]) + "," + str(ts[2]) + ") q[" + qid + "];\n"
			statement += s
		return statement

def mergeU(para1,para2):
	#Argument:
	#para1 & para2 stand for the parameter of the two gate, and the type of them is list
	#u1 has one parameter, u2 two, u3 three	
	tmpPara = copy.deepcopy(para1)
	para1 = copy.deepcopy(para2)
	para2 = tmpPara
	first_U = len(para1)
	second_U = len(para2)
	resPara = []
	if first_U == 1 and second_U == 1:
		#u1 and u1
		tmp = []
		tmp.append(para1[0] + para2[0])
		resPara.append(tmp)
	elif first_U == 1 and second_U == 2:
		#u1 and u2
		tmp = []
		tmp.append(para1[0] + para2[0])
		tmp.append(para2[1])
		resPara.append(tmp)
	elif first_U == 1 and second_U == 3:
		#u1 and u3
		tmp = []
		tmp.append(para2[0])
		tmp.append(para1[0] + para2[1])
		tmp.append(para2[2])
		resPara.append(tmp)
	elif first_U == 2 and second_U == 1:
		#u2 and u1
		tmp = []
		tmp.append(para1[0])
		tmp.append(para2[0] + para1[1])	
		resPara.append(tmp)	
	elif first_U == 2 and second_U == 2:
		#u2 and u2
		tmp = []
		tmp.append(3.14159265359-para1[1]-para2[0])
		tmp.append(para1[0] + 1.570796326795)
		tmp.append(para2[1] + 1.570796326795)
		resPara.append(tmp)
	elif first_U == 2 and second_U == 3:
		#u2 and u3
		#change u2 to u3 and this case will be same with u3 and u3
		para1 = change2to3(para1)
	elif first_U == 3 and second_U == 1:
		#u3 and u1
		tmp = []
		tmp.append(para1[0])
		tmp.append(para1[1])
		tmp.append(para1[2] + para2[0])
		resPara.append(tmp)
	elif first_U == 3 and second_U == 2:
		#u3 and u2
		#change u2 to u3 and this case will be same with u3 and u3
		para2 = change2to3(para2)
	else:
		pass

	#merge the special case: u3 and u3
	if len(para1) == 3 and len(para2) == 3:
		#u3 and u3

		#combine the U3 gates
		tmp = mergeU3(para1,para2)
		resPara.append(tmp)

		#don't combine the U3 gates
		#resPara.append(para2)
		#resPara.append(para1)

	#return the parameter of the new single gate
	return resPara

def change2to3(para:list):
	#change u2 to u3
	res = []
	res.append(1.570796326795)
	res.append(para[0])
	res.append(para[1])
	return res

def mergeU3(para1:list,para2:list):
	#Argument:
	#para1 and para2 stands for the parameter list of the two u3 gate

	a1 = para1[0]
	b1 = para1[1]
	c1 = para1[2]
	a2 = para2[0]
	b2 = para2[1]
	c2 = para2[2]
	#merge u3(a1,b1,c1) and u3(a2,b2,c2)

	#compute a3,b3,c3
	res = yzyTOzyz(a1/2,(c1+b2)/2,a2/2)	
	a3 = 2*res[0]
	b3 = b1 + 2*res[1]
	c3 = 2*res[2] + c2
	# print("!!!!!!!!!!")
	# print(a1/2)
	# print((c1+b2)/2)
	# print(a2/2)
	# print(res)
	# print("!!!!!!!!!!!")
	# print("\n")
	#new method
	#res = [[0,(a1+a2),(c1+b2)],[0,(a1+a2),-(c1+b2)],[(c1+b2),(a1+a2),0],[-(c1+b2),(a1+a2),0]]
	#a3 = res[0][0]
	#b3 = b1 + res[0][1]
	#c3 = res[0][2] + c2
	
	#return a3,b3,c3 stands for three parameters of the combined u3
	return [a3,b3,c3]

def yzyTOzyz(a,b,c):
    eps = 1e-8
    xi = b
    theta1 = a
    theta2 = c
    solutions = []  
	
	#restore the result of the computing to reduce the time
    if sympy.Abs(sympy.cos(xi)).evalf() < eps:
        return [theta2 - theta1, xi, 0]
    elif sympy.Abs(sympy.sin(theta1 + theta2)).evalf() < eps:
        phi_minus_lambda = [
            sympy.pi / 2,
            3 * sympy.pi / 2,
            sympy.pi / 2,
            3 * sympy.pi / 2]
        stheta_1 = sympy.asin(sympy.sin(xi) * sympy.sin(-theta1 + theta2))
        #stheta_2 = sympy.asin(-sympy.sin(xi) * sympy.sin(-theta1 + theta2))
        stheta_2 = -stheta_1
        stheta_3 = sympy.pi - stheta_1
        stheta_4 = sympy.pi - stheta_2
        stheta = [stheta_1, stheta_2, stheta_3, stheta_4]
        phi_plus_lambda = list(map(lambda x:
                                   sympy.acos(sympy.cos(theta1 + theta2) *
                                              sympy.cos(xi) / sympy.cos(x)),
                                   stheta))
        sphi = [(term[0] + term[1]) / 2 for term in
                zip(phi_plus_lambda, phi_minus_lambda)]
        slam = [(term[0] - term[1]) / 2 for term in
                zip(phi_plus_lambda, phi_minus_lambda)]
        solutions = list(zip(stheta, sphi, slam))
    elif sympy.Abs(sympy.cos(theta1 + theta2)).evalf() < eps:
        phi_plus_lambda = [
            sympy.pi / 2,
            3 * sympy.pi / 2,
            sympy.pi / 2,
            3 * sympy.pi / 2]
        stheta_1 = sympy.acos(sympy.sin(xi) * sympy.cos(theta1 - theta2))
        #stheta_2 = sympy.acos(-sympy.sin(xi) * sympy.cos(theta1 - theta2))
        stheta_2 = stheta_1
        stheta_3 = -stheta_1
        stheta_4 = -stheta_2
        stheta = [stheta_1, stheta_2, stheta_3, stheta_4]
        phi_minus_lambda = list(map(lambda x:
                                    sympy.acos(sympy.sin(theta1 + theta2) *
                                               sympy.cos(xi) / sympy.sin(x)),
                                    stheta))
        sphi = [(term[0] + term[1]) / 2 for term in
                zip(phi_plus_lambda, phi_minus_lambda)]
        slam = [(term[0] - term[1]) / 2 for term in
                zip(phi_plus_lambda, phi_minus_lambda)]
        solutions = list(zip(stheta, sphi, slam))
    else:
        sinxi = sympy.sin(xi)
        cosxi = sympy.cos(xi)
        costheta12 = sympy.cos(theta1 + theta2)
        phi_plus_lambda = sympy.atan(sinxi * sympy.cos(theta1 - theta2) /
                                     (cosxi * costheta12))
        phi_minus_lambda = sympy.atan(sinxi * sympy.sin(-theta1 +
                                                                theta2) /
                                      (cosxi * sympy.sin(theta1 +
                                                                 theta2)))
        sphi = (phi_plus_lambda + phi_minus_lambda) / 2
        slam = (phi_plus_lambda - phi_minus_lambda) / 2
        cossphislam = sympy.cos(sphi + slam)
        acos = sympy.acos(cosxi * costheta12 / cossphislam)
        solutions.append((acos,sphi,slam))
        solutions.append((acos,sphi + sympy.pi / 2,slam + sympy.pi / 2))
        solutions.append((acos,sphi + sympy.pi / 2,slam - sympy.pi / 2))
        solutions.append((acos,sphi + sympy.pi,slam))
    # Select the first solution with the required accuracy
    for ss in solutions:
        if testSolution(ss[0],ss[1],ss[2],xi,theta1,theta2,eps):
            return ss

def testSolution(theta, phi, lamb, xi, theta1, theta2 , eps):
	sinxi = sympy.sin(xi)
	cosxi = sympy.cos(xi)
	sintheta = sympy.sin(theta)
	costheta = sympy.cos(theta)

	#the max value of the four deltas must be less than eps
	delta1 = sympy.Abs(sympy.cos(phi + lamb) * costheta - cosxi * sympy.cos(theta1 + theta2))
	if delta1.evalf() > eps:
		return False

	delta2 = sympy.Abs(sympy.sin(phi + lamb) * costheta - sinxi * sympy.cos(theta1 - theta2))
	if delta2.evalf() > eps:
		return False
	
	delta3 = sympy.Abs(sympy.cos(phi - lamb) * sintheta - cosxi * sympy.sin(theta1 + theta2))
	if delta3.evalf() > eps:
		return False
	
	delta4 = sympy.Abs(sympy.sin(phi - lamb) * sintheta - sinxi * sympy.sin(-theta1 + theta2))
	if delta4.evalf() > eps:
		return False
	return True

def getCNOT(qasm):
    cnotSet = []
    cnotTimes = []
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
    return cnotTimes
