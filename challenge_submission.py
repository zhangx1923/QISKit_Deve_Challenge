# -*- coding: utf-8 -*-

#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =============================================================================

"""
---------------> please fill out this section <---------------

Your Name : Zhangxin

Your E-Mail :zhang.x@cqu.edu.cn

Description of the algorithm :

- How does the algorithm work?
- Did you use any previously published schemes? Cite relevant papers.
- What packages did you use in your code and for what part of the algorithm?
- How general is your approach? Does it work for arbitrary coupling layouts (qubit number)?
- Are there known situations when the algorithm fails?


---------------> please fill out this section <---------------
"""

# Include any Python modules needed for your implementation here
import re

# The following class is the input and output circuit representation for a
# QISKit compiler
from qiskit.dagcircuit import DAGCircuit
from challenge_evaluation import *


def compiler_function(dag_circuit,coupling_map=None,gate_costs=None):
    #the original DAGCircuit
    oriQASM = dag_circuit.qasm()
    print(oriQASM)
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!")
    cnotSets = getCNOT(oriQASM)

    #the QASM can be executed on the targat platform
    tarQASM = ""

    #judge how to compiled the original qasm so that it can satisfies the requirement of the target platform
    notAllowCNOT = getIll(cnotSets , dicToList(coupling_map))

    if notAllowCNOT == []:
        #the original QASM don't need to compiled
        tarQASM = oriQASM       
    else:
        allowSets = getAllow(coupling_map)
        needHandle = getIll(cnotSets , allowSets)

        if needHandle == []:
            #we only need to convert the direction of the cnot in notAllowCNOT
            tarQASM = adjustDire(oriQASM , coupling_map)
        else:
            #some other skills
            #there are two situations:
            #1. the original qasm code can be exeucted on the platform after adjusted the id of qubit
            #2. the original qasm can't be exeucted by adjusting the id
            badCnot = copy.deepcopy(needHandle)

            #judge whether can be executed
            idMap = adjustID(cnotSets , allowSets , dag_circuit)

            if idMap == None:
                #can't be changed
                #make the impossible cnot gate possible with the help of SWAP
                tarQASM = useSwap(oriQASM , badCnot , allowSets)
                tarQASM = adjustDire(tarQASM , coupling_map)
            else:
                #change the qubit's id according to the map genereated by adjustID
                tarQASM = changeID(oriQASM , idMap)

    #generate the instance of DAGCircuit
    resCircuit = qasm_to_dag_circuit(tarQASM)
    print("更新后的qasm：")
    print(resCircuit.qasm())
    return resCircuit

#convert a dic to list
def dicToList(d):
    l = []
    for k in d.keys():
        for i in range(0,len(d[k])):
            tmp = [int(k),int(d[k][i])]
            if tmp not in l:
                l.append(tmp)
    return l

#get the cnot connectivity of this QP
def getCNOT(qasm):
    cnotSet = []
    for q in qasm.split(";\n"):
        m = re.match(r'cx q\[(\d)\],q\[(\d)\]',q)
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
    return cnotSet

#get the allowed cnot
def getAllow(coupling_map):
    allowCnot = []
    for k in coupling_map.keys():
        for i in range(0,len(coupling_map[k])):
            tmp = [int(k),int(coupling_map[k][i])]
            tmpReser = [int(coupling_map[k][i]),int(k)]
            if tmp not in allowCnot:
                allowCnot.append(tmp)
            if tmpReser not in allowCnot:
                allowCnot.append(tmpReser)
    return allowCnot

#get the illegal cnot
def getIll(cnotSets,allowSets):
    import copy
    bakSets = copy.deepcopy(cnotSets)
    for cnot in bakSets:
        if cnot in allowSets:
            cnotSets.remove(cnot)
    return cnotSets

def judgeReverse(qid1,qid2,cmap):
    #Argument:
    #qid1 & qid2: judge whether "cx q[qid1],q[qid2]" can be executed by add H to reverse the direction
    #cmap: the coupling map of the layout
    if (qid1 not in cmap.keys()) or (qid1 in cmap.keys() and qid2 not in cmap[qid1]): 
        if qid2 in cmap.keys() and qid1 in cmap[qid2]:
            return True
    return False

#adjust the direction of the cnot in this QP
def adjustDire(qasm,cmap):
    #Argument:
    #qasm: the code of the current QP
    #cmap: the coupling map of this layout

    result = ""
    cnotStr = []

    for q in qasm.split(";\n"):
        m = re.match(r'cx q\[(\d)\],q\[(\d)\]',q)
        if m != None and judgeReverse(int(m.group(1)),int(m.group(2)),cmap):
            #add H gate to reverse the direction of the cnot gate
            tmpStr = "u3(1.5707963268,0.0,3.14159265359) q[" + str(m.group(2)) + "];\n"
            tmpStr += "u3(1.5707963268,0.0,3.14159265359) q[" + str(m.group(1)) + "];\n"
            tmpStr += "cx q[" + str(m.group(2)) + "],q[" + str(m.group(1)) + "];\n"
            tmpStr += "u3(1.5707963268,0.0,3.14159265359) q[" + str(m.group(2)) + "];\n"
            tmpStr += "u3(1.5707963268,0.0,3.14159265359) q[" + str(m.group(1)) + "];\n"                
            result += tmpStr
        elif q == "":
            pass
        else:
            result += q
            result += ";\n"

    return result

#adjust the id of qubits so that it can satisfies the requirement of the device
def adjustID(cnots,allowcnots,dagCircuit):
    #the function "DAGCircuit.width()" will return the number of the qubit used in the circuit
    qNum = dagCircuit.width()

    #the solution of adjusting
    solution = [-1] * qNum

    #construct the solution space
    solution_space = {}
    #print(allowcnots)

    #record the root node
    rootNode = []
    for index in range(0,qNum):
        rootNode.append(index)
    solution_space[0] = rootNode

    #record the node whose depth is no less than 1
    for ids in range(0,qNum):
        tmp = []    
        for cnot in allowcnots:
            cq = int(cnot[0])
            tq = int(cnot[1])
            if ids == cq and tq not in tmp:
                #record the id of the qubit
                tmp.append(tq)
            elif ids == tq and cq not in tmp:
                tmp.append(cq)
            else:
                pass
        solution_space[ids+1] = tmp
    #the solution_space is a 2-dimens list, and ss[0] means the qubit's id
    #which can connect to '0' qubit.

    a = backTrace(0,qNum,solution,cnots,solution_space,allowcnots)
    return a

#choose the feasible solution by backtracing
def backTrace(depth,N,solution,cnots,ss,allowcnots):
    #Argument:
    #depth: the current depth of the solution space tree
    #N: the max depth of the solution space
    #solution: the current solution, which is a dict
    #cnots: the cnot connectivity in the QP
    #allowcnots: the cnot connectivity of the current platform
    #ss: the solution space tree, which is actually a dict

    #clean the redundant data
    for j in range(depth,N):
        solution[j] = -1

    #if we have found the solution
    if depth >= N:
        #construct the map
        maps = {}
        for i in range(0,N):
            maps[i] = int(solution[i])
        if checkConstraint(cnots,allowcnots,maps):
            return maps
        else:
            return None
    else:
        for i in range(0,len(ss[depth])):
            if ss[depth][i] in solution[0:depth]:
                continue
            solution[depth] = ss[depth][i]
            #print(solution)
            #pruning 
            if pruning(solution,depth,cnots,allowcnots):
                backTrace(depth + 1,N,solution,cnots,ss,allowcnots)

#the pruning function
def pruning(solution,depth,cnots,acnots):
    #Argument:
    #solution: the solution in this case
    #depth: the current depth
    #cnots: the cnots connectivity appeared in the QP
    #acnots: the allowed cnots connectivity in this platform

    for c in cnots:
        if c[0] <= depth and c[1] <= depth:
            tmp = [solution[c[0]],solution[c[1]]]
            if tmp not in acnots:
                solution[depth] = -1
                return False            
    return True
def checkConstraint(cnots,acnots,maps):
    #Argument:
    #cnots: the cnot connectivity appeared in this QP
    #acnots: the allowed cnot connectivity
    #maps: the adjusting map

    for c in cnots:
        c[0] = maps[int(c[0])]
        c[1] = maps[int(c[1])]
        if c not in acnots:
            return False
    return True

#change the qubit's id according to the map
def changeID(qasm,maps):
    #Argument:
    #qasm: QASM code
    #maps: a dict, and the key is the original id, the value is new id.
    newQasm = ""
    for q in qasm.split(";\n"):
        qids = re.findall(r'\[(\d)\]',q)
        elementStr = re.split(r'\[\d\]',q)
        tmpQASM = elementStr[0]
        for t in range(0,len(qids)):
            try:
                qid = str(maps[int(qids[t])])
            except KeyError as ke:
                qid = t
            tmpQASM += "["
            tmpQASM += qid
            tmpQASM += "]"
            tmpQASM += elementStr[t+1]
        tmpQASM += ";\n"
        newQasm += tmpQASM
    return newQasm

#use swap gates to make bad cnot feasible 
def useSwap(qasm,cnots,acnots):
    #Argument:
    #qasm : QASM code
    #cnots: a list which records the bad cnot(can't be supplied by the target platform no matter how to change the qubit's id)
    #acnots: the allowed cnot connectivity of the layout
    newQasm = ""
    #print("进入了useSwap函数中")
    for q in qasm.split(";\n"):
        m = re.match(r'cx q\[(\d)\],q\[(\d)\]',q)
        if m != None:
            #cnot
            cq = int(m.group(1))
            tq = int(m.group(2))
            if [cq,tq] in cnots or [tq,cq] in cnots:
                #print(q)
                #the current cnot is the bad cnot
                #choose a qubit(interQ) which is adjacent to cq and tq
                interQs_cq = []
                interQs_tq = []
                for ac in acnots:
                    if cq in ac:
                        if ac[0] == cq:
                            interQs_cq.append(ac[1])
                        else:
                            interQs_cq.append(ac[0])
                    if tq in ac:
                        if ac[0] == tq:
                            interQs_tq.append(ac[1])
                        else:
                            interQs_tq.append(ac[0])

                interQ = -1
                #compute the intersection of the two list
                for iq in interQs_cq:
                    if iq in interQs_tq:
                        interQ = iq
                        break

#需要补充如果两个量子比特之间没有中间节点的情况
                if interQ == -1:
                    #there is no intersection of the two qubit
                    #we should introduce another interQ
                    pass

                #the code depends on one hypothesis that there is no isolated qubit
                newQasm += getSwapCode(interQ,tq)

                #add a new cnot statement 
                newcnot = "cx q[" + str(cq) + "],q[" + str(interQ) + "];\n"
                newQasm += newcnot

                #restore the state of the qubits
                newQasm += getSwapCode(interQ,tq)

        elif q == "":
            pass
        else:   
            newQasm += q
            newQasm += ";\n"       
    #choose the appropriate direction of the swap gate
    #so that we can use less H gates to adjust the direction of the elemental cnot gates
    return newQasm

#get the code of the swap gate     
def getSwapCode(qid1,qid2):
    #Argument:
    #qid1 & qid2: the two qubits need to swap

    swap = "cx q[" + str(qid1) + "],q[" + str(qid2) + "];\n"
    #cnot qid2,qid1
    swap += "u3(1.5707963268,0.0,3.14159265359) q[" + str(qid1) + "];\n"
    swap += "u3(1.5707963268,0.0,3.14159265359) q[" + str(qid2) + "];\n"
    swap += "cx q[" + str(qid1) + "],q[" + str(qid2) + "];\n"
    swap += "u3(1.5707963268,0.0,3.14159265359) q[" + str(qid1) + "];\n"
    swap += "u3(1.5707963268,0.0,3.14159265359) q[" + str(qid2) + "];\n"
    swap += "cx q[" + str(qid1) + "],q[" + str(qid2) + "];\n"
    return swap

# def compiler_function(dag_circuit, coupling_map=None, gate_costs=None):
#     """
#     Modify a DAGCircuit based on a gate cost function.
    
#     Instructions:
#         Your submission involves filling in the implementation
#         of this function. The function takes as input a DAGCircuit
#         object, which can be generated from a QASM file by using the
#         function 'qasm_to_dag_circuit' from the included 
#         'submission_evaluation.py' module. For more information
#         on the DAGCircuit object see the or QISKit documentation
#         (eg. 'help(DAGCircuit)').

#     Args:
#         dag_circuit (DAGCircuit): DAGCircuit object to be compiled.
#         coupling_circuit (list): Coupling map for device topology.
#                                  A coupling map of None corresponds an
#                                  all-to-all connected topology.
#         gate_costs (dict) : dictionary of gate names and costs.

#     Returns:
#         A modified DAGCircuit object that satisfies an input coupling_map
#         and has as low a gate_cost as possible.
#     """

#     #####################
#     # Put your code here
#     #####################
#     # Example using mapper passes in Qiskit
#     import copy
#     from qiskit.mapper import swap_mapper, direction_mapper, cx_cancellation, optimize_1q_gates, Coupling
#     from qiskit import qasm, unroll    
#     initial_layout = None
#     coupling = Coupling(coupling_map)
#     print(coupling)
#     compiled_dag, final_layout = swap_mapper(copy.deepcopy(dag_circuit), coupling, initial_layout, trials=40, seed=19)
#     # Expand swaps
#     basis_gates = "u1,u2,u3,cx,id"  # QE target basis
#     program_node_circuit = qasm.Qasm(data=compiled_dag.qasm()).parse()
#     unroller_circuit = unroll.Unroller(program_node_circuit,
#                                        unroll.DAGBackend(
#                                            basis_gates.split(",")))
#     compiled_dag = unroller_circuit.execute()
#     # Change cx directions
#     compiled_dag = direction_mapper(compiled_dag, coupling)
#     # Simplify cx gates
#     cx_cancellation(compiled_dag)
#     # Simplify single qubit gates
#     compiled_dag = optimize_1q_gates(compiled_dag)
#     #####################
#     # Put your code here
#     #####################
#     # Return the compiled dag circuit
#     return compiled_dag


# def basicAL():
#     QASM = "OPENQASM 2.0;include \"qelib1.inc\";qreg q[5];creg c[5];u3(0.634346715653683,0.880910487587674,0.827499019150487) q[1];"
#     c = qasm_to_dag_circuit(QASM)
#     return c


if __name__ == "__main__":
    #a = basicAL()
    #print(a.qasm())
    #help(DAGCircuit)
    print(score(compiler_function))
    # obj = open("circuits/random0_n16_d16.qasm","r")
    # qasm = obj.read()
    # DC = qasm_to_dag_circuit(qasm)
    # compiler_function(DC)