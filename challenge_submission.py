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

Your Name : Zhang xin, Xiang hong, Li yanhong, Shen ruping

Your E-Mail :zhang.x@cqu.edu.cn

Description of the algorithm :

- How does the algorithm work?
- Did you use any previously published schemes? Cite relevant papers.
- What packages did you use in your code and for what part of the algorithm?
- How general is your approach? Does it work for arbitrary coupling layouts (qubit number)?
- Are there known situations when the algorithm fails?


The optimization scheme we designed is divided into three steps: 
the whole qubit adjustment, the local qubit adjustment and the mergence of the single qubit gates. 
In the first, we extracte all the CNOT operations in the current quantum program and sorting them according to the order they perform in the quantum program.
Then, traversing the CNOT operations fromfront to back, and conducting qubit ID adjustment to each illegal CNOT encountered, since each adjustment will change all of the qubit ID, so we should make sure that this adjustment does not make the previous CNOT operation illegal. And then, continue to adjust backwards until we have no way to modify the current illegal CNOT operation by adjusting the qubit ID without affecting the previous CNOT operation. 
Secondly, we do the following design in the algorithm: After switching the state of the qubit by using SWAP gates, it is not necessary to use these SWAP gates to restore it again, instead, we use the Qubit ID that participates in transformation as a mapping to adjust the qubit ID to the subsequent code(This is why we call it local qubit adjustment). 
In the end, we combine the single qubits which act on a same qubit in possible.
In order to write and execute programs more conveniently, we introduced the following three Python modules:
1. re
2. copy
3. sympy
It should be noted that all of these Python modules can be installed directly through PIP execution or CONDA instructions, and the code can be successfully executed on Windows and MAC systems.
Our method is general and works well for arbitrary coupling layouts and all qubits number. For more details, you can read the IBM Developer challenge code manual.

---------------> please fill out this section <---------------
"""

# Include any Python modules needed for your implementation here
import re
from Optimize import *

# The following class is the input and output circuit representation for a
# QISKit compiler
from qiskit.dagcircuit import DAGCircuit
from challenge_evaluation import *


def compiler_function(dag_circuit,coupling_map=None,gate_costs=None):
    #Instantiate an object of that class
    opt = Optimize(dag_circuit,coupling_map)
    #print(opt.qasm)
    #print(dag_circuit.qasm())

    #there are two cases:
    #1. the original qasm code can be exeucted on the platform after adjusted the id of qubit
    #2. the original qasm can't be exeucted by adjusting the id, then we have to use SWAP gates
    oriQASM = dag_circuit.qasm()
    opt.setQASM(oriQASM)

    #the QASM can be executed on the targat platform
    tarQASM = opt.qasm
    if opt.badCnot != []:
        oriQASM = opt.adjustID()
        #generate the final qasm code
    tarQASM = opt.generateQASM(oriQASM)
    tarQASM = opt.reduceSameGate(tarQASM)
    #generate the instance of DAGCircuit
    resCircuit = qasm_to_dag_circuit(tarQASM)
    return resCircuit


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
#    
#     initial_layout = None
#     coupling = Coupling(coupling_map)
#     #print(coupling)
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


if __name__ == "__main__":
    dics = {0:"circle",1:"linear",2:"neighbour",3:"center"}
    for q in range(9,10):
        for cm in range(3,4):
            for d in range(5,15):
                res = score(d,cm,q,compiler_function,'local_qasm_simulator')
                name= "d" + str(d) + ".txt"
                location = "ha/" + str(q) + "-" + dics[cm] + "/"
                files = open(location+name,"w")
                strs = ""
                for k in res[1].keys():
                    strs += str(res[1][k])
                    strs += "\n"
                strs += str(res[0])
                files.write(strs)
                files.close()