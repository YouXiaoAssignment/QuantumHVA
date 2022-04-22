import numpy as np
from qiskit import *
import pybobyqa

pi = np.pi
runtime_count = 0

def parametric_state(par):
    """
    one layer construct
    :param v1: varies of optimization of \hat{x}
    :param v2: varies of optimization of \hat{x}\hat{x}
    :param v3: ...
    :param v4: ...
    :param v5: ...
    :param v6: ...
    :return: Ansatz-target state
    """
    q = QuantumRegister(qubit_number)
    c = ClassicalRegister(qubit_number)
    circuit = QuantumCircuit(q, c)
    number = int((par.size)/6)
    #circuit.x(1)#初始条件的设定，以下6行可以去掉
    #circuit.x(3)
    #circuit.x(5)
    #circuit.x(7)
    #circuit.x(9)
    #circuit.x(11)
    for j in range(number):
        """
        The hamiltonian for X part
        """
        circuit.rx(2 * par[6 * j], range(qubit_number))
        for px in range(nx):
            circuit.h(xp[2*px])
            circuit.h(xp[2*px+1])
            circuit.cnot(xp[2*px+1], xp[2*px])
            circuit.rz(2 * par[(6 * j) + 1], xp[2*px])
            circuit.cnot(xp[2 * px + 1], xp[2 * px])
            circuit.h(xp[2 * px + 1])
            circuit.h(xp[2 * px])
        """
        The hamiltonian for Y part
        """
        circuit.ry(2 * par[(6 * j) + 2], range(qubit_number))
        for py in range(ny):
            circuit.rx(pi/2, yp[2*py])
            circuit.rx(pi/2, yp[2*py+1])
            circuit.cnot(yp[2*py+1], yp[2*py])
            circuit.rz(2 * par[(6 * j) + 3], yp[2*py])
            circuit.cnot(yp[2 * py + 1], yp[2 * py])
            circuit.rx(-pi / 2, yp[2 * py + 1])
            circuit.rx(-pi / 2, yp[2 * py])
        """
        The hamiltonian for Z part
        """
        circuit.rz(2 * par[(6 * j) + 4], range(qubit_number))
        for pz in range(nz):
            circuit.cnot(zp[2*pz+1], zp[2*pz])
            circuit.rz(2 * par[(6 * j) + 5], zp[2*pz])
            circuit.cnot(zp[2*pz+1], zp[2*pz])
    return circuit


def XX(par, i, j):
    """
    Gives expectation value of XX
    sub-hamiltonian from measurement
    on parametric state

    :return: expectation value of XX
    """
    circuit = parametric_state(par)
    ###########  Transformation on XX ###########
    q = circuit.qregs[0]
    c = circuit.cregs[0]
    # circuit.ry(-np.pi / 2, q[i])
    # circuit.ry(-np.pi / 2, q[j])
    circuit.h(q[i])
    circuit.h(q[j])
    circuit.cnot(q[i], q[j])
    circuit.measure(q, c)
    ############  XX measurement  ##############
    exp_XX = double_measurement(circuit, i, j)
    return exp_XX


def YY(par, i, j):
    """
    Gives expectation value of YY
    sub-hamiltonian from measurement
    on parametric state.
    :param theta: angle in radian
    :return: expectation value of YY
    """
    circuit = parametric_state(par)
    ########### Transformation on YY #########
    q = circuit.qregs[0]
    c = circuit.cregs[0]
    # circuit.rx(np.pi / 2, q[i])
    # circuit.rx(np.pi / 2, q[j])
    circuit.sdg(q[i])
    circuit.sdg(q[j])
    circuit.cnot(q[i], q[j])
    circuit.measure(q, c)
    ############ YY Measurement ###########
    exp_YY = double_measurement(circuit, i, j)
    return exp_YY

def ZZ(par, i, j):
    """
    Gives expectation value of ZZ
    sub-hamiltonian from measurement
    on parametric state.
    :return: expectation value of ZZ
    """
    circuit = parametric_state(par)
    ##########################################
    q = circuit.qregs[0]
    c = circuit.cregs[0]
    circuit.measure(q, c)
    ###########  ZZ measurement ###############
    exp_ZZ = double_measurement(circuit, i, j)
    return exp_ZZ

def X(par, i):
    """
    Gives expectation value of X
    sub-hamiltonian from measurement
    on parametric state.
    :return: expectation value of X
    """
    circuit = parametric_state(par)
    q = circuit.qregs[0]
    c = circuit.cregs[0]
    # circuit.ry(-np.pi / 2, q[i])
    circuit.h(q[i])
    circuit.measure(q, c)
    exp_X = single_measurement(circuit, i)
    return exp_X

def Y(par, i):
    """
    Gives expectation value of Y
    sub-hamiltonian from measurement
    on parametric state.
    :return: expectation value of Y
    """
    circuit = parametric_state(par)
    q = circuit.qregs[0]
    c = circuit.cregs[0]
    # circuit.rx(np.pi / 2, q[i])
    circuit.sdg(q[i])
    circuit.h(q[i])
    circuit.measure(q, c)
    exp_Y = single_measurement(circuit, i)
    return exp_Y


def Z(par, i):
    """
    Gives expectation value of Z
    sub-hamiltonian from measurement
    on parametric state.
    :return: expectation value of Z
    """
    circuit = parametric_state(par)
    q = circuit.qregs[0]
    c = circuit.cregs[0]
    circuit.measure(q, c)
    exp_Z = single_measurement(circuit, i)
    return exp_Z


def vqe(par):  # ------------------- creates ansatz measures and performs addition to get the eigen value
    """
    Contain the complete Hamiltonian
    :return: expectation value of whole hamiltonian
    """
    Jx = -1
    Jy = -1
    Jz = -1
    hx = 0
    hy = 0
    hz = 0
    E = 0
    Mz = 0
    Mx = 0

    plaquette=0

    for jx in range(nx):
        E = E - Jx * (XX(par, xp[2*jx], xp[2*jx+1]))
    for jy in range(ny):
        E = E - Jy * (YY(par, yp[2*jy], yp[2*jy+1]))
    for jz in range(nz):
        E = E - Jz * (ZZ(par, zp[2*jz], zp[2*jz+1]))
    for i in range(qubit_number):
        E = E + hx * X(par, i) + hy * Y(par, i) + hz * Z(par, i)
    for j in range(qubit_number):
        Mz = Mz + Z(par, j)
    for j in range(qubit_number):
        Mx = Mx + X(par, j)


    "拓扑量子数plaquette，下面仅适用于论文中的晶格构型"
    for i in range(1,qubit_number-1,2):
        loop = X(par,i-1) * Z(par,i) * Y(par,(i+1)-6*((i+1)%6 == 0)) * X(par,((i+5)%qubit_number+1)) * Z(par,(i+5)%qubit_number) * Y(par,(i+5)%qubit_number-1+6*((((i+5)%qubit_number)%6==0)))
        plaquette += loop

    plaquette = abs((int(qubit_number/2)-plaquette)/qubit_number-1)
        

    Mz = Mz / qubit_number
    Mx = Mx / qubit_number
    global runtime_count
    runtime_count += 1
    if runtime_count % 5 == 0:
        print('Optimization time:{}'.format(runtime_count))
        print('Present Circuit Parameters:{}'.format(par))
        print('Current Energy:{}'.format(E))
        print('Plaquette:{}'.format(plaquette))
        print('normalized total magnetization in z direction:{}'.format(Mz))
        print('normalized total magnetization in x direction:{}'.format(Mx))
    return E


def doukey_check(my_dict: dict, my_key: str, i, j):
    """
    If key is missing returns 0
    otherwise the corresponding value.
    :param my_dict: the dictionary
    :param my_key: key (string)
    :return: 0 or value corresponding to key
    """
    response = 0
    for everyKey in my_dict:
        if everyKey[qubit_number - 1 - i] == my_key[qubit_number - 1 - i] and everyKey[qubit_number - 1 - j] == my_key[qubit_number - 1 - j]:
            response = response + my_dict[everyKey]
    return response

def sigkey_check(my_dict:dict, my_key:str, i):
    response = 0
    for everyKey in my_dict:
        if everyKey[qubit_number - 1 - i] == my_key[qubit_number - 1 - i]:
            response = response + my_dict[everyKey]
    return response

def double_measurement(circuit, i, j):  # ------------------ Takes a quantum circuit and perform measurement
    """
    Takes the quantum circuit as
    input to perform measurement
    :param circuit: quantumm circuit
    :return: expectation value
    """
    shots = 5000
    backend = BasicAer.get_backend('qasm_simulator')
    job = execute(circuit, backend, shots=shots)
    result = job.result()
    counts = result.get_counts()
    p11 = ''
    p00 = ''
    for t in range(qubit_number):
        if t == qubit_number - 1 - i or t == qubit_number - 1 - j:
            p11 = p11 + '1'
        else:
            p11 = p11 + '0'
        p00 = p00 + '0'
    p01 = ''
    for t in range(qubit_number):
        if t == qubit_number - 1 - j:
            p01 = p01 + '1'
        else:
            p01 = p01 + '0'
    p10 = ''
    for t in range(qubit_number):
        if t == qubit_number - 1 - i:
            p10 = p10 + '1'
        else:
            p10 = p10 + '0'
    n00 = doukey_check(counts, p00, i, j)
    n01 = doukey_check(counts, p01, i, j)
    n10 = doukey_check(counts, p10, i, j)
    n11 = doukey_check(counts, p11, i, j)
    expectation_value = ((n00 + n11) - (n01 + n10)) / shots
    return expectation_value


def single_measurement(circuit, i):
    """
        Takes the quantum circuit as
        input to perform measurement
        :param circuit: quantumm circuit
        :return: expectation value
        """
    shots = 5000
    backend = BasicAer.get_backend('qasm_simulator')
    job = execute(circuit, backend, shots=shots)
    result = job.result()
    counts = result.get_counts()
    p1 = ''
    p2 = ''
    for t in range(qubit_number):
        if t == qubit_number - 1 - i:
            p1 = p1 + '1'
        else:
            p1 = p1 + '0'
        p2 = p2 + '0'
    n1 = sigkey_check(counts, p1, i)
    n0 = sigkey_check(counts, p2, i)
    expectation_value = (n0 - n1) / shots
    return expectation_value

qubit_number = 12
xp = np.array([1, 2, 3, 4, 7, 8, 9, 10])
yp = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
zp = np.array([2, 7, 4, 9])
nx = int((xp.size)/2)
ny = int((yp.size)/2)
nz = int((zp.size)/2)
par_array = np.array([1, 1, 0.1, 1, 0.1, 0.1])
result = pybobyqa.solve(vqe, par_array)
print(result)