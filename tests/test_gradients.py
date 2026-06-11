import numpy as np
from scipy.sparse.linalg import spsolve
from ElecSolver import TemporalSystemBuilder, FrequencySystemBuilder


def test_backpropagation():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([1,1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[0],[1]],dtype=int)
    mutuals_data = np.array([1],dtype=float)


    res_mutuals_coords=np.array([[0],[1]],dtype=int)
    res_mutuals_data = np.array([1],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.add_current_source(10,1,0)
    elec_sys.set_ground(0)
    elec_sys.build_system()

    S_i,b = elec_sys.get_init_system(sparse_rhs=True)
    S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)

    dS1 = np.array([3**i for i in range(len(S1.data))])
    dS2 = np.array([3**i for i in range(len(S2.data))])
    dS_i = np.array([3**i for i in range(len(S_i.data))])
    drhs = np.array([3**i for i in range(len(rhs.data))])
    elec_sys.backpropagate_gradients(dS1=dS1)
    elec_sys.backpropagate_gradients(dS2=dS2)
    elec_sys.backpropagate_gradients(dS_init=dS_i)
    elec_sys.backpropagate_gradients(drhs=drhs)
    dS1 = np.ones_like(S1.data)*0.1
    dS2 = np.ones_like(S2.data)*0.001
    dS_i = np.ones_like(S_i.data)*0.01
    drhs = np.ones_like(rhs.data,dtype=float)
    elec_sys.backpropagate_gradients(dS1=dS1,dS2=dS2,dS_init=dS_i,drhs=drhs)


def test_optim_res():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([1,1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.add_current_source(10,1,0)
    elec_sys.set_ground(0)
    elec_sys.build_system()

    S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
    sol =spsolve(S_i.tocsr(),rhs.todense())
    ## Solution when res_data = np.array([4,1],dtype=float)
    sol_target = [ 0.,0,-2.,-8.,8.,-8.,0.,8.,0.,8.]
    for i in range(1000):
        ## computing gradients
        db = 2*spsolve(S_i.tocsr().T, sol - sol_target)
        drhs = db[rhs.row]
        dS_i = -( db[S_i.row]*sol[S_i.col])
        gradients = elec_sys.backpropagate_gradients(dS_init=dS_i)# ,drhs=drhs)
        elec_sys.res_data = elec_sys.res_data - 0.01*gradients.res_data
        elec_sys.res_data[1] = 1.
        elec_sys.build_system()
        S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
        sol = spsolve(S_i.tocsr(),rhs.todense())
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(elec_sys.res_data, np.array([4,1],dtype=float))

def test_optim_capa():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([1,1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.add_current_source(10,1,0)
    elec_sys.set_ground(0)
    elec_sys.build_system()

    S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
    S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
    sol_init =spsolve(S_i.tocsr(),rhs.todense())
    dt=0.8
    B = rhs*dt+S2@sol_init
    A = S2+dt*S1
    sol = spsolve(A,B)
    ## Solution when capa_data = np.array([0.1,1],dtype=float) for first timestep
    sol_target = [ 3.24786325, -1.1965812,  -6.16809117,  0.61253561,  0.58404558, -2.63532764, 0., 6.16809117,  2.10826211,  1.4957265 ]
    for i in range(1000):
        ## computing gradients
        dB = 2*spsolve(A.T, sol - sol_target)
        ## chain rule for gradients of capa_data (S2 appears twice in the computation graph)
        dS2 = -( dB[S2.row]*sol[S2.col])+(dB[S2.row]*sol_init[S2.col])
        gradients = elec_sys.backpropagate_gradients(dS2=dS2)# ,drhs=drhs)
        elec_sys.capa_data = elec_sys.capa_data - 0.01*gradients.capa_data
        elec_sys.capa_data[1] = 1.
        elec_sys.build_system()
        S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
        B = rhs*dt+S2@sol_init
        A = S2+dt*S1
        sol = spsolve(A,B)
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(elec_sys.capa_data, np.array([0.1,1],dtype=float))

def test_optim_coil():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([1,1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.add_current_source(10,1,0)
    elec_sys.set_ground(0)
    elec_sys.build_system()

    S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
    S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
    sol_init =spsolve(S_i.tocsr(),rhs.todense())
    dt=0.8
    B = rhs*dt+S2@sol_init
    A = S2+dt*S1
    sol = spsolve(A,B)
    ## Solution when coil_data = np.array([3,1],dtype=float) for first timestep
    sol_target = np.array([ 1.02628285,-2.27784731,-5.57015714,-1.1257127,3.40356001,-2.15199555,0.,5.57015714, 1.72159644, 2.84730914])
    for i in range(1000):
        ## computing gradients
        dB = 2*spsolve(A.T, sol - sol_target)
        ## chain rule for gradients of coil_data (S2 appears twice in the computation graph)
        dS2 = -( dB[S2.row]*sol[S2.col])+(dB[S2.row]*sol_init[S2.col])
        gradients = elec_sys.backpropagate_gradients(dS2=dS2)# ,drhs=drhs)
        elec_sys.coil_data = elec_sys.coil_data - 0.05*gradients.coil_data
        elec_sys.coil_data[1] = 1.
        elec_sys.build_system()
        S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
        B = rhs*dt+S2@sol_init
        A = S2+dt*S1
        sol = spsolve(A,B)
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(elec_sys.coil_data, np.array([3,1],dtype=float))

def test_optim_mutual():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([1,1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    inductive_mutuals_coords=np.array([[0],[1]],dtype=int)
    inductive_mutuals_data = np.array([0.1],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,inductive_mutuals_coords,inductive_mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.add_current_source(10,1,0)
    elec_sys.set_ground(0)
    elec_sys.build_system()

    S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
    S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
    sol_init =spsolve(S_i.tocsr(),rhs.todense())
    dt=0.8
    B = rhs*dt+S2@sol_init
    print(B)
    A = S2+dt*S1
    sol = spsolve(A,B)
    ## Solution when mutuals_data = np.array([0.5],dtype=float) for first timestep
    sol_target = np.array([ 3.07692308, -3.07692308, -4.14529915,  0.2991453,   2.77777778, -2.77777778, 0., 4.14529915,  2.22222222,  1.92307692])
    for i in range(1000):
        ## computing gradients
        dB = 2*spsolve(A.T, sol - sol_target)
        ## chain rule for gradients of coil_data (S2 appears twice in the computation graph)
        dS2 = -( dB[S2.row]*sol[S2.col])+(dB[S2.row]*sol_init[S2.col])
        gradients = elec_sys.backpropagate_gradients(dS2=dS2)# ,drhs=drhs)
        elec_sys.inductive_mutual_data = elec_sys.inductive_mutual_data - 0.001*gradients.inductive_mutual_data
        elec_sys.build_system()
        S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
        B = rhs*dt+S2@sol_init
        A = S2+dt*S1
        sol = spsolve(A,B)
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(elec_sys.inductive_mutual_data, np.array([0.5],dtype=float))

def test_optim_res_mutual():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([1,1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    inductive_mutuals_coords=np.array([[],[]],dtype=int)
    inductive_mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[0],[1]],dtype=int)
    res_mutuals_data = np.array([0.1],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,inductive_mutuals_coords,inductive_mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.add_current_source(10,1,0)
    elec_sys.set_ground(0)
    elec_sys.build_system()

    S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
    S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
    sol_init =spsolve(S_i.tocsr(),rhs.todense())
    dt=0.8
    B = rhs*dt+S2@sol_init
    A = S2+dt*S1
    sol = spsolve(A,B)
    ## Solution when res mutuals_data = np.array([0.5],dtype=float) for first timestep
    sol_target = np.array([ 2.85714286, -2.85714286, -4.36507937,  0.07936508,  2.77777778, -2.77777778, 0., 4.36507937,  2.22222222,  2.14285714])
    for i in range(2000):
        ## computing gradients
        dB = 2*spsolve(A.T, sol - sol_target)
        ## chain rule for gradients of coil_data (S1 appears once in the computation graph)
        dS1 = -( dB[S1.row]*sol[S1.col])*dt
        gradients = elec_sys.backpropagate_gradients(dS1=dS1)# ,drhs=drhs)
        elec_sys.res_mutual_data = elec_sys.res_mutual_data - 0.05*gradients.res_mutual_data
        elec_sys.build_system()
        S1,S2,rhs = elec_sys.get_system(sparse_rhs=True)
        B = rhs*dt+S2@sol_init
        A = S2+dt*S1
        sol = spsolve(A,B)
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(elec_sys.res_mutual_data, np.array([0.5],dtype=float))

def test_optim_current_sources():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([1,1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.add_current_source(10,1,0)
    elec_sys.set_ground(0)
    elec_sys.build_system()

    S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
    sol =spsolve(S_i.tocsr(),rhs.todense())
    ## Solution when current_source_data = np.array([8],dtype=float)
    sol_target = np.array([ 0.,  0., -4., -4.,  4., -4.,  0.,  4.,  0.,  4.])
    for i in range(1000):
        ## computing gradients (source information is in the rhs of the system)
        db = 2*spsolve(S_i.tocsr().T, sol - sol_target)
        drhs = db[rhs.row]
        ## backpropagation of gradients to current_source_data
        gradients = elec_sys.backpropagate_gradients(drhs=drhs)
        elec_sys.current_source_data = elec_sys.current_source_data - 0.01*gradients.current_source_data
        elec_sys.build_system()
        S_i,rhs = elec_sys.get_init_system(sparse_rhs=True)
        sol = spsolve(S_i.tocsr(),rhs.todense())
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(elec_sys.current_source_data, np.array([8],dtype=float))

def test_backpropagation_frequency():
    ## Simple circuit with one coil and one res
    impedence_coords = np.array([[0,0],[1,2]],dtype=int)
    impedence_data = np.array([1,1j],dtype=complex)

    mutuals_coords=np.array([[0],[1]],dtype=int)
    mutuals_data = np.array([0.1j],dtype=complex)

    electric_sys = FrequencySystemBuilder(impedence_coords,impedence_data,mutuals_coords,mutuals_data)
    electric_sys.add_current_source(intensity=10,input_node=1,output_node=0)
    electric_sys.set_ground(0)
    electric_sys.build_system()

    S,b = electric_sys.get_system(sparse_rhs=True)
    dS = np.array([3**i+1j for i in range(len(S.data))])
    drhs = np.array([3**i+1j for i in range(len(b.data))])
    electric_sys.backpropagate_gradients(dS=dS,drhs=drhs)

def test_optim_mutual_complex():
    ## sparse python res matrix
    impedence_coords = np.array([[0,0,1],[1,2,2]],dtype=int)
    impedence_data = np.array([1,1,1],dtype=complex)

    ## mutuals
    mutuals_coords=np.array([[0],[1]],dtype=int)
    mutuals_data = np.array([2.j],dtype=complex)


    electric_sys = FrequencySystemBuilder(impedence_coords,impedence_data,mutuals_coords,mutuals_data)
    electric_sys.add_voltage_source(voltage=10,input_node=1,output_node=0)
    # setting the ground
    electric_sys.set_ground(0)
    electric_sys.build_system()

    ## Need to evaluate the system because it was altered when calling the second member
    sys,b = electric_sys.get_system(sparse_rhs=True)
    sol = spsolve(sys.tocsr(),b.todense())
    ## Target solution when impedence_data = np.array([1+1j,1,1],dtype=complex)
    sol_target = np.array([-2.+4.j, -1.+2.j,  1.-2.j,  3.-6.j,  0.+0.j, 10.+0.j,  9.+2.j])
    for i in range(3000):
        ## computing gradients
        db = 2*spsolve(sys.tocsr().conj().T, sol - sol_target)
        dS = -( db[sys.row]*np.conj(sol[sys.col]))
        gradients = electric_sys.backpropagate_gradients(dS=dS)# ,drhs=drhs)
        electric_sys.impedence_data = electric_sys.impedence_data - 0.01*gradients.impedence_data
        electric_sys.build_system()
        sys,b = electric_sys.get_system(sparse_rhs=True)
        sol = spsolve(sys.tocsr(),b.todense())
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(electric_sys.impedence_data, np.array([1+1j,1,1],dtype=complex))

def test_optim_source_complex():
    ## sparse python res matrix
    impedence_coords = np.array([[0,0,1],[1,2,2]],dtype=int)
    impedence_data = np.array([1,1,1],dtype=complex)

    ## mutuals
    mutuals_coords=np.array([[0],[1]],dtype=int)
    mutuals_data = np.array([2.j],dtype=complex)


    electric_sys = FrequencySystemBuilder(impedence_coords,impedence_data,mutuals_coords,mutuals_data)
    electric_sys.add_voltage_source(voltage=10,input_node=1,output_node=0)
    # setting the ground
    electric_sys.set_ground(0)
    electric_sys.build_system()

    ## Need to evaluate the system because it was altered when calling the second member
    sys,b = electric_sys.get_system(sparse_rhs=True)
    sol = spsolve(sys.tocsr(),b.todense())
    ## Target solution when voltage_source_data = np.array([5],dtype=complex)
    sol_target = np.array([-1.66666667+1.66666667j,-0.83333333+1.66666667j,0.83333333-1.66666667j,2.5       -3.33333333j,  0.        +0.j,          5.        +0.j, 4.16666667+1.66666667j])
    for i in range(3000):
        ## computing gradients
        db = 2*spsolve(sys.tocsr().conj().T, sol - sol_target)
        drhs = db[b.row]
        gradients = electric_sys.backpropagate_gradients(drhs=drhs)
        electric_sys.voltage_source_data = electric_sys.voltage_source_data - 0.01*gradients.voltage_source_data
        electric_sys.build_system()
        sys,b = electric_sys.get_system(sparse_rhs=True)
        sol = spsolve(sys.tocsr(),b.todense())
    ## Checking whether we converged to the right solution
    np.testing.assert_allclose(electric_sys.voltage_source_data, np.array([5],dtype=complex))


if __name__ == "__main__":
    # test_backpropagation()
    # test_optim_res()
    # test_optim_capa()
    # test_optim_coil()
    # test_optim_mutual()
    # test_optim_res_mutual()
    # test_optim_current_sources()
    # test_backpropagation_frequency()
    # test_optim_mutual_complex()
    test_optim_source_complex()