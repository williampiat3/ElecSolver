import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix
from TemporalSystemBuilder import TemporalElectricSystemBuilder
from utils import build_big_temporal_system
from mumps import DMumpsContext



def test_temporal():
    ## Simple tetrahedron
    res_coords  = np.array([[0,0,0,2],[1,2,3,3]],dtype=int)
    res_data = np.array([1,(1+np.sqrt(5))/2,(1+np.sqrt(5))/2,1],dtype=float)

    coil_coords  = np.array([[1],[2]],dtype=int)
    coil_data = np.array([1],dtype=float)

    capa_coords = np.array([[1],[3]],dtype=int)
    capa_data = np.array([1],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalElectricSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_mass(0)
    elec_sys.build_system()
    elec_sys.build_second_member_intensity(10,1,0)
    S1,S2,S_i = elec_sys.S1,elec_sys.S2,elec_sys.S_init
    rhs = elec_sys.rhs
    S1 = coo_matrix(S1,shape=(10,10))
    S2 = coo_matrix(S2,shape=(10,10))


    S_i = coo_matrix(S_i)
    b = np.zeros(S_i.shape[0])
    b[rhs[1][0]]=rhs[0]

    sol = spsolve(S_i,b)
    print(elec_sys.build_intensity_and_voltage_from_vector(sol))
    dt=0.08
    vals = []
    vals_capa = []
    for i in range(500):
        currents_coil,currents_res,currents_capa,voltages,_ = elec_sys.build_intensity_and_voltage_from_vector(sol)
        # print("_______________________")
        # print(currents_coil)
        # print(currents_res)
        # print(currents_capa)
        # print(voltages)
        vals.append(currents_coil[0])
        vals_capa.append(currents_capa[0])
        sol = spsolve(S2+dt*S1,b*dt+S2@sol)
        if i ==299:
            print(elec_sys.build_intensity_and_voltage_from_vector(sol))



    # import matplotlib.pyplot as plt
    # plt.xlabel("Time")
    # plt.ylabel("Intensity")
    # plt.plot(vals,label="intensity coil")
    # plt.plot(vals_capa,label="intensity capa")
    # plt.legend()
    # plt.savefig("test")

def test_temporal2():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([0.1,0.1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([1,1],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalElectricSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_mass(0)
    elec_sys.build_system()
    elec_sys.build_second_member_intensity(10,1,0)
    S1,S2,S_i = elec_sys.S1,elec_sys.S2,elec_sys.S_init
    rhs = elec_sys.rhs
    S1 = coo_matrix(S1,shape=(10,10))
    S2 = coo_matrix(S2,shape=(10,10))

    S_i = coo_matrix(S_i)

    b = np.zeros(S_i.shape[0])
    b[rhs[1][0]]=rhs[0]
    print((S1+1j*S2).todense())
    sol = spsolve(S_i,b)
    print(elec_sys.build_intensity_and_voltage_from_vector(sol))
    dt=0.08
    vals = []
    vals_capa = []
    for i in range(50):
        currents_coil,currents_res,currents_capa,voltages,_ = elec_sys.build_intensity_and_voltage_from_vector(sol)
        vals.append(currents_res[1])
        vals_capa.append(currents_res[0])
        sol = spsolve(S2+dt*S1,b*dt+S2@sol)
    # import matplotlib.pyplot as plt
    # plt.xlabel("Time")
    # plt.ylabel("Intensity")
    # plt.plot(vals,label="intensity res 1")
    # plt.plot(vals_capa,label="intensity res 2")
    # plt.legend()
    # plt.savefig("test")


def test_one_shot_temporal():
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

    elec_sys = TemporalElectricSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_mass(0)
    elec_sys.build_system()
    elec_sys.build_second_member_intensity(10,1,0)
    S_i,b = elec_sys.get_init_system()
    sol = spsolve(S_i.tocsr(),b)
    S1,S2,rhs = elec_sys.get_system()
    dt=0.08
    nb_timesteps=100
    ## build big system for no for loop!
    S,RHS = build_big_temporal_system(S1,S2,dt,rhs,sol,nb_timesteps)
    ctx = DMumpsContext()
    import time
    if ctx.myid == 0:
        ctx.set_centralized_sparse(S)
        x = RHS.copy()
        ctx.set_rhs(x) # Modified in place
    ctx.run(job=6) # Analysis + Factorization + Solve
    ctx.destroy() # Cleanup

    ## adding initial conditions
    sols = np.concatenate([sol[np.newaxis],x.reshape(S.shape[0]//sol.shape[0],sol.shape[0])],axis=0)


    currents_coil,currents_res,currents_capa,voltages,_ = elec_sys.build_intensity_and_voltage_from_vector(sols)

    import matplotlib.pyplot as plt
    plt.xlabel("Time")
    plt.ylabel("Intensity")
    plt.plot(currents_res[:,0],label="intensity res 1")
    plt.plot(currents_res[:,1],label="intensity res 2")
    plt.legend()
    plt.savefig("test")

def test_tension():
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

    elec_sys = TemporalElectricSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_mass(0)
    elec_sys.build_system()
    elec_sys.build_second_member_tension(10,1,0)
    S_i,b = elec_sys.get_init_system()
    sol = spsolve(S_i.tocsr(),b)







if __name__ == "__main__":
    test_one_shot_temporal()
    test_tension()
    test_temporal()
    test_temporal2()
