import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix
from ElecSolver import TemporalSystemBuilder
from ElecSolver.utils import build_big_temporal_system
from mumps import Context



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

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_ground(0)
    elec_sys.build_system()
    elec_sys.build_second_member_intensity(10,1,0)
    S1,S2,S_i = elec_sys.S1,elec_sys.S2,elec_sys.S_init
    rhs = elec_sys.rhs
    S1 = coo_matrix(S1,shape=(10,10))
    S2 = coo_matrix(S2,shape=(10,10))


    S_i = coo_matrix(S_i)
    b = np.zeros(S_i.shape[0])
    b[rhs[1][0]]=rhs[0]

    sol = spsolve(S_i.tocsr(),b)
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
    # plt.clf()

def test_temporal2():
    ## Simple tetrahedron
    res_coords  = np.array([[0,2],[1,3]],dtype=int)
    res_data = np.array([1,1],dtype=float)

    coil_coords  = np.array([[1,0],[2,3]],dtype=int)
    coil_data = np.array([0.1,0.1],dtype=float)

    capa_coords = np.array([[1,2],[3,0]],dtype=int)
    capa_data = np.array([2,2],dtype=float)

    ## total impedance
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_ground(0)
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
    sol = spsolve(S_i.tocsr(),b)
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
    # plt.clf()


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

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_ground(0)
    elec_sys.build_system()
    elec_sys.build_second_member_intensity(10,1,0)
    S_i,b = elec_sys.get_init_system()
    sol = spsolve(S_i.tocsr(),b)
    S1,S2,rhs = elec_sys.get_system()
    dt=0.08
    nb_timesteps=100
    ## build big system for no for loop!
    S,RHS = build_big_temporal_system(S1,S2,dt,rhs,sol,nb_timesteps)
    ctx = Context()
    ctx.set_matrix(S)
    ctx.analyze()
    ctx.factor()
    x = ctx.solve(RHS)

    ## adding initial conditions
    sols = np.concatenate([sol[np.newaxis],x.reshape(S.shape[0]//sol.shape[0],sol.shape[0])],axis=0)


    currents_coil,currents_res,currents_capa,voltages,_ = elec_sys.build_intensity_and_voltage_from_vector(sols)

    # import matplotlib.pyplot as plt
    # plt.xlabel("Time")
    # plt.ylabel("Intensity")
    # plt.plot(currents_res[:,0],label="intensity res 1")
    # plt.plot(currents_res[:,1],label="intensity res 2")
    # plt.legend()
    # plt.savefig("test")
    # plt.clf()

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

    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_ground(0)
    elec_sys.build_system()
    elec_sys.build_second_member_tension(10,1,0)
    S_i,b = elec_sys.get_init_system()
    sol = spsolve(S_i.tocsr(),b)


def test_big_grid():
    import networkx as nx

    size=31
    G = nx.grid_2d_graph(size, size)
    pos = {(x,y):(y,-x) for x,y in G.nodes()}
    center = size**2//2


    adjacency = nx.adjacency_matrix(G).tocoo()
    mask = adjacency.row>adjacency.col
    impedence_coords = np.array([adjacency.row,adjacency.col],dtype=int)[:,mask]
    impedence_data = adjacency.data[mask]
    mask_coil = (np.abs(impedence_coords[0]-impedence_coords[1])==size)*((impedence_coords[0]%size!=(size-1)//2)*(impedence_coords[1]%size!=(size-1)//2))
    mask_res = ~mask_coil
    coords_coil = impedence_coords[:,mask_coil]
    data_coil = impedence_data[mask_coil]
    coords_res = impedence_coords[:,mask_res]
    data_res = impedence_data[mask_res]

    coords_capa = np.array([[],[]],dtype=int)
    data_capa = np.array([],dtype=float)

    mutual_coords = [[],[]]
    mutual_data = []




    electric_sys = TemporalSystemBuilder(coords_coil,data_coil,coords_res,data_res,coords_capa,data_capa,mutual_coords,mutual_data,mutual_coords,mutual_data)
    electric_sys.set_ground(0)

    electric_sys.build_second_member_intensity(intensity=2,input_node=0,output_node=center)
    electric_sys.build_second_member_intensity(intensity=2,input_node=size-1,output_node=center)
    electric_sys.build_second_member_intensity(intensity=2,input_node=size**2-1,output_node=center)
    electric_sys.build_second_member_intensity(intensity=2,input_node=size**2-size,output_node=center)

    ## building system
    electric_sys.build_system()
    S_i,b = electric_sys.get_init_system()

    ctx = Context()
    ## set scotch ordering instead of METIS

    ctx.set_matrix(S_i)
    ctx.analyze()
    ctx.factor(ordering = "scotch")
    sol = ctx.solve(b)


    S1,S2,rhs = electric_sys.get_system()

    dt=0.8
    ## Using Implicit Euler integration scheme
    A = (S2+dt*S1).tocoo()
    B = rhs*dt

    ctx.set_matrix(A)
    ctx.analyze()
    ctx.factor(ordering = "scotch")

    for i in range(100):
        b = B+S2@sol
        sol = ctx.solve(b)


    currents_coil,currents_res,currents_capa,voltages,_  = electric_sys.build_intensity_and_voltage_from_vector(sol)
    intensities_sparse = coo_matrix((np.concatenate([currents_coil,currents_res],axis=0),(np.concatenate([coords_coil[0],coords_res[0]],axis=0),np.concatenate([coords_coil[1],coords_res[1]],axis=0))),shape=(size**2,size**2))
    intensities_sparse = intensities_sparse - intensities_sparse.T

    graph =  nx.from_scipy_sparse_array(intensities_sparse)
    # Layout
    weights = [np.abs(graph[u][v]['weight']**(1/5)) for u, v in graph.edges()]
    # # Normalize weights for colormap

    # import matplotlib.pyplot as plt
    # import matplotlib
    # norm = matplotlib.colors.Normalize(vmin=min(weights), vmax=max(weights))
    # cmap = matplotlib.cm.viridis  # You can also try plasma, inferno, coolwarm, etc.
    # edge_colors = [cmap(norm(w)) for w in weights]

    # pos = {n:(y,-x) for (x,y),n in zip(G.nodes(),graph.nodes)}
    # # Plot the graph
    # plt.figure(figsize=(8, 6))
    # nx.draw_networkx_edges(graph, pos, edge_color=edge_colors, width=np.array(weights)*4)
    # plt.title(f"Input node 0,{size-1},{size**2-size},{size**2-1} Output node {center}")
    # plt.axis('off')
    # plt.savefig("resistance_grid.png")
    # plt.clf()

def freq_simulation():
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

    fvect = [10,20,100,1000]
    elec_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    elec_sys.set_ground(0)
    elec_sys.build_system()
    elec_sys.build_second_member_intensity(10,1,0)
    sys1,sys2,rhs_ref = elec_sys.get_system()
    fake_sys,rhs = (sys1+10j*sys2).tocoo(),rhs_ref
    ctx = Context()
    ctx.set_matrix(fake_sys)
    ctx.analyze()
    for f in fvect:
        frequency_system,_ = (sys1+1j*np.pi*2*f*sys2).tocoo(),rhs_ref
        ctx.factor(frequency_system,reuse_analysis=True)
        sol = ctx.solve(rhs_ref)

        currents_coil,currents_res,currents_capa,voltages,current_source= elec_sys.build_intensity_and_voltage_from_vector(sol)

def test_hydraulic():
    ## Defining resistances
    res_coords  = np.array([[0,2,1,0,1,3],[1,3,3,2,2,0]],dtype=int)
    res_data = np.array([2,3,1,1,1,1],dtype=float)
    ## Defining coils
    coil_coords  = np.array([[],[]],dtype=int)
    coil_data = np.array([],dtype=float)
    ## Defining capacities
    capa_coords = np.array([[],[]],dtype=int)
    capa_data = np.array([],dtype=float)

    ## Defining empty mutuals here
    mutuals_coords=np.array([[],[]],dtype=int)
    mutuals_data = np.array([],dtype=float)


    res_mutuals_coords=np.array([[],[]],dtype=int)
    res_mutuals_data = np.array([],dtype=float)

    ## initializing system
    hydraulic_sys = TemporalSystemBuilder(coil_coords,coil_data,res_coords,res_data,capa_coords,capa_data,mutuals_coords,mutuals_data,res_mutuals_coords,res_mutuals_data)
    ## Seting ground at point 0
    hydraulic_sys.set_ground(0)
    ## Build second member
    hydraulic_sys.build_system()
    ## enforcing a pressure delta of 10 Pa
    hydraulic_sys.build_second_member_tension(10,1,0)
    # get system (S1 is real part, S2 derivative part)
    # the problem is only resitive thus S2 =0
    S1,S2,rhs = hydraulic_sys.get_system()

    sol = spsolve(S1.tocsr(),rhs)
    solution = hydraulic_sys.build_intensity_and_voltage_from_vector(sol)
    # After you computed the solution of the system

    pressure_input=10000
    pressure_node=0
    # Rescaling the potential to the new reference
    potentials = solution.potentials - solution.potentials[pressure_node] + pressure_input
    print("Pressures in the system:", potentials)
    ## get the flux passing through the source
    print("Debit through the system",solution.intensities_sources[0])

if __name__ == "__main__":
    test_hydraulic()
    freq_simulation()
    test_big_grid()
    test_one_shot_temporal()
    test_tension()
    test_temporal()
    test_temporal2()
