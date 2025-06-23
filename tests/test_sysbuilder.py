import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix
from SystemBuilder import ElectricSystemBuilder


def test_sys_general_coo_build():
    jw = 5j
    ## sparse python res matrix
    impedence_coords = np.array([[0,0,1],[1,2,2]],dtype=int)
    impedence_data = np.array([1,1,1],dtype=complex)

    ## total impedance
    mutuals_coords=np.array([[0],[1]],dtype=int)
    mutuals_data = np.array([2.j],dtype=complex)


    electric_sys = ElectricSystemBuilder(impedence_coords,impedence_data,mutuals_coords,mutuals_data)
    # setting the mass
    electric_sys.set_mass(0)
    electric_sys.build_system()
    electric_sys.build_second_member_tension(tension=10,input_node=1,output_node=0)
    ## Need to evaluate the system because it was altered when calling the second member
    sys,b = electric_sys.get_system()
    print(sys.todense())
    print(b)
    sol = spsolve(sys.tocsr(),b)
    intensities,potentials = electric_sys.build_intensity_and_voltage_from_vector(sol)
    print(intensities)
    print(potentials)


def test_sys_general_mutual_intensity():
    jw = 5j
    ## sparse python res matrix
    impedence_coords = np.array([[0,0,1,3],[1,2,2,4]],dtype=int)
    impedence_data = np.array([1,1j,1,1j],dtype=complex)

    ## total impedance
    mutuals_coords=np.array([[1],[3]],dtype=int)
    mutuals_data = np.array([2.j],dtype=complex)


    electric_sys = ElectricSystemBuilder(impedence_coords,impedence_data,mutuals_coords,mutuals_data)
    # setting the mass
    electric_sys.set_mass(0,3)
    electric_sys.build_system()
    electric_sys.build_second_member_intensity(intensity=10,input_node=1,output_node=0)
    ## Need to evaluate the system because it was altered when calling the second member
    sys,b = electric_sys.get_system()
    print(sys.todense())
    print(b)
    sol = spsolve(sys.tocsr(),b)
    intensities,potentials = electric_sys.build_intensity_and_voltage_from_vector(sol)
    print(intensities)
    print(potentials)



def test_res_grid():
    import networkx as nx
    import matplotlib.pyplot as plt
    import matplotlib
    G = nx.grid_2d_graph(7, 7)
    pos = {(x,y):(y,-x) for x,y in G.nodes()}

    adjacency = nx.adjacency_matrix(G).tocoo()
    mask = adjacency.row>adjacency.col
    impedence_coords = np.array([adjacency.row,adjacency.col],dtype=int)[:,mask]
    impedence_data = adjacency.data[mask]
    mutual_coords = [[],[]]
    mutual_data = []



    electric_sys = ElectricSystemBuilder(impedence_coords=impedence_coords,impedence_data=impedence_data,mutuals_coords=mutual_coords,mutuals_data=mutual_data)
    electric_sys.set_mass(0)
    electric_sys.build_second_member_intensity(intensity=2,input_node=0,output_node=24)
    electric_sys.build_second_member_intensity(intensity=2,input_node=42,output_node=24)
    electric_sys.build_second_member_intensity(intensity=2,input_node=48,output_node=24)
    electric_sys.build_second_member_intensity(intensity=2,input_node=6,output_node=24)

    ## building system
    electric_sys.build_system()
    sys,b = electric_sys.get_system()
    sol = spsolve(sys.tocsr(),b)
    intensities,voltages = electric_sys.build_intensity_and_voltage_from_vector(sol)
    intensities_sparse = coo_matrix((intensities,(impedence_coords[0],impedence_coords[1])),shape=(49,49))
    intensities_sparse = intensities_sparse - intensities_sparse.T
    graph =  nx.from_scipy_sparse_array(intensities_sparse)
    # Layout
    weights = [np.abs(graph[u][v]['weight']**(1/5)) for u, v in graph.edges()]
    # Normalize weights for colormap
    norm = matplotlib.colors.Normalize(vmin=min(weights), vmax=max(weights))
    cmap = matplotlib.cm.viridis  # You can also try plasma, inferno, coolwarm, etc.
    edge_colors = [cmap(norm(w)) for w in weights]

    pos = {n:(y,-x) for (x,y),n in zip(G.nodes(),graph.nodes)}
    # Plot the graph
    plt.figure(figsize=(8, 6))
    nx.draw_networkx_nodes(graph, pos, node_color='lightblue', node_size=700)
    nx.draw_networkx_labels(graph, pos)
    # Draw edges with color and width
    nx.draw_networkx_edges(graph, pos, edge_color=edge_colors, width=np.array(weights)*10)

    plt.title("Input node 0,6,42,48 Output node 24")
    plt.axis('off')
    plt.savefig("resistance_grid.png")
    plt.show()




if __name__ == "__main__":
    test_res_grid()
    test_sys_general_coo_build()
    test_sys_general_mutual_intensity()