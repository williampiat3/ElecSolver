import numpy as np
from scipy.sparse import coo_matrix
from collections import namedtuple

SolutionFrequency = namedtuple("ComplexSolution",["intensities","potentials","intensities_sources"])
SolutionTemporal = namedtuple("RealSolution",["intensities_coil","intensities_res","intensities_capa","potentials","intensities_sources"])

GradientsParametersFrequency= namedtuple("GradientsParametersFrequency",["impedence_data","mutual_data","voltage_source_data","current_source_data"])
GradientsParametersTemporal = namedtuple("GradientsParametersTemporal",["coil_data","res_data","capa_data","inductive_mutual_data","res_mutual_data","voltage_source_data","current_source_data"])

def parallel_sum(*impedences):
    """Combine any number of impedance graphs in parallel.

    Parameters
    ----------
    *impedences : scipy.sparse.coo_matrix
        Sparse impedance matrices to combine.

    Returns
    -------
    scipy.sparse.coo_matrix
        Resulting impedance matrix.
    """
    coords_tot = np.concatenate([impedence.coords for impedence in impedences],axis=1)
    data_tot = np.concatenate([impedence.data for impedence in impedences])
    current_indexes = np.arange(0,data_tot.shape[0],dtype=int)

    uniques,indexes,inverse,counts = np.unique(coords_tot,return_index=True,return_inverse=True,return_counts=True,axis=1)
    new_coords = uniques
    new_data = data_tot[indexes].astype(complex)

    remaining_coords = coords_tot.copy()
    remaining_data = data_tot.copy()
    reverse_counts = counts[inverse]
    reverse_counts_ref = reverse_counts.copy()

    reverse_counts[indexes]=0
    reverse_indexes = indexes[inverse]


    while np.max(reverse_counts)>1:
        mask = (reverse_counts>1)
        remaining_coords = remaining_coords[:,mask]

        remaining_data = remaining_data[mask]
        current_indexes=current_indexes[mask]

        uniques,indexes,inverse,counts = np.unique(remaining_coords,return_index=True,return_inverse=True,return_counts=True,axis=1)
        new_data[reverse_indexes[current_indexes[indexes]]] = new_data[reverse_indexes[current_indexes[indexes]]]*remaining_data[indexes]/(new_data[reverse_indexes[current_indexes[indexes]]]+ remaining_data[indexes])
        reverse_counts = counts[inverse]
        reverse_counts[indexes]=0

    indexed_data = np.zeros(impedences[0].shape[0]**2,dtype=complex)
    indexed_data[new_coords[0]*impedences[0].shape[0]+new_coords[1]]=new_data
    return coo_matrix((new_data,(new_coords[0],new_coords[1])))


def serie_sum(*impedences):
    """Combine any number of impedance matrices in series.

    Parameters
    ----------
    *impedences : scipy.sparse.coo_matrix
        Sparse impedance matrices to combine.

    Returns
    -------
    scipy.sparse.coo_matrix
        Resulting impedance matrix.
    """
    return sum(impedences).tocoo()


def cast_complex_system_in_real_system(sys,b):
    """Convert an ``n``-dimensional complex system to a real system.

    The equivalent real system has dimension ``2n``. The original complex
    solution can be reconstructed as
    ``sol_real[:n] + 1.0j * sol_real[n:]``.

    Parameters
    ----------
    sys : scipy.sparse.coo_matrix, shape (n, n)
        System matrix with complex data.
    b : numpy.ndarray, shape (n,)
        Right-hand side with real or complex values.

    Returns
    -------
    sys_comp : scipy.sparse.coo_matrix, shape (2n, 2n)
        Real system equivalent to the complex system.
    new_b : numpy.ndarray, shape (2n,)
        Real right-hand side equivalent to the complex right-hand side.
    """
    coords = np.stack((sys.row,sys.col),axis=0)
    data = np.array(sys.data,dtype=complex)
    b= b.astype(complex)
    new_coords = np.concatenate((coords,coords+[[0],[sys.shape[0]]],coords+[[sys.shape[0]],[0]],coords+[[sys.shape[0]],[sys.shape[0]]]),axis=1)
    new_data = np.concatenate((data.real,-data.imag,data.imag,data.real),axis=0)
    sys_comp = coo_matrix((new_data,(new_coords[0],new_coords[1])))
    new_b = np.concatenate((b.real,b.imag),axis=0)
    return sys_comp,new_b


def constant_block_diag(A,repetitions):
    """Repeat a matrix along the diagonal.

    This is faster than a general block-diagonal construction when every block
    is identical.

    Parameters
    ----------
    A : scipy.sparse.coo_matrix, shape (n, n)
        Block to repeat along the diagonal.
    repetitions : int
        Number of repetitions.

    Returns
    -------
    scipy.sparse.coo_matrix, shape (repetitions * n, repetitions * n)
        Block-diagonal sparse matrix.
    """
    size = A.shape[0]
    indexes = A.data.shape[0]
    rows = np.tile(A.row,(repetitions,))+size*np.repeat(np.arange(0,repetitions,dtype=int),indexes)
    cols = np.tile(A.col,(repetitions,))+size*np.repeat(np.arange(0,repetitions,dtype=int),indexes)
    data = np.tile(A.data,(repetitions,))
    return coo_matrix((data,(rows,cols)),shape=(repetitions*size,repetitions*size))

def compute_graph_components(all_coords):
    """Compute graph components from edge coordinates.

    Parameters
    ----------
    all_coords : numpy.ndarray, shape (2, n_edges)
        Node indices for each graph edge.

    Returns
    -------
    list of numpy.ndarray
        Node indices in each connected component.
    """
    number_of_nodes = np.max(all_coords)+1
    number_of_edges = all_coords.shape[1]
    # Initialize the root of each node to itself.
    roots = np.arange(0,number_of_nodes,dtype=int)
    # Make each point refer to its lowest-index neighbor.
    for i in range(number_of_edges):
        min_value = min(all_coords[0,i],all_coords[1,i])
        roots[all_coords[0,i]] = min(min_value,roots[all_coords[0,i]])
        roots[all_coords[1,i]] = min(min_value,roots[all_coords[1,i]])

    converged = False
    # Iteratively update each node to point to the root of its component.
    while not converged:
        previous_roots = roots.copy()
        roots = roots[roots]
        # Converged when no roots have changed in this iteration.
        converged = np.all(previous_roots==roots)

    return [np.arange(0,number_of_nodes,dtype=int)[roots==i] for i in np.unique(roots)]




def build_big_temporal_system(S1,S2,dt,rhs,sol,nb_timesteps):
    """Build a temporal system over multiple time steps.

    The solution concatenates all time steps except the initial one. Reshaping
    it to ``(nb_timesteps, sol.shape[0])`` produces an array indexed by time
    step.

    Parameters
    ----------
    S1 : scipy.sparse.coo_matrix, shape (n, n)
        Real part of the temporal system.
    S2 : scipy.sparse.coo_matrix, shape (n, n)
        Derivative part of the temporal system.
    dt : float
        Simulation time step.
    rhs : numpy.ndarray, shape (n,)
        Right-hand side of the system.
    sol : numpy.ndarray, shape (n,)
        Initial-condition solution.
    nb_timesteps : int
        Number of time steps.

    Returns
    -------
    S : scipy.sparse.coo_matrix, shape (nb_timesteps * n, nb_timesteps * n)
        Left-hand side of the expanded temporal system.
    RHS : numpy.ndarray, shape (nb_timesteps * n,)
        Right-hand side of the expanded temporal system.
    """
    A = constant_block_diag((S2+dt*S1).tocoo(),nb_timesteps)
    B = constant_block_diag(-S2,nb_timesteps-1)
    B = coo_matrix((B.data,(B.row+S1.shape[0],B.col)),shape=A.shape)
    S = A+B
    RHS = np.concatenate([rhs*dt+S2@sol]+[rhs*dt]*(nb_timesteps-1),axis=0)
    return S,RHS


