import numpy as np
from scipy.sparse import coo_matrix

def parallel_sum(*impedences):
    """Function to compute the graph of impedences resulting from // graphs
    works for any number of impedence graphs

    Returns
    -------
    scipy.sparse.coo_matrix
        resulting impedence
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
    """Function to compute the impedences that go serial

    Returns
    -------
    scipy.sparse.coo_matrix
        resulting impedence

    """
    return sum(impedences).tocoo()


def cast_complex_system_in_real_system(sys,b):
    """Function to cast an n dimensional complex system into an
    equivalent 2n dimension real system
    the solution of the initial system is the concatenation of the real and
    imaginary part of the solution of this system:
    sol_comp = sol_real[:n]+1.0j*sol_real[n:]

    Parameters
    ----------
    sys : scipy.sparse.coo_matric
        system with complex data
    b : np.array
        second member with real or complex values

    Returns
    -------
    sys_comp
        real system equivalent to complex system
    new_b
        real second member equivalent to complex system
    """
    coords = np.stack((sys.row,sys.col),axis=1)
    data = np.array(sys.data,dtype=complex)
    b= b.astype(complex)
    new_coords = np.concatenate((coords,coords+[[0],[sys.shape[0]]],coords+[[sys.shape[0]],[0]],coords+[[sys.shape[0]],[sys.shape[0]]]),axis=1)
    new_data = np.concatenate((data.real,-data.imag,data.imag,data.real),axis=0)
    sys_comp = coo_matrix((new_data,(new_coords[0],new_coords[1])))
    new_b = np.concatenate((b.real,b.imag),axis=0)
    return sys_comp,new_b


def constant_block_diag(A,repetitions):
    """Function to repeat the matrix A multiple times along the diagonal
    This is a faster version of block_diag in the case of having always the same block

    Parameters
    ----------
    A : scpiy.sparse.coo_matrix
        block to repeat multiple times on the diagonal
    repetitions : int
        number of repetitions to perform

    Returns
    -------
    coo_matrix
        block diagonal sparse matrix
    """
    size = A.shape[0]
    indexes = A.data.shape[0]
    rows = np.tile(A.row,(repetitions,))+size*np.repeat(np.arange(0,repetitions,dtype=int),indexes)
    cols = np.tile(A.col,(repetitions,))+size*np.repeat(np.arange(0,repetitions,dtype=int),indexes)
    data = np.tile(A.data,(repetitions,))
    return coo_matrix((data,(rows,cols)),shape=(repetitions*size,repetitions*size))



def build_big_temporal_system(S1,S2,dt,rhs,sol,nb_timesteps):
    """Function to build the temporal system for nb_timesteps
    The solution of this system is the concatenation of the all the timsteps
    except for the initial timestep that the user is free to concatenate with the solutions
    Tip: reshaping the solution of the systme with shape (nb_timesteps,sol.shape[0]) provides
    the solution array indexed by the timestep

    Parameters
    ----------
    S1 : coo_matrix
        real part of the temporal system
    S2 : coo_matrix
        derivative part of the temporal system
    dt : float
        timestep for the simulation
    rhs : np.array
        second member of the system
    sol : initial condition
        solution of initial condition system
    nb_timesteps : int
        number of timesteps

    Returns
    -------
    S : coo_matrix
        system left hand side for the temporal
    """
    A = constant_block_diag((S2+dt*S1).tocoo(),nb_timesteps)
    B = constant_block_diag(-S2,nb_timesteps-1)
    B = coo_matrix((B.data,(B.row+S1.shape[0],B.col)),shape=A.shape)
    S = A+B
    RHS = np.concatenate([rhs*dt+S2@sol]+[rhs*dt]*(nb_timesteps-1),axis=0)
    return S,RHS


