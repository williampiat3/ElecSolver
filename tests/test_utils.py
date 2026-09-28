from scipy.sparse import coo_matrix
from ElecSolver.utils import *


def test_block_diag():
    A = coo_matrix(([1,2,3,4],([0,0,1,1],[0,1,0,1])))
    constant_block_diag(A,3).todense()

def test_graph_components():
    all_coords = np.array([[0,1],[1,2],[3,4]])
    components = compute_graph_compontents(all_coords)
    assert len(components)==2
    assert np.all(components[0]==np.array([0,1,2]))
    assert np.all(components[1]==np.array([3,4]))