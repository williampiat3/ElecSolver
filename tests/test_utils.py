from scipy.sparse import coo_matrix
from ElecSolver.utils import *


def test_block_diag():
    A = coo_matrix(([1,2,3,4],([0,0,1,1],[0,1,0,1])))
    constant_block_diag(A,3).todense()

def test_graph_components():
    all_coords = np.array([[0,1,3],[1,2,4]],dtype=int)
    components = compute_graph_components(all_coords)
    print(components)
    assert len(components)==2
    assert np.all(components[0]==np.array([0,1,2]))
    assert np.all(components[1]==np.array([3,4]))

if __name__ == "__main__":
    test_block_diag()
    test_graph_components()