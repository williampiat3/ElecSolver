from scipy.sparse import coo_matrix
from ElecSolver.utils import *


def test_block_diag():
    A = coo_matrix(([1,2,3,4],([0,0,1,1],[0,1,0,1])))
    constant_block_diag(A,3).todense()