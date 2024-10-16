import numpy as np


def is_equiv_lattice(lat1, lat2, eps):
    """
    Determine if two lattices are equivalent.

    Parameters:
    lat1 : numpy.ndarray
        The first 3x3 lattice matrix.
    lat2 : numpy.ndarray
        The second 3x3 lattice matrix.
    eps : float
        Tolerance value for checking equivalence.

    Returns:
    bool
        True if the lattices are equivalent, False otherwise.
    """
    atol = 5e-4
    is_equiv_lattice = False
    
    lat1inv=np.linalg.inv(lat1)
    
    S = np.matmul(lat1inv, lat2)

    
    def equal(mat1, mat2, eps, atol):
        return np.allclose(mat1, mat2, atol=atol, rtol=eps)
    
    det_S = np.linalg.det(S)
   

    if equal(abs(det_S), 1.0, eps, atol) and equal(S, np.round(S), eps, atol):
    #if equal(abs(det_S), 1.0, eps, atol):
        is_equiv_lattice = True
        #print(S)
        #print(np.round(S))
        #print("I AM HERE 111")
    
    return is_equiv_lattice

# Example usage:
#lat1 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
#lat2 = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
#eps = 1e-5

#result = is_equiv_lattice(lat1, lat2, eps)
#print("Are the lattices equivalent?", result)

