import numpy as np

def get_HNF_diagonals(volume):
    id = 0  # Number of diagonals found
    tempDiag = []

    for i in range(1, volume + 1):  # Loop over possible first factors
        if volume % i != 0:
            continue
        quotient = volume // i
        for j in range(1, quotient + 1):  # Loop over possible second/third factors
            if quotient % j != 0:
                continue
            id += 1
            tempDiag.append([i, j, quotient // j])  # Construct the factor triplet

    diagonals = np.array(tempDiag).T  # Convert list to NumPy array and transpose
    return diagonals

def get_all_HNFs(volume):
    d = get_HNF_diagonals(volume)
    N = d.shape[1]

    # Count the total number of HNF matrices for the given determinant (volume)
    Nhnf = 0
    for i in range(N):
        Nhnf += d[1, i] * d[2, i] ** 2

    hnf = np.zeros((Nhnf, 3, 3), dtype=int)
    ihnf = 0

    for i in range(N):  # Loop over the permutations of the diagonal elements of the HFNs
        for j in range(d[1, i]):  # Loop over possible values of row 2, element 1
            for k in range(d[2, i]):  # Loop over possible values of row 3, element 1
                for l in range(d[2, i]):  # Loop over possible values of row 3, element 2
                    #hnf[ihnf, :, :] = np.array([[d[0, i], j, k], [0, d[1, i], l], [0, 0, d[2, i]]])
                     #FORTRAN from enumlib
                     #hnf(:,:,ihnf) = reshape((/ d(1,i),      j,     k,        &
                     #0, d(2,i),     l,        &
                     #0,      0, d(3,i)  /), (/3,3/))
                     # The reshape in fortan creat a lower trigonal matrix; it takes each three number and 
                     # put them in the columns, ...So to generate the same, we need to transpose the matrix

                    hnf[ihnf, :, :] = np.array([[d[0, i], j, k], [0, d[1, i], l], [0, 0, d[2, i]]]).T
                    ihnf += 1

    if ihnf != Nhnf:
        raise ValueError("HNF: not all the matrices were generated... (bug!)")

    return hnf

def get_HNF_2D_diagonals(volume):
    diagonals = []
    for i in range(1, volume + 1):  # Loop over possible first factors
        if volume % i != 0:
            continue
        quotient = volume // i
        diagonals.append([1, i, quotient])  # Construct the factor triplet

    diagonals = np.array(diagonals).T  # Convert to numpy array and transpose to match Fortran style
    return diagonals


def get_all_2D_HNFs(volume):
    d = get_HNF_2D_diagonals(volume)
    N = d.shape[1]

    # Count the total number of HNF matrices for the given determinant (volume)
    Nhnf = sum(d[2, :])

    hnf = np.zeros((Nhnf, 3, 3), dtype=int)
    ihnf = 0

    for i in range(N):  # Loop over the permutations of the diagonal elements of the HFNs
        for j in range(d[2, i]):  # Loop over possible values of row 2, element 1
            hnf[ihnf,:, :] = np.array([
                [d[2, i], j, 0],
                [0, d[1, i], 0],
                [0, 0, d[0, i]]
            ]).T
            ihnf += 1  # Count the HNFs and construct the next one

    if ihnf != Nhnf:
        raise ValueError("HNF: not all the matrices were generated...(bug!)")

    return hnf



if __name__ == "__main__":
    print("HNFs for vol=4 for 3D")
    print(get_all_HNFs(4))
    print("HNFs for vol=4 for 2D")
    print(get_all_2D_HNFs(4))
