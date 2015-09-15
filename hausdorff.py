__author__ = 'p0054421'



def distanceBetweenCurves(C1, C2):
    D = scipy.spatial.distance.cdist(C1, C2, 'euclidean')

    #none symmetric Hausdorff distances
    H1 = np.max(np.min(D, axis=1))
    H2 = np.max(np.min(D, axis=0))

    return (H1 + H2) / 2.

def distanceMatrixOfCurves(Curves):
    numC = len(Curves)

    D = np.zeros((numC, numC))
    for i in range(0, numC-1):
        for j in range(i+1, numC):
            D[i, j] = D[j, i] = distanceBetweenCurves(Curves[i], Curves[j])

    return D
