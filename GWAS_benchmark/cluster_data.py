"""
module to hierarchically cluster SNP data
"""

import sys
from pysnptools.snpreader import Bed
import os.path
import pylab
import scipy.cluster.hierarchy as sch
import fastcluster as fc
import logging



def cluster_data(snpreader):
    """
    compute hierarchical clustering of snp data set in bed_fn
    """

    if isinstance(snpreader,str):
        snpreader = Bed(snpreader)
    G = snpreader.read().standardize().val

    # Generate distance matrix
    from sklearn.metrics.pairwise import euclidean_distances
    D = euclidean_distances(G, G)

    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = fc.linkage(D, method='average') #method="centroid" is cubic!
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    #Y = sch.linkage(D, method='single')
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    #dx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx1]
    axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    pylab.show()


def test():
    cluster_data(bed_fn)


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)
    bed_fn = sys.argv[1]
    cluster_data(bed_fn)
