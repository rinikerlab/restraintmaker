from matplotlib import pyplot as plt


# PLOT DATA
def _plot_data(coords, eigen_vectors, title):
    ##original data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(coords[:, 0], coords[:, 1], c="k", label="data")
    ax.plot((-eigen_vectors[0, 0], eigen_vectors[0, 0]), (-eigen_vectors[0, 1], eigen_vectors[0, 1]), label="pc1")
    ax.plot((- eigen_vectors[1, 0], eigen_vectors[1, 0]), (-eigen_vectors[1, 1], eigen_vectors[1, 1]), label="pc2")
    ax.legend()
    ax.set_title("orig data - " + title)

    fig.set_tight_layout("thight")
    fig.show()


def _plot_pc_projection(projected_coords, title):
    ##porjection of data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(projected_coords[0], projected_coords[1], c="k", label="data")
    ax.scatter(projected_coords[0], [0 for x in projected_coords[0]], c="r", label="prPC1")
    ax.scatter([0 for x in projected_coords[1]], projected_coords[1], c="b", label="prPC2")
    ax.legend()
    ax.set_title("pca_projection - " + title)

    fig.set_tight_layout("thight")
    fig.show()


def _plot_2D_PCA(coords, selected_arr, mol_name="", eig_vecs=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(coords[:, 0], coords[:, 1], c="k", alpha=0.4)
    ax.scatter(selected_arr[:, 0], selected_arr[:, 1], c="r", label="selected")

    if (not eig_vecs):
        for ind, vec in enumerate(eig_vecs):
            ax.plot((vec[0], vec[1]), label="pca" + str(ind))

    ax.set_xlabel('X ')
    ax.set_ylabel('Y ')
    ax.set_title("2DPCA " + mol_name + " molecule")
    ax.legend()

    fig.set_tight_layout("thight")
    fig.show()


def _plot_3D_PCA(coords, eig_vecs=None, title=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ##vis mol
    ax.scatter(list(coords[:, 0]), list(coords[:, 1]),
               list(coords[:, 2]), c="k", alpha=0.4)
    ##vis selected
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c="r", label="selected")

    ##vis pca vector
    if (type(eig_vecs) != None and len(eig_vecs) > 0):
        for ind, vec in enumerate(eig_vecs):
            ax.plot(xs=[0, 1 * vec[0]], ys=[0, 1 * vec[1]], zs=[0, 1 * vec[2]], label="pca" + str(ind))

    ##vis tuning
    ax.set_xlabel('X ')
    ax.set_ylabel('Y ')
    ax.set_zlabel('Z ')
    ax.set_title("3DPCA " + title + " molecule")
    # ax.legend()

    fig.set_tight_layout("thight")
    fig.show()
