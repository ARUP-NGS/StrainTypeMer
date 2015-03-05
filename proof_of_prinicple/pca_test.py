__author__ = 'ksimmon'
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

kmer_data = []

for line in open('/home/ksimmon/Documents/atrib_1000_remove_shared.txt', "r"):
    kmer_data.append( map(int, line.strip().split(",")))



# np.random.seed(4294967295)
#
# mu_vec1 = np.array([0,0,0])
# cov_mat1 = np.array([[1,0,0],[0,1,0],[0,0,1]])
# class1_sample = np.random.multivariate_normal(mu_vec1, cov_mat1, 20).T
#
# assert class1_sample.shape == (3,20), "The matrix has not the dimensions 3X20"
#
# mu_vec2 = np.array([1,1,1])
# cov_mat2 = np.array([[1,0,0],[0,1,0],[0,0,1]])
# class2_sample = np.random.multivariate_normal(mu_vec2, cov_mat2, 20).T
# assert class1_sample.shape == (3,20), "The matrix has not the dimensions 3x20"
#
# print class1_sample
# print class2_sample
#
#
#
#
#
#
# fig = plt.figure(figsize=(8,8))
# ax = fig.add_subplot(111, projection='3d')
# plt.rcParams['legend.fontsize'] = 10
# ax.plot(class1_sample[0,:], class1_sample[1,:],\
#     class1_sample[2,:], 'o', markersize=8, color='blue', alpha=0.5, label='class1')
# ax.plot(class2_sample[0,:], class2_sample[1,:],\
#     class2_sample[2,:], '^', markersize=8, alpha=0.5, color='red', label='class2')
#
# plt.title('Samples for class 1 and class 2')
# ax.legend(loc='upper right')
# #plt.show()
#
# all_samples = np.concatenate((class1_sample, class2_sample), axis=1)
# assert all_samples.shape == (3,40), "The matrix has not the dimensions 3x40"
#
# mean_x = np.mean(all_samples[0,:])
# mean_y = np.mean(all_samples[1,:])
# mean_z = np.mean(all_samples[2,:])
#
# mean_vector = np.array([[mean_x],[mean_y],[mean_z]])
#
# print('Mean Vector:\n', mean_vector)
#
# scatter_matrix = np.zeros((3,3))
# for i in range(all_samples.shape[1]):
#     scatter_matrix += (all_samples[:,i].reshape(3,1)\
#         - mean_vector).dot((all_samples[:,i].reshape(3,1) - mean_vector).T)
# print('Scatter Matrix:\n', scatter_matrix)
#
# cov_mat = np.cov([all_samples[0,:],all_samples[1,:],all_samples[2,:]])
# print('Covariance Matrix:\n', cov_mat)
#
#
#
# # eigenvectors and eigenvalues for the from the scatter matrix
# eig_val_sc, eig_vec_sc = np.linalg.eig(scatter_matrix)
#
# # eigenvectors and eigenvalues for the from the covariance matrix
# eig_val_cov, eig_vec_cov = np.linalg.eig(cov_mat)
#
# for i in range(len(eig_val_sc)):
#     eigvec_sc = eig_vec_sc[:,i].reshape(1,3).T
#     eigvec_cov = eig_vec_cov[:,i].reshape(1,3).T
#     assert eigvec_sc.all() == eigvec_cov.all(), 'Eigenvectors are not identical'
#
#     print('Eigenvector {}: \n{}'.format(i+1, eigvec_sc))
#     print('Eigenvalue {} from scatter matrix: {}'.format(i+1, eig_val_sc[i]))
#     print('Eigenvalue {} from covariance matrix: {}'.format(i+1, eig_val_cov[i]))
#     print('Scaling factor: ', eig_val_sc[i]/eig_val_cov[i])
#     print(40 * '-')
#
#
#
#
#
#
#
#
#
#
#
#
# from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d import proj3d
# from matplotlib.patches import FancyArrowPatch
#
#
# class Arrow3D(FancyArrowPatch):
#     def __init__(self, xs, ys, zs, *args, **kwargs):
#         FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
#         self._verts3d = xs, ys, zs
#
#     def draw(self, renderer):
#         xs3d, ys3d, zs3d = self._verts3d
#         xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
#         self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
#         FancyArrowPatch.draw(self, renderer)
#
# fig = plt.figure(figsize=(7,7))
# ax = fig.add_subplot(111, projection='3d')
#
# ax.plot(all_samples[0,:], all_samples[1,:],\
#     all_samples[2,:], 'o', markersize=8, color='green', alpha=0.2)
# ax.plot([mean_x], [mean_y], [mean_z], 'o', \
#     markersize=10, color='red', alpha=0.5)
# for v in eig_vec_sc.T:
#     a = Arrow3D([mean_x, v[0]], [mean_y, v[1]],\
#         [mean_z, v[2]], mutation_scale=20, lw=3, arrowstyle="-|>", color="r")
#     ax.add_artist(a)
# ax.set_xlabel('x_values')
# ax.set_ylabel('y_values')
# ax.set_zlabel('z_values')
#
# plt.title('Eigenvectors')
#
# for ev in eig_vec_sc:
#     np.testing.assert_array_almost_equal(1.0, np.linalg.norm(ev))
#     # instead of 'assert' because of rounding errors
#
# # Make a list of (eigenvalue, eigenvector) tuples
# eig_pairs = [(np.abs(eig_val_sc[i]), eig_vec_sc[:,i]) for i in range(len(eig_val_sc))]
#
# # Sort the (eigenvalue, eigenvector) tuples from high to low
# eig_pairs.sort()
# eig_pairs.reverse()
#
# # Visually confirm that the list is correctly sorted by decreasing eigenvalues
# for i in eig_pairs:
#     print(i[0])
#
# matrix_w = np.hstack((eig_pairs[0][1].reshape(3,1), eig_pairs[1][1].reshape(3,1)))
# print('Matrix W:\n', matrix_w)
#
#
# transformed = matrix_w.T.dot(all_samples)
# assert transformed.shape == (2,40), "The matrix is not 2x40 dimensional."
#
#
# plt.plot(transformed[0,0:20], transformed[1,0:20],\
#      'o', markersize=7, color='blue', alpha=0.5, label='class1')
# plt.plot(transformed[0,20:40], transformed[1,20:40],
#      '^', markersize=7, color='red', alpha=0.5, label='class2')
# plt.xlim([-4,4])
# plt.ylim([-4,4])
# plt.xlabel('x_values')
# plt.ylabel('y_values')
# plt.legend()
# plt.title('Transformed samples with class labels')

#plt.show()


from matplotlib.mlab import PCA as mlabPCA

kmr_nd = np.array(kmer_data)

print kmr_nd[0]

#all_samples = np.concatenate((class1_sample, class2_sample), axis=1)



print kmr_nd.shape
mlab_pca = mlabPCA(kmr_nd)


print mlab_pca.Y.shape

# print('PC axes in terms of the measurement axes'\
#         ' scaled by the standard deviations:\n',\
#           mlab_pca.Wt)


_key = { 1: "A1",
          2 : "A2",
         3: "A3",
         4: "A4_2",
         5: "A4",
         6: "A5_2",
         7: "A5",
         8 :"A6",
         9 : "A7",
         10 : "A8",
         11 : "A9",
         12 : "AC1",
         13 : "AC2",
         14 : "AC3",
         15 : "AC4",
         16 : "AC5",
         17 : "AC6",
         18 : "AC7",
         19 : "AC8",
         20 : "AC9_2",
         21 : "AC9" }
for i in range(21):
    plt.plot(mlab_pca.Y[i:i+1,0],mlab_pca.Y[i:i+1,1], 'o', markersize=10,\
            color='blue', alpha=0.2, label='class1')
    plt.text(mlab_pca.Y[i:i+1,0],mlab_pca.Y[i:i+1,1], _key[i+1], fontsize=16)


# plt.plot(mlab_pca.Y[20:40,0], mlab_pca.Y[20:40,1], '^', markersize=7,\
#          color='red', alpha=0.5, label='class2')

plt.xlabel('x_values')
plt.ylabel('y_values')
plt.xlim([-4,4])
plt.ylim([-4,4])
#plt.legend()

plt.show()