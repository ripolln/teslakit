import numpy as np
np.set_printoptions(precision=4)
import sys

WT_probs_fut = [0.19, 0.6, 0.21]
trans_matrix = [[0, .75, .25], [2/9.0, 4/9.0, 3/9.0], [1/6.0, 2/6.0, 3/6.0]]
# trans_matrix = [[0.2, 0.7, 0.1], [0.9, 0, 0.1], [0.2, 0.8, 0]]
num_clusters = 3
# print(trans_matrix)
# print()

#---------------------
A = np.append(np.transpose(trans_matrix)-np.identity(num_clusters), np.ones((1,num_clusters)), axis=0)
b = np.append(np.zeros((num_clusters,1)),[[1]])
WT_probs_his = np.linalg.solve(np.transpose(A).dot(A), np.transpose(A).dot(b))
print(WT_probs_his)
print()
#---------------------


alfa_p = np.identity(num_clusters)
beta_p = trans_matrix

alfa_q = [[beta_p[0][0]*WT_probs_fut[0], beta_p[1][0]*WT_probs_fut[1], beta_p[2][0]*WT_probs_fut[2]],
          [beta_p[0][1]*WT_probs_fut[0], beta_p[1][1]*WT_probs_fut[1], beta_p[2][1]*WT_probs_fut[2]],
          [beta_p[0][2]*WT_probs_fut[0], beta_p[1][2]*WT_probs_fut[1], beta_p[2][2]*WT_probs_fut[2]]]


beta_q = [[np.sum(alfa_q, axis=1)[0], 0, 0],
          [0, np.sum(alfa_q, axis=1)[1], 0],
          [0, 0, np.sum(alfa_q, axis=1)[2]]]


b = [0, 0, 0, WT_probs_fut[0]- beta_q[0][0],
              WT_probs_fut[1]- beta_q[1][1],
              WT_probs_fut[2]- beta_q[2][2]]


A = np.concatenate((alfa_p, beta_p), axis=1)
B = np.concatenate((alfa_q, beta_q), axis=1)
coef = np.concatenate((A,B), axis=0)


alfa_beta = np.linalg.solve(coef, b)
alfa = alfa_beta[:3]
beta = alfa_beta[3:]

trans_matrix_modif = np.zeros(np.shape(trans_matrix))*np.nan

for i in range(num_clusters):
     for j in range(num_clusters):
          trans_matrix_modif[i][j] = trans_matrix[i][j] * (1 + alfa[i] + beta[j])

# print(trans_matrix)
# print()
# print(trans_matrix_modif)



#---------------------
trans_matrix = trans_matrix_modif
A = np.append(np.transpose(trans_matrix)-np.identity(num_clusters), np.ones((1,num_clusters)), axis=0)
b = np.append(np.zeros((num_clusters,1)),[[1]])
WT_probs_his = np.linalg.solve(np.transpose(A).dot(A), np.transpose(A).dot(b))

print('-------')
print(WT_probs_fut)
print()
print(WT_probs_his)
print()


# kk = np.transpose(trans_matrix)-np.identity(num_clusters)
# A = np.append(kk[:-1,:], np.ones((1,num_clusters)), axis=0)
# b = np.append(np.zeros((2,1)),[[1]])
#
# WT_probs_his = np.linalg.solve(A, b)
#
# print(WT_probs_his)



