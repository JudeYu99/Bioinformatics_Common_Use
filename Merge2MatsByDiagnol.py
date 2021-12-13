# @Author: Yu Zhu
# @Email: yzhu99@stu.suda.edu.cn

# @ Usage:
#   Merge two square matrices by diagonal and set diagonal elements to 0.
#   Numpy package required!


def Merge2MatsByDiagnol(matrix1, matrix2):
    
    import numpy as np
    #import numpy.matlib
    
    if matrix1.shape != matrix2.shape:

        return "Sizes of the two matrices are not same!"

    else:
        #new_mat = np.matlib.empty(matrix1.shape)
        
        new_mat = matrix2
        new_mat[range(matrix2.shape[0]), range(matrix2.shape[0])] = 0
        
        for i in range(matrix1.shape[0]):
            for j in range(i):
                new_mat[i, j] = matrix1[i, j]

        return new_mat
    
    
if __name__ == '__main__':

    a = np.arange(25).reshape(5,5)

    b = np.arange(100, 125).reshape(5,5)

    Merge2MatsByDiagnol(a, b)
