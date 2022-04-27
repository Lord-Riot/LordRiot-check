## strange_classic_revenge

### 出题idea

> 在strange_classic中，主要应用的性质即是矩阵的特征值在相似变换下不变。于是在revenge中，想要更进一步，通过某些变换下不变的格的性质恢复格。

​		对于由格基矩阵 $B_{n\times n}$ 张成的格 $L$，考虑一个可逆矩阵 $U$, $B' = U \sdot B$，则以 $B'$ 为格基矩阵的格 $L'$ 为格 $L$ 的子格。子格的行列式必然为原格行列式的倍数，而如果取 $n$ 个线性无关的不同子格的向量 $v_i$，则以 $v_i$ 为基的格 $L''$ 同样为 $L$ 的子格。

​		理论上，只要泄露了足够多随机向量 $r_i \sdot B$ 的结果，则可以通过组合不同的子格 $L_i$ 求行列式，再求最大公因数来获得原格的行列式。之后可以不断的给子格中添加向量，并计算Hermite标准型，然后判断行列式是否等于原格行列式，从而得到原格的信息，但直接考这个未免太简单了。

​		考虑如果得到了与 $B$ 中部分向量正交的向量集，通过组合则可以获得以格 $L$ 为子格的格。例如，对于格基矩阵 $B = \begin{bmatrix}\vec b_0\\ \vec b_1\\ \vdots\\ \vec b_{n-1} \end{bmatrix}$ ，可以获得了 $B_0 = \begin{bmatrix}\vec b_0\\ \vec b_1\\ \vdots\\ \vec b_{j} \end{bmatrix}$ 的解空间的基 $\vec v_0, \vec v_1, \dots, \vec v_j$， 则可以求出 $\vec v_0, \vec v_1, \dots, \vec v_j$ 构成齐次方程的解空间的基 $\vec b_0', \dots, \vec b_j'$。同样对于 $B_1 = \begin{bmatrix}\vec b_{j+1}\\ \vec b_{j+2}\\ \vdots\\ \vec b_{n-1} \end{bmatrix}$ 可以求出其解空间的解空间的基 $\vec b_{j+1}', \dots, \vec b_{n-1}'$ 。以 $\vec b_0', \dots, \vec b_{n-1}' $ 为基构成的格 $\mathbb{L}$ 以 $L$ 为子格。由此，如果我们可以得到一个格 $L$ 的多个子格 $L_i$ 的解空间的基 $\vec v_0, \vec v_1, \dots, \vec v_{n-1}$，如果格 $L$ 的行列式为素数，则可以恢复出格 $L$。

​		如果有多个随机向量乘格基矩阵的结果 $r_i$， 给出 $f(r_i)$，需要利用 $f(r_i)$ 求出$r_i$ 对应解空间的基，再利用上述方法恢复出格 $L$ 。最终决定构造 $f(r_i, r_j) = (\sqrt{|r_{i,0}\sdot r_{j,0}|}, \sqrt{|r_{i,1}\sdot r_{j,1}|}, \dots, \sqrt{|r_{i,n-1}\sdot r_{j,n-1}|}) + \sum w_k r_k$ ，其中 $w_i$ 为随机数。以多个 $f(r_i, r_j)$ 构造垂直格，规约出的向量即落在了 $r_i, \dots, r_{n-1}$ 的解空间上，得到足够多个即得到了解空间，再应用上面方法，即可得到原格。

### exp

```python
from sage.all import *
from hashlib import sha1


Bits = 16
m = 32
K = 2


def get_orthogonal_basis(B):
    M = block_matrix(ZZ, [B, identity_matrix(B.nrows())], ncols=2)
    return M.LLL()[:m//2, -m:]


sub_lattice = []
cipher = 
for i in range(len(cipher)):
    mix0, mix1 = cipher[i]
    my_A0k = get_orthogonal_basis(Matrix(mix0).transpose())
    my_A1k = get_orthogonal_basis(Matrix(mix1).transpose())
    sub_lattice.append(block_matrix([Matrix(my_A0k.right_kernel().basis()), Matrix(my_A1k.right_kernel().basis())], nrows=2))

L2 = block_matrix(ZZ, sub_lattice, nrows=K)
print('MRCTF{%s}' % sha1(str(L2.hermite_form()[:m]).encode()).hexdigest())

```

