from sage.all import *
from Crypto.Util.number import *
# from secret import flag
from hashlib import sha1
from time import time
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


Bits = 16
m = 32
K = 2


def Rotor(x):
    return (vector(x) * random_matrix(ZZ, m, x=0, y=2**Bits) * L0).change_ring(ZZ).list()


def get_mixed(A):
    mix = []
    for j in range(len(A)):
        for jj in range(j + 1, len(A)):
            w = [getrandbits(200) for _ in range(len(A))]
            # mix.append(vector([int(sqrt(abs(A[j][ii] * A[jj][ii]))) + vector(w) * vector([A[jjj][ii] for jjj in range(len(A))]) for ii in range(len(A[0]))]))
            mix.append(vector([getRandomNBitInteger(16) + vector(w) * vector([A[jjj][ii] for jjj in range(len(A))]) for ii in range(len(A[0]))]))
    return mix[:m//2]


def get_orthogonal_basis(B):
    M = block_matrix(ZZ, [B, identity_matrix(B.nrows())], ncols=2)
    return M.LLL()[:m//2, -m:]


while True:
    while True:
        L0 = random_matrix(ZZ, m, x=0, y=2**Bits)
        if is_prime(L0.det()):
            break
    print('found L0')
    vectors = [Rotor(random_vector(ZZ, m, x=0, y=2**Bits).list()) for _ in range(m+K)]
    orthogonal_vectors = []
    sub_lattice = []
    for i in range(K):
        V0 = Matrix(Matrix(Matrix(ZZ, vectors[i:i+m//2]).right_kernel().basis()).right_kernel().basis())
        V1 = Matrix(Matrix(Matrix(ZZ, vectors[i+m//2:i+m]).right_kernel().basis()).right_kernel().basis())
        sub_lattice.append(block_matrix([V0, V1], nrows=2))
    L1 = block_matrix(ZZ, sub_lattice, nrows=K)
    if not L1.hermite_form()[:m] == L0.hermite_form():
        continue

    sub_lattice = []
    cipher = []
    for i in range(K):
        mix0 = get_mixed(vectors[i:i+m//2])
        my_A0k = get_orthogonal_basis(Matrix(mix0).transpose())
        mix1 = get_mixed(vectors[i+m//2:i+m])
        my_A1k = get_orthogonal_basis(Matrix(mix1).transpose())
        sub_lattice.append(block_matrix([Matrix(my_A0k.right_kernel().basis()), Matrix(my_A1k.right_kernel().basis())], nrows=2))
        cipher.append([mix0, mix1])

    L2 = block_matrix(ZZ, sub_lattice, nrows=K)
    assert L2.hermite_form()[:m] == L0.hermite_form()

    flag = 'MRCTF{%s}' % sha1(str(L0.hermite_form()).encode()).hexdigest()
    assert flag.lstrip('MRCTF{').rstrip('}') == sha1(str(L0.hermite_form()).encode()).hexdigest()
    open('cipher.txt', 'w').write('cipher = ' + str(cipher))
    print(flag)
    # print(t)
    break
