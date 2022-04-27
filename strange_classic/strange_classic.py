from sage.all import *
from Crypto.Util.number import *
from secret import flag, t, vectors, find_G
from hashlib import sha1


N = 43
m = 9
while True:
    F = random_matrix(GF(N), m)
    if F.is_invertible():
        break


def Rotor(x):
    y = x.copy()
    for i in range(m):
        y.append(int(sum([y[j+i] * t[j] % N for j in range(m)]) % N))
    return [int(vi) for vi in (vector(y[-m:]) * G_).list()]


def Reflector(x, inv=False):
    return (F.inverse() * vector(x)).change_ring(ZZ).list() if inv else (F * vector(x)).change_ring(ZZ).list()


G = find_G(t)
C = Matrix(GF(N), t[1:]).transpose()
tmp_matrix = block_matrix([identity_matrix(len(t) - 1), C], ncols=2, subdivdie=False)
C = block_matrix([Matrix([0] * (len(t) - 1) + [t[0]]), tmp_matrix], nrows=2, subdivdie=False)
assert G * C == C * G
check = G.eigenvalues()
for eigenvalue in check:
    assert eigenvalue in GF(N)
G_ = G ** (N * (N-1)//3)

assert flag.lstrip('MRCTF{').rstrip('}') == sha1(str(t).encode()).hexdigest()
vectors += [random_vector(GF(N), m).change_ring(ZZ).list() for _ in range(len(vectors)*4)]
shuffle(vectors)

print('plain =', [vec[:m//2] for vec in vectors])
print('cipher =', [Reflector(Rotor(Reflector(vi)), inv=True)[:m//2] for vi in vectors])
print('hash =', sha1(flag.encode()).hexdigest())
