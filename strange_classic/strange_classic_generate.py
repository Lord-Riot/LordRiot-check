from sage.all import *
from Crypto.Util.number import *
# from secret import flag, t, vectors
from hashlib import sha1
from time import time
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


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
    return [int(vi) for vi in (vector(x) * F.inverse())] if inv else [int(vi) for vi in (vector(x) * F)]


def leak_part(vec, num):
    return vec[:num]


def find_random_roots():
    roots = []
    while len(roots) < m:
        tmp = getRandomRange(1, N - 1)
        if tmp in roots:
            continue
        else:
            roots.append(tmp)
    return roots


def find_GC():
    guess_num = m//2
    PR = PolynomialRing(GF(N), names=['x%d_%d' % (i, j) for i in range(m) for j in range(m)])
    PRx = PolynomialRing(GF(N), 'x')
    x = PRx.gens()[0]
    a = PR.gens()

    while True:
        t_roots = find_random_roots()
        while True:
            fx = 1
            for i in range(m):
                fx *= (x - t_roots[i])
            if len(fx.coefficients()) == m + 1:
                break
        t = [int(-i % N) for i in fx.coefficients()[:-1]]
        C = Matrix(GF(N), t[1:]).transpose()
        tmp_matrix = block_matrix([identity_matrix(len(t) - 1), C], ncols=2, subdivdie=False)
        C = block_matrix([Matrix([0] * (len(t) - 1) + [t[0]]), tmp_matrix], nrows=2, subdivdie=False)
        M = Matrix([[a[i+j*m] for i in range(m)] for j in range(m)])
        judge = M*C - C*M
        Ideal = ideal(judge.list())
        GB = Ideal.groebner_basis()
        Xma = [a[i] - GB[i] for i in range(len(GB))]
        Xm = Xma + list(a)[-m:]
        M = Matrix([[Xm[i+j*m] for i in range(m)] for j in range(m)])
        char_poly_coefficients = M.characteristic_polynomial().coefficients()[:-1][::-1]

        PRR = PolynomialRing(GF(N), names=['y%d' % i for i in range(m)] + ['y'])
        y = PRR.gens()
        while True:
            expected_fx = GF(N)(1)
            roots = find_random_roots()
            for i in range(m):
                expected_fx *= (y[-1] - roots[i])
            expected_poly_coefficients = expected_fx.coefficients()[1:]
            if len(expected_poly_coefficients) == m:
                break
        ok = False
        for guess in range(N**guess_num):
            guess_y = []
            for i in range(guess_num):
                guess_y.append(guess % N)
                guess //= N
            condition_poly = [char_poly_coefficients[i](a[:-m]+y[:-1]) - expected_poly_coefficients[i] for i in range(m)]
            Ideal = ideal(condition_poly+[y[i] - guess_y[i] for i in range(guess_num)])
            GB = Ideal.groebner_basis()
            if GB != [1]:
                print(GB, str(t))
                ok = True
                break
        if ok:
            Mx = [GF(N)(GB[i]-y[i]) for i in range(m)]
            MT = M([0] * m * (m-1) + Mx)
            if MT * C == C * MT:
                return MT, C, t_roots, True
            else:
                return GB, C, t_roots, False


def is_k(v0, v1):
    judge = v0[0] * inverse(v1[0], N) % N
    for w in range(1, len(v1)):
        if v0[w] * inverse(v1[w], N) % N != judge:
            return False
    return judge


def find_k(plain_texts, cipher_texts):
    k = []
    for i in range(len(plain_texts)):
        if plain_texts[i] in cipher_texts and cipher_texts[i] in plain_texts:
            xi = plain_texts[cipher_texts.index(plain_texts[i])]
            xj = cipher_texts[plain_texts.index(cipher_texts[i])]
            tmp_k = is_k(xj, xi)
            if tmp_k:
                k.append(tmp_k)
    return k


PR = PolynomialRing(GF(N), 'x')
x = PR.gens()[0]
PRR = PolynomialRing(GF(N), names=['x%d_%d' % (i, j) for i in range(m) for j in range(m)])
a = PRR.gens()
# while True:
#     G, C, roots, judge = find_GC()
#     if judge:
#         break
roots = [31, 17, 11, 39, 24, 1, 25, 3, 5]
fx = 1
for i in range(m):
    fx *= (x - roots[i])
t = [int(-i % N) for i in fx.coefficients()[:-1]]
C = Matrix(GF(N), t[1:]).transpose()
tmp_matrix = block_matrix([identity_matrix(len(t) - 1), C], ncols=2, subdivdie=False)
C = block_matrix([Matrix([0] * (len(t) - 1) + [t[0]]), tmp_matrix], nrows=2, subdivdie=False)
M = Matrix([[a[i+j*m] for i in range(m)] for j in range(m)])
judge = M*C - C*M
Ideal = ideal(judge.list())
GB = Ideal.groebner_basis()
Xma = [a[i] - GB[i] for i in range(len(GB))]
Xm = Xma + list(a)[-m:]
M = Matrix([[Xm[i+j*m] for i in range(m)] for j in range(m)])
Mx = [1, 0, 0, 0, 8, 7, 8, 3, 5]
G = M([0] * m * (m-1) + Mx)

n = 3
G_ = G ** (N * (N-1)//n)
assert G_ * C == C * G_
assert G_ ** n == identity_matrix(GF(N), m)

index = n * m
target_matrix = C**index
vectors = []
for i in range(m):
    vectors.append((getRandomRange(1, N - 1) * (target_matrix - roots[i]**index * identity_matrix(GF(N), m)).kernel().basis()[0]))
vectors = [Reflector(vi, inv=True) for vi in vectors]
real_vectors = []
for i in range(m):
    v0 = Reflector(Rotor(Reflector(vectors[i])), inv=True)
    v1 = Reflector(Rotor(Reflector(v0)), inv=True)
    real_vectors.append(v0)
    real_vectors.append(v1)
init_vectors = vectors
vectors = vectors + real_vectors
vectors += [random_vector(GF(N), m).change_ring(ZZ).list() for _ in range(len(vectors)*4)]
shuffle(vectors)

flag = 'MRCTF{%s}' % sha1(str(t).encode()).hexdigest()
assert flag.lstrip('MRCTF{').rstrip('}') == sha1(str(t).encode()).hexdigest()
print('plain =', [vec[:m//2] for vec in vectors])
print('cipher =', [Reflector(Rotor(Reflector(vi)), inv=True)[:m//2] for vi in vectors])
print('hash =', sha1(flag.encode()).hexdigest())

all_k = list(set(find_k([vec[:m//2] for vec in vectors],  [Reflector(Rotor(Reflector(vi)), inv=True)[:m//2] for vi in vectors])))
if len(all_k) == m:
    print('Ok')
