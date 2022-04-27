from Crypto_tools import *
import traceback


# display matrix picture with 0 and X
def matrix_overview(BB, bound):
    for ii in range(BB.dimensions()[0]):
        a = ('%02d ' % ii)
        for jj in range(BB.dimensions()[1]):
            a += '0' if BB[ii, jj] == 0 else 'X'
            if BB.dimensions()[0] < 60:
                a += ' '
        if BB[ii, ii] >= bound:
            a += '~'
        print(a)


beta = 0.32
delta = 0.04
amplification = 2048

# e = 2953544268002866703872076551930953722572317122777861299293407053391808199220655289235983088986372630141821049118015752017412642148934113723174855236142887
# N = 6006128121276172470274143101473619963750725942458450119252491144009018469845917986523007748831362674341219814935241703026024431390531323127620970750816983
# c = 4082777468662493175049853412968913980472986215497247773911290709560282223053863513029985115855416847643274608394467813391117463817805000754191093158289399


Xp = int(N**(delta + beta))
Yp = int(N**beta)
Yq = N//Yp
assert p < Xp and dp < Yp

modulus = e
mm = 5
ss = 0
tt = 3

P.<x, y, z> = PolynomialRing(ZZ)
Q = P.quotient(N - y * z)
pol = x * (N - y) + N
pol = Q(pol).lift()

# x-z-shifts
gg = []
monomials = []
for ii in range(mm + 1):
    for jj in range(mm - ii + 1):
        x_z_shift = z ^ ss * x ^ jj * modulus ^ (mm - ii) * pol ^ ii
        x_z_shift = Q(x_z_shift).lift()
        gg.append(x_z_shift)

# y-z-shifts (selected by Herrman and May)
for ii in range(mm + 1):
    for jj in range(1, tt + 1):
        y_z_shift = z ^ ss * y ^ jj * pol ^ ii * modulus ^ (mm - ii)
        y_z_shift = Q(y_z_shift).lift()
        gg.append(y_z_shift)

# list of monomials
for polynomial in gg:
    for monomial in polynomial.monomials():
        if monomial not in monomials:
            monomials.append(monomial)


print(monomials)
print('N =', N)
print('e =', e)

# construct lattice B
nn = len(monomials)
BB = Matrix(ZZ, nn)
for ii in range(nn):
    for jj in range(0, nn):
        if monomials[jj] in gg[ii].monomials():
            BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](Xp, Yp, Yq)

matrix_overview(BB, modulus ^ mm)

det = abs(BB.det())
bound = modulus ^ (mm * nn)
print('Bound check:', det < bound)
print(int(det).bit_length(), int(bound).bit_length())

# LLL
BB = BB.LLL()
print('LLL done')
matrix_overview(BB, modulus ^ mm)

PR.<xp, yp, zp> = PolynomialRing(ZZ)
PRQ = PR.quotient(N - yp * zp)
all_pol = []

for pol1_idx in tqdm(range(nn)):
    pol1 = 0
    for jj in range(nn):
        pol1 += monomials[jj](xp, yp, zp) * BB[pol1_idx, jj] / monomials[jj](Xp, Yp, Yq)
    all_pol.append(pol1)

# for i in range(4):
#     print(all_pol[i])
# tmp = all_pol[0].resultant(all_pol[1])
# PRR.<w> = PolynomialRing(ZZ)
# print(tmp(w, w, w).roots())

I = ideal(all_pol[:5])
GB = I.groebner_basis()
print('Groebner basis:')
print(GB)
print('-' * 32)

xv, yv, zv = var("xp,yp,zp")
print('roots:')
res = solve([h_i(xv, yv, zv) for h_i in GB], xv, yv, zv)

PRRR.<w> = PolynomialRing(ZZ)
for part_res in res:
    then_res = PRRR(part_res[1](w))
    p = abs(then_res.coefficients()[0].numerator())
    q = N // p
    if p == 1 and q == 1:
        continue
    assert p * q == N
    print(p)
    print(q)
    print(long_to_bytes(pow(c, inverse(e, (p-1)*(q-1)), N)))