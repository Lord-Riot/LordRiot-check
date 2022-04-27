from Crypto.Util.number import *
from time import time
from random import *
from time import time
from hashlib import sha256


def h2p(hex_str, E):
    return E((int(hex_str[:64], 16), int(hex_str[64:], 16)))


def mov_atk(G, xG, T, n):
    e0 = T.weil_pairing(G, n)
    e1 = T.weil_pairing(xG, n)
    return e1.log(e0)


def ph_mov_atk(p, a, b, G, xG, limit):
    N = G.order()
    for k in range(2, 6):
        if (p**k-1) % N == 0:
            break
    Ek = EllipticCurve(GF(p**k, 'w'), [a, b])
    EkN = Ek.gens()[0].order()
    factor_list = N.factor()
    a = []
    n = []
    judge = 1
    for pi, index in factor_list:
        tmp_n = pi ** index
        if tmp_n < 2**10:
            tmp_a = bsgs(N // tmp_n * G, N // tmp_n * xG, (2**5, 2**10), operation='+')
        else:
            while True:
                T = EkN//tmp_n * Ek.random_element()
                if T.order() == tmp_n:
                    break
            tmp_a = mov_atk(Ek(N // tmp_n * G), Ek(N // tmp_n * xG), T, tmp_n)
        n.append(tmp_n)
        a.append(tmp_a)
        judge *= tmp_n
        if judge > limit:
            break
    return CRT(a, n)


def small_nonce_attack(_n, h, s, r, bits, block_size=20, private_key=None, alpha1=None, alpha2=None):
    num = len(h)
    if alpha1 is None:
        alpha2 = []
        alpha1 = []
        for j in range(num):
            alpha2.append(int(h[j] * inverse(s[j], _n) % _n))
            alpha1.append(int(r[j] * inverse(s[j], _n) % _n))

    adjust = 2 ** (_n.bit_length() - bits)
    add_adjust = -6 * 2 ** (_n.bit_length() - 2)
    L = Matrix(ZZ, num + 2)
    for j in range(num):
        L[j, j] = _n * adjust
        L[-2, j] = alpha1[j] * adjust
        L[-1, j] = alpha2[j] * adjust + add_adjust
    L[-2, num] = 1
    L[-1, num + 1] = _n

    start_time = time()

    if block_size == 0:
        LL = L.LLL()
    else:
        LL = L.BKZ(block_size=block_size)
    end_time = time()
    print("use total %d seconds to reduce Lattice" % int(end_time - start_time))

    if debug:
        T = [int(alpha2[j] + alpha1[j] * private_key) // _n for j in range(num)]
        expected_sol = Matrix(ZZ, [T + [-private_key, -1]])
        expected_res = (expected_sol * L)[0]
        print('expected vector:', list(expected_res))
        print('expected norm:', int(expected_res.norm()).bit_length(), int(expected_res.norm()))
        j = 0
        for i in expected_res[:num+2]:
            if i > 0:
                j += 1
        print('balance: %d/%d' % (j, num+2))
    for vec in LL:
        if vec[0] == 0:
            continue
        if abs(vec[num + 1]) == _n:
            possible_private_keys = [abs(vec[num]), _n - abs(vec[num])]
            for possible_private_key in possible_private_keys:
                if debug:
                    if possible_private_key == private_key:
                        print('Got it', possible_private_key)
                        return
                if possible_private_key.nbits() == 256:
                    print(possible_private_key)
                    return
            print("Attacked Unsuccessfully")
    return


class Curve:
    def __init__(self, ecc_table, private_key=None, public_key=None):
        self.ecc_table = ecc_table
        self.para_len = len(ecc_table['p'])
        self.ecc_a3 = (int(ecc_table['a'], base=16) + 3) % int(ecc_table['p'], base=16)
        if 'g' not in ecc_table.keys():
            while True:
                gx = getPrime(127)
                gy = sqrt_mod(gx ** 3 + gx, int(self.ecc_table['p'], 16))
                if gy is None:
                    continue
                g = hex(gx)[2:].rjust(64, '0') + hex(gy)[2:].rjust(64, '0')
                if self.kg((int(self.ecc_table['p'], 16) - 1) // 2, g) is not None:
                    self.ecc_table['g'] = g
                    break

        if private_key is None:
            private_key = getRandomRange(int(self.ecc_table['p'], 16)//2, int(self.ecc_table['p'], 16)-1)
        if public_key is None:
            public_key = self.kg(private_key, self.ecc_table['g'])
        self.private_key = private_key
        self.public_key = public_key

    def kg(self, k, Point):
        Point = '%s%s' % (Point, '1')
        mask_str = '8'
        for i in range(self.para_len - 1):
            mask_str += '0'
        mask = int(mask_str, 16)
        Temp = Point
        flag = False
        for n in range(self.para_len * 4):
            if flag:
                Temp = self._double_point(Temp)
            if (k & mask) != 0:
                if flag:
                    Temp = self._add_point(Temp, Point)
                else:
                    flag = True
                    Temp = Point
            k = k << 1
        return self._convert_jacb_to_nor(Temp)

    def _double_point(self, Point):
        l = len(Point)
        len_2 = 2 * self.para_len
        if l < self.para_len * 2:
            return None
        else:
            x1 = int(Point[0:self.para_len], 16)
            y1 = int(Point[self.para_len:len_2], 16)
            if l == len_2:
                z1 = 1
            else:
                z1 = int(Point[len_2:], 16)

            T6 = (z1 * z1) % int(self.ecc_table['p'], base=16)
            T2 = (y1 * y1) % int(self.ecc_table['p'], base=16)
            T3 = (x1 + T6) % int(self.ecc_table['p'], base=16)
            T4 = (x1 - T6) % int(self.ecc_table['p'], base=16)
            T1 = (T3 * T4) % int(self.ecc_table['p'], base=16)
            T3 = (y1 * z1) % int(self.ecc_table['p'], base=16)
            T4 = (T2 * 8) % int(self.ecc_table['p'], base=16)
            T5 = (x1 * T4) % int(self.ecc_table['p'], base=16)
            T1 = (T1 * 3) % int(self.ecc_table['p'], base=16)
            T6 = (T6 * T6) % int(self.ecc_table['p'], base=16)
            T6 = (self.ecc_a3 * T6) % int(self.ecc_table['p'], base=16)
            T1 = (T1 + T6) % int(self.ecc_table['p'], base=16)
            z3 = (T3 + T3) % int(self.ecc_table['p'], base=16)
            T3 = (T1 * T1) % int(self.ecc_table['p'], base=16)
            T2 = (T2 * T4) % int(self.ecc_table['p'], base=16)
            x3 = (T3 - T5) % int(self.ecc_table['p'], base=16)

            if (T5 % 2) == 1:
                T4 = (T5 + ((T5 + int(self.ecc_table['p'], base=16)) >> 1) - T3) % int(self.ecc_table['p'], base=16)
            else:
                T4 = (T5 + (T5 >> 1) - T3) % int(self.ecc_table['p'], base=16)

            T1 = (T1 * T4) % int(self.ecc_table['p'], base=16)
            y3 = (T1 - T2) % int(self.ecc_table['p'], base=16)

            form = '%%0%dx' % self.para_len
            form = form * 3
            return form % (x3, y3, z3)

    def _add_point(self, P1, P2):
        len_2 = 2 * self.para_len
        l1 = len(P1)
        l2 = len(P2)
        if (l1 < len_2) or (l2 < len_2):
            return None
        else:
            X1 = int(P1[0:self.para_len], 16)
            Y1 = int(P1[self.para_len:len_2], 16)
            if l1 == len_2:
                Z1 = 1
            else:
                Z1 = int(P1[len_2:], 16)
            x2 = int(P2[0:self.para_len], 16)
            y2 = int(P2[self.para_len:len_2], 16)

            T1 = (Z1 * Z1) % int(self.ecc_table['p'], base=16)
            T2 = (y2 * Z1) % int(self.ecc_table['p'], base=16)
            T3 = (x2 * T1) % int(self.ecc_table['p'], base=16)
            T1 = (T1 * T2) % int(self.ecc_table['p'], base=16)
            T2 = (T3 - X1) % int(self.ecc_table['p'], base=16)
            T3 = (T3 + X1) % int(self.ecc_table['p'], base=16)
            T4 = (T2 * T2) % int(self.ecc_table['p'], base=16)
            T1 = (T1 - Y1) % int(self.ecc_table['p'], base=16)
            Z3 = (Z1 * T2) % int(self.ecc_table['p'], base=16)
            T2 = (T2 * T4) % int(self.ecc_table['p'], base=16)
            T3 = (T3 * T4) % int(self.ecc_table['p'], base=16)
            T5 = (T1 * T1) % int(self.ecc_table['p'], base=16)
            T4 = (X1 * T4) % int(self.ecc_table['p'], base=16)
            X3 = (T5 - T3) % int(self.ecc_table['p'], base=16)
            T2 = (Y1 * T2) % int(self.ecc_table['p'], base=16)
            T3 = (T4 - X3) % int(self.ecc_table['p'], base=16)
            T1 = (T1 * T3) % int(self.ecc_table['p'], base=16)
            Y3 = (T1 - T2) % int(self.ecc_table['p'], base=16)

            form = '%%0%dx' % self.para_len
            form = form * 3
            return form % (X3, Y3, Z3)

    def _convert_jacb_to_nor(self, Point):
        len_2 = 2 * self.para_len
        x = int(Point[0:self.para_len], 16)
        y = int(Point[self.para_len:len_2], 16)
        z = int(Point[len_2:], 16)
        z_inv = pow(z, int(self.ecc_table['p'], base=16) - 2, int(self.ecc_table['p'], base=16))
        z_invSquar = (z_inv * z_inv) % int(self.ecc_table['p'], base=16)
        z_invQube = (z_invSquar * z_inv) % int(self.ecc_table['p'], base=16)
        x_new = (x * z_invSquar) % int(self.ecc_table['p'], base=16)
        y_new = (y * z_invQube) % int(self.ecc_table['p'], base=16)
        z_new = (z * z_inv) % int(self.ecc_table['p'], base=16)
        if z_new == 1:
            form = '%%0%dx' % self.para_len
            form = form * 2
            return form % (x_new, y_new)
        else:
            return None

    def sign(self, data, K):
        E = data.hex()
        e = int(E, 16)

        d = self.private_key
        k = K

        P1 = self.kg(k, self.ecc_table['g'])

        x = int(P1[0:self.para_len], 16)
        R = ((e + x) % int(self.ecc_table['n'], base=16))
        if R == 0 or R + k == int(self.ecc_table['n'], base=16):
            return None
        d_1 = pow(d + 1, int(self.ecc_table['n'], base=16) - 2, int(self.ecc_table['n'], base=16))
        S = (d_1 * (k + R) - R) % int(self.ecc_table['n'], base=16)
        if S == 0:
            return None
        else:
            return '%064x%064x' % (R, S)


# flag1
for i in range(1, 5):
    p = int(input().strip(), 16)
    E = EllipticCurve(GF(p), [1, 0])
    G = h2p(input().strip(), E)
    xG = h2p(input().strip(), E)
    print(ph_mov_atk(p, 1, 0, G, xG, 2**(30*i)))

# flag2
num = 60
h = [bytes_to_long(sha256(b'Lord Riot %d' % i).digest()) for i in range(num)]
my_ecc_table = {
    'n': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123',
    'p': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF',
    'g': '32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7'
         'bc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0',
    'a': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC',
    'b': '28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93',
}
SM2 = Curve(ecc_table=my_ecc_table)
power = 4
debug = False
if debug:
    samples = [SM2.sign(sha256(b'Lord Riot %d' % i).digest(), randrange(2 ** (255 - power), 2 ** (256 - power))) for i in range(num)]
    r = [int(i[:64], 16) for i in samples]
    s = [int(i[64:], 16) for i in samples]
else:
    r = []
    s = []
    for i in range(num):
        data = input()
        data = data[data.index(':')+2:]
        r.append(int(data[:64], 16))
        s.append(int(data[64:], 16))

n = int(my_ecc_table['n'], 16)
if debug:
    small_nonce_attack(n, h, r, s, 255 - power, block_size=35, private_key=SM2.private_key,
                       alpha1=[s[j] + r[j] for j in range(num)],
                       alpha2=[s[j] for j in range(num)])
    print(SM2.private_key)
else:
    small_nonce_attack(n, h, r, s, 255 - power, block_size=30,
                       alpha1=[s[j] + r[j] for j in range(num)],
                       alpha2=[s[j] for j in range(num)])
