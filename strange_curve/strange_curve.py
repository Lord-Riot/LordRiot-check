from Crypto.Util.number import *
from sympy import sqrt_mod
from hashlib import sha256
from string import ascii_letters
from signal import alarm
from private import flag1, flag2, get_prime


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


def proof_of_work():
    m = [ascii_letters[getRandomRange(1, len(ascii_letters))] for _ in range(20)]
    check = sha256(''.join(m).encode()).hexdigest()
    print('SHA256(XXXX+%s) = %s' % (''.join(m[4:]), check))
    if sha256((input().strip() + ''.join(m[4:])).encode()).hexdigest() == check:
        return True
    return False


def check_power(secret):
    alarm(100)
    for i in range(1, 5):
        x = getRandomNBitInteger(30 * i)
        my_ecc_table = {
            'p': hex(get_prime(128))[2:].rjust(64, '0'),
            'a': '0' * 63 + '1',
            'b': '0' * 64,
        }
        E = Curve(private_key=x, ecc_table=my_ecc_table)
        g = E.ecc_table['g']
        xg = E.kg(x, g)
        print('p:', E.ecc_table['p'])
        print('g:', g)
        print('xg:', xg)
        if int(input('your guess:').strip()) == x:
            print('here your gift:', secret[17 * i * (i - 1):17 * i * (i + 1)])
        else:
            break
    return i


def use_power(power):
    my_ecc_table = {
        'n': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123',
        'p': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF',
        'g': '32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7'
             'bc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0',
        'a': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC',
        'b': '28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93',
    }
    E = Curve(ecc_table=my_ecc_table)
    for i in range(60):
        print('signature %d:' % (i+1), E.sign(sha256(b'Lord Riot %d' % i).digest(), getRandomRange(2 ** (255 - power), 2 ** (256 - power))))
    alarm(10)
    if int(input('your guess:', ).strip()) == E.private_key:
        print('here your flag:', flag2)


def main():
    if not proof_of_work():
        return
    secret = bin(bytes_to_long(flag1))[2:].ljust(340, '0')
    use_power(check_power(secret))


if __name__ == '__main__':
    main()
