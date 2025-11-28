"""
Decodificador Reed-Solomon RS(15,9) GF(16), polinomio x^4 + x + 1.
"""

import os, sys

# =====================================================================
# GF(16) = GF(2)[x] / (x^4 + x + 1)
# =====================================================================

class GF16:
    def __init__(self):
        self.exp = [0]*32
        self.log = [0]*16
        self.mul = [[0]*16 for _ in range(16)]
        self.div = [[0]*16 for _ in range(16)]
        self._build()

    def _build(self):
        self.exp[0] = 1
        for i in range(1, 15):
            prev = self.exp[i-1]
            nxt = (prev << 1) & 0x0F
            if prev & 0x08:
                nxt ^= 0x03
            self.exp[i] = nxt

        self.exp[15] = 1
        for i in range(16, 32):
            self.exp[i] = self.exp[i % 15]

        self.log[0] = -1
        for i in range(15):
            self.log[self.exp[i]] = i

        for a in range(16):
            for b in range(16):
                if a == 0 or b == 0:
                    self.mul[a][b] = 0
                else:
                    self.mul[a][b] = self.exp[(self.log[a] + self.log[b]) % 15]

        for a in range(16):
            for b in range(1, 16):
                if a == 0:
                    self.div[a][b] = 0
                else:
                    self.div[a][b] = self.exp[(self.log[a] - self.log[b] + 15) % 15]

    def add(self, a, b): return a ^ b
    def mul2(self, a, b): return self.mul[a][b]
    def div2(self, a, b):
        if b == 0:
            raise ValueError("División por cero en GF(16)")
        return self.div[a][b]

gf16 = GF16()

# =====================================================================
# Polinomios sobre GF(16)
# =====================================================================

def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def deg(p):
    return len(trim(p[:])) - 1

def poly_add(p, q):
    m = max(len(p), len(q))
    r = [0] * m
    for i in range(m):
        a = p[i] if i < len(p) else 0
        b = q[i] if i < len(q) else 0
        r[i] = gf16.add(a, b)
    return trim(r)

def poly_scale(p, c):
    return trim([gf16.mul2(a, c) for a in p])

def poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i+j] = gf16.add(r[i+j], gf16.mul2(a, b))
    return trim(r)

def poly_divmod(num, den):
    num = num[:]
    den = trim(den[:])

    q = [0] * max(1, (len(num) - len(den) + 1))

    while deg(num) >= deg(den) and num != [0]:
        d = deg(num) - deg(den)
        coef = gf16.div2(num[-1], den[-1])
        q[d] = coef
        num = poly_add(num, [0]*d + poly_scale(den, coef))

    return trim(q), trim(num)

def poly_eval(p, x):
    y = 0
    power = 1
    for c in p:
        y = gf16.add(y, gf16.mul2(c, power))
        power = gf16.mul2(power, x)
    return y

def poly_derivative(p):
    d = []
    for i in range(1, len(p)):
        d.append(p[i] if i % 2 == 1 else 0)
    return trim(d)

# =====================================================================
# Reed-Solomon RS(15,9) por Euclides
# =====================================================================

class ReedSolomonEuclides:
    def __init__(self, n=15, k=9):
        self.n = n
        self.k = k
        self.t = (n - k) // 2

    def syndromes(self, w):
        return [poly_eval(w, gf16.exp[i]) for i in range(1, 2*self.t + 1)]

    def euclides(self, S):
        r_prev = [0]*(2*self.t) + [1]
        r_curr = S[:] 

        t_prev = [0]
        t_curr = [1]

        while deg(r_curr) >= self.t:
            q, r_next = poly_divmod(r_prev, r_curr)
            t_next = poly_add(t_prev, poly_mul(q, t_curr))
            r_prev, r_curr = r_curr, r_next
            t_prev, t_curr = t_curr, t_next

        return trim(t_curr), trim(r_curr)

    def normalizar(self, L_raw, O_raw):
        c = L_raw[0]
        inv = gf16.div2(1, c)
        return poly_scale(L_raw, inv), poly_scale(O_raw, inv)

    def chien(self, Lambda):
        pos = []
        for i in range(self.n):
            x = gf16.exp[(15 - i) % 15]
            if poly_eval(Lambda, x) == 0:
                pos.append(i)
        return pos

    def forney(self, Omega, Lambda, positions):
        Lp = poly_derivative(Lambda)
        mags = []
        for p in positions:
            x_inv = gf16.exp[(15 - p) % 15]
            num = poly_eval(Omega, x_inv)
            den = poly_eval(Lp, x_inv)
            mags.append(gf16.div2(num, den))
        return mags

    def decodificar_palabra(self, w):
        S = self.syndromes(w)

        if all(s == 0 for s in S):
            return w, w[self.n-self.k:], [], []

        L_raw, O_raw = self.euclides(S)
        Lambda, Omega = self.normalizar(L_raw, O_raw)

        pos = self.chien(Lambda)

        if len(pos) > self.t:
            raise ValueError("Más de 3 errores — palabra irrecuperable.")

        mags = self.forney(Omega, Lambda, pos)

        w_corr = w[:]
        errores = []
        for p, m in zip(pos, mags):
            orig = w_corr[p]
            corr = gf16.add(orig, m)
            errores.append((p, orig, m, corr))
            w_corr[p] = corr

        return w_corr, w_corr[self.n-self.k:], pos, errores

# =====================================================================
# Lectura A2 en 1 línea o varias
# =====================================================================

def asciihex_to_symbol(b):
    if 48 <= b <= 57: return b - 48
    if 65 <= b <= 70: return b - 55
    raise ValueError("Símbolos hex inválidos.")

def leer_A2(path, n=15):
    raw = open(path, "rb").read()
    syms = [asciihex_to_symbol(b) for b in raw if b not in b"\r\n\t "]

    if len(syms) % n != 0:
        raise ValueError("El archivo no es múltiplo de 15 símbolos.")

    return [syms[i:i+n] for i in range(0, len(syms), n)]

# =====================================================================
# Reconstrucción de bytes
# =====================================================================

def combinar_nibbles(h, l):
    return ((h & 0xF) << 4) | (l & 0xF)

def reconstruir_bytes(info_blocks):
    nibbles = []
    for info in info_blocks:
        nibbles.extend(info[:9])

    out = []
    for i in range(0, len(nibbles)-1, 2):
        out.append(combinar_nibbles(nibbles[i], nibbles[i+1]))
    return bytes(out)

# =====================================================================
# MAIN DECODIFICACIÓN
# =====================================================================

def decodificar_archivo(A2_path, A1_out):
    rs = ReedSolomonEuclides(15, 9)
    palabras = leer_A2(A2_path, 15)

    infos = []

    for idx, w in enumerate(palabras):
        print(f"\n--- Palabra RS #{idx} ---")
        try:
            w_corr, info, pos, errores = rs.decodificar_palabra(w)
            if (len(errores) > 3):
                raise ValueError("Más de 3 errores en una palabra")
            if pos:
                print(f"Errores detectados: {pos}")
                for (p, o, m, c) in errores:
                    print(f"  pos {p}: {o:X} -> {c:X} (e={m:X})")
            else:
                print("Sin errores.")

            infos.append(info)

        except ValueError as e:
            print(">>> PALABRA IRRECUPERABLE:", e)
            return

    data = reconstruir_bytes(infos)
    with open(A1_out, "wb") as f:
        f.write(data)
    return data

def elegir_archivo_txt(carpeta):
    txts = [f for f in os.listdir(carpeta) if f.lower().endswith(".txt")]
    txts.sort()

    if not txts:
        print("No hay archivos .txt en la carpeta:", carpeta)
        return None

    print("\nArchivos .txt disponibles en:", carpeta)
    for i, name in enumerate(txts):
        print(f"  {i+1}) {name}")

    while True:
        op = input("Elegí un archivo por número (ENTER para cancelar)> ").strip()
        if op == "":
            return None
        if op.isdigit():
            k = int(op)
            if 1 <= k <= len(txts):
                return os.path.join(carpeta, txts[k-1])
        print("Opción inválida. Probá de nuevo.")

# =====================================================================
# MAIN
# =====================================================================

if __name__ == "__main__":

    # carpeta por defecto = donde está el script
    base = os.path.dirname(os.path.abspath(__file__))

    # si querés elegir otra carpeta, descomentá esto:
    # carpeta = input(f"Carpeta donde buscar .txt (ENTER = {base})> ").strip()
    # if carpeta:
    #     base = carpeta

    A2_path = elegir_archivo_txt(base)
    if A2_path is None:
        print("Cancelado.")
        sys.exit(0)

    out = os.path.join(base, "A1_decodificado.txt")

    try:
        data = decodificar_archivo(A2_path, out)
        if (data):
            print("\nASCII:", data.decode("ascii", errors="replace"))
            print("Salida guardada en:", out)
    except Exception as e:
        print("✗ Error general:", e)
        raise