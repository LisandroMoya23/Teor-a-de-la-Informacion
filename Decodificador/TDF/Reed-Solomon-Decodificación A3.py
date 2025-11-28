"""
Decodificador Reed-Solomon RS(15,9) para A3
GF(16), polinomio x^4 + x + 1.
"""

import os
import sys
import random

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

    def decodificar_palabra(self, w, verbose=False):
        S = self.syndromes(w)
        
        if verbose:
            print("  Síndromes:", [f"{s:X}" for s in S])

        if all(s == 0 for s in S):
            return w, w[self.n-self.k:], [], []

        L_raw, O_raw = self.euclides(S)
        Lambda, Omega = self.normalizar(L_raw, O_raw)
        pos = self.chien(Lambda)

        if len(pos) > self.t:
            raise ValueError(f"Más de {self.t} errores – palabra irrecuperable.")

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
# DESENTRELAZADO ORIGINAL (columna por columna de TODAS las palabras)
# =====================================================================

def desentrelazar_codigos_original(datos_entrelazados, num_palabras, n=15):
    """
    Invierte el entrelazado ORIGINAL que lee columna por columna.
    
    Formato original:
    - Lee columna 0 de palabra 0, 1, 2, ..., M-1
    - Luego columna 1 de palabra 0, 1, 2, ..., M-1
    - etc.
    
    Entrada: lista plana [col0_w0, col0_w1, ..., col0_wM, col1_w0, ...]
    Salida: matriz [num_palabras x n] con las palabras RS originales.
    """
    palabras = [[0]*n for _ in range(num_palabras)]
    idx = 0
    
    # Recorremos columnas primero, luego filas
    for j in range(n):              # columnas 0..14
        for i in range(num_palabras):  # palabras 0..M-1
            if idx < len(datos_entrelazados):
                palabras[i][j] = datos_entrelazados[idx]
                idx += 1
    
    return palabras


# =====================================================================
# Lectura de A3
# =====================================================================

def asciihex_to_symbol(b):
    if 48 <= b <= 57: return b - 48
    if 65 <= b <= 70: return b - 55
    raise ValueError("Símbolos hex inválidos.")

def leer_A3(path):
    """Lee A3 y devuelve lista plana de símbolos entrelazados"""
    raw = open(path, "rb").read()
    syms = [asciihex_to_symbol(b) for b in raw if b not in b"\r\n\t "]
    return syms

# =====================================================================
# Inserción de Ráfagas de Errores
# =====================================================================

def insertar_rafaga(simbolos, pos_inicio, longitud):
    """Inserta una ráfaga de errores consecutivos"""
    datos = simbolos[:]
    cambios = []
    
    for i in range(longitud):
        pos = pos_inicio + i
        if pos >= len(datos):
            break
        
        viejo = datos[pos]
        k = random.randint(1, 14)
        e = gf16.exp[k]
        nuevo = gf16.add(viejo, e)
        
        datos[pos] = nuevo
        cambios.append((pos, viejo, nuevo, e, k))
    
    return datos, cambios

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
# MAIN DECODIFICACIÓN A3
# =====================================================================

def decodificar_A3(A3_path, A1_out, insertar_errores=False, longitud_rafaga=0, pos_rafaga=0):
    """
    Decodifica A3 (entrelazado ORIGINAL) y reconstruye A1.
    Opcionalmente inserta ráfaga de errores para testear.
    """
    print("\n" + "="*70)
    print("DECODIFICADOR A3 (ENTRELAZADO ORIGINAL - COLUMNA POR COLUMNA)")
    print("="*70)
    
    # Leer A3
    print(f"\n[1] Leyendo A3: {os.path.basename(A3_path)}")
    simbolos_A3 = leer_A3(A3_path)
    num_palabras = len(simbolos_A3) // 15
    print(f"    Total símbolos: {len(simbolos_A3)}")
    print(f"    Palabras RS: {num_palabras}")
    
    # Insertar ráfaga si se solicita
    if insertar_errores:
        print(f"\n[2] Insertando ráfaga de {longitud_rafaga} errores en pos {pos_rafaga}...")
        simbolos_A3, cambios = insertar_rafaga(simbolos_A3, pos_rafaga, longitud_rafaga)
        print(f"    ✓ {len(cambios)} errores insertados")
        if cambios:
            print(f"    Ejemplo: pos {cambios[0][0]}: {cambios[0][1]:X} → {cambios[0][2]:X} (e=α^{cambios[0][4]})")
    
    # Desentrelazar (VERSIÓN ORIGINAL)
    print(f"\n[3] Desentrelazando códigos (método ORIGINAL)...")
    palabras = desentrelazar_codigos_original(simbolos_A3, num_palabras, n=15)
    print(f"    ✓ {len(palabras)} palabras reconstruidas")
    
    # Decodificar
    print(f"\n[4] Decodificando palabras RS...")
    rs = ReedSolomonEuclides(15, 9)
    
    infos = []
    total_errores_corregidos = 0
    palabras_con_errores = 0
    palabras_irrecuperables = []
    
    for idx, w in enumerate(palabras):
        try:
            w_corr, info, pos, errores = rs.decodificar_palabra(w, verbose=False)
            
            if errores:
                palabras_con_errores += 1
                total_errores_corregidos += len(errores)
                if idx < 5:  # Mostrar primeras 5 palabras con errores
                    print(f"    W{idx:02d}: {len(errores)} error(es) en pos {pos}")
            
            infos.append(info)

        except ValueError as e:
            print(f"    W{idx:02d}: ✗ IRRECUPERABLE")
            palabras_irrecuperables.append(idx)
            infos.append([0]*9)  # Padding para no romper estructura
    
    # Estadísticas
    print(f"\n[5] Estadísticas de decodificación:")
    print(f"    Total palabras:          {num_palabras}")
    print(f"    Palabras sin errores:    {num_palabras - palabras_con_errores - len(palabras_irrecuperables)}")
    print(f"    Palabras con errores (corregidas): {palabras_con_errores} ✓")
    print(f"    Palabras irrecuperables: {len(palabras_irrecuperables)} ✗")
    print(f"    Total errores corregidos: {total_errores_corregidos}")
    
    if palabras_irrecuperables:
        print(f"    Palabras perdidas: {palabras_irrecuperables[:20]}")
        if len(palabras_irrecuperables) > 20:
            print(f"                       ... y {len(palabras_irrecuperables)-20} más")
    
    # Reconstruir A1
    print(f"\n[6] Reconstruyendo A1...")
    data = reconstruir_bytes(infos)
    
    with open(A1_out, "wb") as f:
        f.write(data)
    
    print(f"    ✓ Guardado: {os.path.basename(A1_out)}")
    print(f"    Tamaño: {len(data)} bytes")
    
    # Mostrar contenido COMPLETO
    try:
        contenido = data.decode("ascii", errors="replace")
        print(f"\n[7] Contenido recuperado COMPLETO:")
        print(f"\n{'='*70}")
        print(contenido)
        print(f"{'='*70}\n")
    except:
        print(f"\n[7] Contenido binario (hex completo):")
        print(f"    {data.hex()}")
    
    print("\n" + "="*70)
    
    return data, palabras_irrecuperables

# =====================================================================
# Interfaz de usuario
# =====================================================================

def elegir_archivo_txt(carpeta):
    txts = [f for f in os.listdir(carpeta) if f.lower().endswith(".txt")]
    txts.sort()

    if not txts:
        print("No hay archivos .txt en la carpeta:", carpeta)
        return None

    print("\nArchivos .txt disponibles:")
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
    
    carpeta = os.path.dirname(os.path.abspath(__file__))
    
    print("="*70)
    print("DECODIFICADOR RS(15,9) PARA A3 (ENTRELAZADO ORIGINAL)")
    print("="*70)
    
    # Seleccionar archivo A3
    A3_path = elegir_archivo_txt(carpeta)
    
    if A3_path is None:
        print("Cancelado.")
        sys.exit(0)
    
    # Verificar que sea A3
    if "A3" not in os.path.basename(A3_path):
        print(f"\n⚠️ Advertencia: El archivo seleccionado no parece ser A3.")
        print(f"   Este decodificador espera archivos ENTRELAZADOS (A3).")
        continuar = input("   ¿Continuar de todos modos? (s/n)> ").lower()
        if not continuar.startswith('s'):
            sys.exit(0)
    
    # Preguntar si insertar errores
    print("\n" + "-"*70)
    test_errores = input("¿Insertar ráfaga de errores para testear? (s/n)> ").lower().startswith('s')
    
    insertar = False
    longitud = 0
    pos = 0
    
    if test_errores:
        simbolos_temp = leer_A3(A3_path)
        max_pos = len(simbolos_temp)
        
        print(f"\nArchivo tiene {max_pos} símbolos.")
        longitud = int(input("Longitud de la ráfaga (ej: 10, 20, 30, 45): ").strip())
        pos = int(input(f"Posición inicial (0-{max_pos-longitud}): ").strip())
        
        if pos + longitud > max_pos:
            print(f"⚠️ Ajustando: ráfaga hasta posición {max_pos}")
            longitud = max_pos - pos
        
        insertar = True
    
    # Archivo de salida
    if insertar:
        A1_out = os.path.join(carpeta, f"A1_desde_A3_rafaga{longitud}.txt")
    else:
        A1_out = os.path.join(carpeta, "A1_desde_A3.txt")
    
    # Decodificar
    try:
        data, irrec = decodificar_A3(A3_path, A1_out, insertar, longitud, pos)
        
        if len(irrec) == 0:
            print("\n✅ ÉXITO: Todas las palabras fueron recuperadas.")
        elif len(irrec) < 5:
            print(f"\n⚠️ PARCIAL: {len(irrec)} palabra(s) no pudieron recuperarse.")
        else:
            print(f"\n✗ ERROR: {len(irrec)} palabra(s) irrecuperables.")
        
    except Exception as e:
        print(f"\n✗ Error general: {e}")
        import traceback
        traceback.print_exc()