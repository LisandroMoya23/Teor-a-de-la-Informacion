import os, random

# ==========================================================
# GF(16) — x^4 + x + 1
# ==========================================================

class GF16:
    def __init__(self):
        self.exp = [0]*32
        self.log = [0]*16
        self._build()

    def _build(self):
        self.exp[0] = 1
        for i in range(1, 15):
            prev = self.exp[i-1]
            nxt = (prev << 1) & 0xF
            if prev & 0x8:
                nxt ^= 0x3   # reducción x^4 + x + 1
            self.exp[i] = nxt

        for i in range(15):
            self.log[self.exp[i]] = i

        # duplicamos ciclo
        for i in range(15, 32):
            self.exp[i] = self.exp[i-15]

    # suma = XOR
    def add(self, a, b): return a ^ b

gf = GF16()

# ==========================================================
# Función para elegir archivo .txt
# ==========================================================

def elegir_archivo_txt(carpeta):
    txts = [f for f in os.listdir(carpeta) if f.lower().endswith(".txt")]
    txts.sort()

    if not txts:
        print("No hay archivos .txt en la carpeta.")
        return None

    print("\nArchivos encontrados:")
    for i, name in enumerate(txts):
        print(f"  {i+1}) {name}")

    while True:
        op = input("\nElegí un archivo por número> ").strip()
        if op.isdigit():
            k = int(op)
            if 1 <= k <= len(txts):
                return os.path.join(carpeta, txts[k-1])
        print("Opción inválida.\n")

# ==========================================================
# Función para cargar y dividir en palabras de 15 símbolos
# ==========================================================

def cargar_palabras_rs(ruta, n=15):
    """
    Lee el archivo A2 y lo divide en palabras de n símbolos.
    Soporta formato en 1 línea o múltiples líneas.
    """
    with open(ruta, "r") as f:
        contenido = f.read()
    
    # Remover espacios, saltos de línea, etc.
    simbolos = "".join(contenido.split())
    
    # Verificar que sea múltiplo de n
    if len(simbolos) % n != 0:
        raise ValueError(f"El archivo no es múltiplo de {n} símbolos. Total: {len(simbolos)}")
    
    # Dividir en palabras de n símbolos
    palabras = [simbolos[i:i+n] for i in range(0, len(simbolos), n)]
    
    return palabras

# ==========================================================
# MAIN
# ==========================================================

if __name__ == "__main__":

    carpeta = os.path.dirname(os.path.abspath(__file__))
    ruta = elegir_archivo_txt(carpeta)

    if ruta is None:
        print("Cancelado.")
        exit()

    # Cargar y dividir en palabras de 15 símbolos
    bloques = cargar_palabras_rs(ruta, n=15)

    print(f"\nA2 cargado OK ({len(bloques)} palabras RS).\n")
    print("Bloques RS (15 símbolos c/u):")
    for i, w in enumerate(bloques):
        print(f"W{i:02d}: {w}")
    print()

    modo = input("Elegí modo:\n  1) Manual (posición + magnitud k)\n  2) Random\nmodo> ").strip()

    cant = int(input("Cantidad de errores a insertar> ").strip())
    mismo = input("¿Todos en la misma palabra? (s/n)> ").lower().startswith("s")

    print("\n==== ERRORES INSERTADOS ====")

    # Elegir palabra inicial si todos van en la misma
    if mismo:
        idx_fija = random.randint(0, len(bloques)-1)
        print(f"Todos los errores irán en la palabra W{idx_fija}\n")

    # Contador de errores por palabra
    errores_por_palabra = {i: set() for i in range(len(bloques))}

    for num_error in range(cant):

        # Elegir palabra
        if mismo:
            idx = idx_fija
        else:
            # Distribuir errores en diferentes palabras
            palabras_disponibles = [i for i in range(len(bloques)) 
                                   if len(errores_por_palabra[i]) < 3]
            
            if not palabras_disponibles:
                print(f"\n⚠️  Ya hay 3 errores en todas las palabras. No se pueden insertar más sin exceder t=3.")
                break
            
            idx = random.choice(palabras_disponibles)

        palabra = list(bloques[idx])

        # Elegir posición (que no tenga error ya)
        posiciones_disponibles = [p for p in range(15) 
                                 if p not in errores_por_palabra[idx]]
        
        if not posiciones_disponibles:
            print(f"\n⚠️  Palabra W{idx} ya tiene errores en todas las posiciones.")
            continue

        pos = random.choice(posiciones_disponibles)
        errores_por_palabra[idx].add(pos)

        viejo = int(palabra[pos], 16)

        # Elegir magnitud α^k
        if modo == "1":
            k = int(input(f"Error #{num_error+1} - k (1..14) para e = α^k: ").strip())
        else:
            k = random.randint(1, 14)

        e = gf.exp[k]

        nuevo = gf.add(viejo, e)
        palabra[pos] = format(nuevo, "X")
        bloques[idx] = "".join(palabra)

        print(f"W{idx:02d} pos {pos:02d}:  {viejo:X} + e=α^{k:02d}({e:X})  ->  {nuevo:X}")

    # Mostrar resumen
    print("\n==== RESUMEN ====")
    total_errores = 0
    for i in range(len(bloques)):
        if errores_por_palabra[i]:
            print(f"W{i:02d}: {len(errores_por_palabra[i])} error(es) en posiciones {sorted(errores_por_palabra[i])}")
            total_errores += len(errores_por_palabra[i])
    
    print(f"\nTotal: {total_errores} errores insertados en {len([e for e in errores_por_palabra.values() if e])} palabra(s)")

    # Guardar como *_err.txt (en múltiples líneas para mejor legibilidad)
    nuevo_path = ruta.replace(".txt", "_err.txt")
    with open(nuevo_path, "w") as f:
        for w in bloques:
            f.write(w + "\n")

    print(f"\nArchivo generado: {nuevo_path}")
    print("Renombralo a A2.txt para usarlo en el decodificador.\n")