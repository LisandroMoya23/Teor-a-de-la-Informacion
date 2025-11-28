# Manual de Usuario - Sistema Reed-Solomon RS(15,9)

Este sistema implementa un codificador/decodificador de control de errores usando **Reed-Solomon RS(15,9)** sobre el campo de Galois **GF(16)**, con soporte para entrelazado de códigos y corrección de hasta **3 errores por palabra**.

Guía de Uso

### **PASO 1: Codificación (A1 → A2 + A3)**

### **PASO 2: Inserción de Errores (Opcional - Para Testing)**

#### ¿Qué hace?
Inserta errores controlados en un archivo A2 o A3 para probar la capacidad de corrección.

#### Cómo usar:
```
python Insertar-Errores.py
```

#### Proceso interactivo:
```
1. Selecciona el archivo (A2.txt o A3.txt)
2. Elige el modo:
   - [1] Manual: Especificas la magnitud del error (α^k)
   - [2] Random: Errores aleatorios
3. Indica cantidad de errores a insertar
4. ¿Todos en la misma palabra? (s/n)
   - s: Concentra errores en una palabra (máx 3)
   - n: Distribuye errores entre palabras
```

### **PASO 3: Decodificación A2 (Sin Entrelazado)**

#### ¿Qué hace?
Recupera el archivo original (A1) desde un archivo A2 con o sin errores.

#### Cómo usar:
```
python "Reed-Solomon-Decodificación A2.py"
```

#### Proceso interactivo:
```
1. Selecciona el archivo A2.txt a decodificar
2. El sistema decodifica cada palabra RS
3. Genera A1_decodificado.txt
```

### **PASO 4: Decodificación A3 (Con Entrelazado)**

#### ¿Qué hace?
Recupera el archivo original (A1) desde un archivo A3 entrelazado, con protección mejorada contra ráfagas de errores.

#### Cómo usar:
```
python "Reed-Solomon-Decodificación A3.py"
```

#### Proceso interactivo:
```
1. Selecciona el archivo A3.txt a decodificar
2. ¿Insertar ráfaga de errores para testear? (s/n)
   - Si eliges 's':
     - Especifica longitud de la ráfaga
     - Especifica posición inicial
3. El sistema:
   - Inserta errores (si se solicitó)
   - Desentrelaza los códigos
   - Decodifica cada palabra RS
   - Reconstruye A1
```