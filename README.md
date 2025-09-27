# Análisis del Genoma SARS-CoV-2

## Descripción del Proyecto

Este proyecto implementa un analizador completo del genoma del virus SARS-CoV-2 en C++, que realiza múltiples tareas de análisis genómico sin utilizar funciones de búsqueda de strings como `str.find()` o similares. El código está organizado de forma modular con clases especializadas para cada tarea.

## Estructura del Proyecto

```
Evidencia1/
├── files/
│   ├── SARS-COV-2-MN908947.3.txt    # Genoma de Wuhan 2019
│   ├── SARS-COV-2-MT106054.1.txt    # Genoma de Texas 2020
│   ├── gen-M.txt                     # Secuencia del gen M
│   ├── gen-S.txt                     # Secuencia del gen S
│   ├── gen-ORF1AB.txt               # Secuencia del gen ORF1AB
│   └── seq-proteins.txt             # Secuencias de proteínas
├── sars_cov2_analyzer.cpp           # Programa principal (versión modular)
├── Makefile                         # Para compilación
├── test_program.sh                  # Script de prueba
├── README.md                        # Documentación completa
├── ESTRUCTURA_MODULAR.md            # Documentación de la estructura modular
├── ARQUITECTURA.md                  # Diagramas y arquitectura
└── INSTRUCCIONES.md                 # Guía de uso
```

## Lógica de Implementación

### 1. Funciones de Comparación de Strings (Sin usar str.find)

**Problema:** Implementar comparación de strings sin usar funciones de búsqueda estándar.

**Solución:**
```cpp
bool stringCompare(const string& str1, const string& str2) {
    if (str1.length() != str2.length()) return false;
    for (size_t i = 0; i < str1.length(); i++) {
        if (str1[i] != str2[i]) return false;
    }
    return true;
}
```

**Lógica:** Comparación carácter por carácter usando un bucle for manual, evitando cualquier función de búsqueda predefinida.

### 2. Búsqueda de Substrings (Sin usar str.find)

**Problema:** Encontrar la posición de una subcadena en una cadena principal.

**Solución:**
```cpp
int findSubstring(const string& text, const string& pattern, int start_pos = 0) {
    for (size_t i = start_pos; i <= text.length() - pattern.length(); i++) {
        bool found = true;
        for (size_t j = 0; j < pattern.length(); j++) {
            if (text[i + j] != pattern[j]) {
                found = false;
                break;
            }
        }
        if (found) return i;
    }
    return -1;
}
```

**Lógica:** Algoritmo de fuerza bruta que compara cada posición posible de la subcadena con el patrón buscado.

### 3. Búsqueda de Palíndromos

**Problema:** Encontrar el palíndromo más largo en una secuencia genética.

**Solución:**
```cpp
bool isPalindrome(const string& str) {
    int left = 0;
    int right = str.length() - 1;
    while (left < right) {
        if (str[left] != str[right]) return false;
        left++;
        right--;
    }
    return true;
}

PalindromeInfo findLongestPalindrome(const string& sequence) {
    // Verificar todos los posibles palíndromos
    for (size_t i = 0; i < sequence.length(); i++) {
        for (size_t j = i + 1; j <= sequence.length(); j++) {
            string substr = sequence.substr(i, j - i);
            if (isPalindrome(substr) && substr.length() > result.length) {
                // Actualizar resultado si es más largo
            }
        }
    }
}
```

**Lógica:** 
- Verificación de palíndromo usando dos punteros (izquierda y derecha)
- Búsqueda exhaustiva de todas las subcadenas posibles
- Los palíndromos son importantes porque son regiones propensas a mutaciones

### 4. Subcadena Común Más Larga (Programación Dinámica)

**Problema:** Encontrar la subcadena común más larga entre dos secuencias.

**Solución:**
```cpp
string findLongestCommonSubstring(const string& str1, const string& str2) {
    int m = str1.length();
    int n = str2.length();
    
    // Crear tabla DP
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    int max_length = 0;
    int end_pos = 0;
    
    // Llenar la tabla DP
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (str1[i-1] == str2[j-1]) {
                dp[i][j] = dp[i-1][j-1] + 1;
                if (dp[i][j] > max_length) {
                    max_length = dp[i][j];
                    end_pos = i - 1;
                }
            } else {
                dp[i][j] = 0;
            }
        }
    }
    
    return str1.substr(end_pos - max_length + 1, max_length);
}
```

**Lógica:** 
- Tabla de programación dinámica 2D
- `dp[i][j]` representa la longitud de la subcadena común más larga que termina en posición `i` de str1 y `j` de str2
- Si los caracteres coinciden, incrementar la longitud; si no, reiniciar a 0

### 5. Marco de Lectura y Agrupación de Codones

**Conceptos Clave:**
- **UTR-5' y UTR-3':** Regiones no traducidas en los extremos 5' y 3' del ARNm
- **Marco de lectura:** Determina cómo se agrupan los nucleótidos en codones
- **Codón de inicio:** ATG (Metionina)
- **Codones de parada:** TAA, TAG, TGA

**Implementación:**
```cpp
string translateDNAToProtein(const string& dna_sequence) {
    string protein = "";
    for (size_t i = 0; i < dna_sequence.length() - 2; i += 3) {
        string codon = dna_sequence.substr(i, 3);
        string amino_acid = codonToAminoAcid(codon);
        if (amino_acid != "STOP") {
            protein += amino_acid;
        } else {
            break;
        }
    }
    return protein;
}
```

**Lógica:** 
- Agrupación en tripletas (codones) comenzando desde diferentes marcos (0, 1, 2)
- Traducción usando tabla de codones estándar
- Parada en codones de terminación

### 6. Tabla de Codones a Aminoácidos

**Implementación:**
```cpp
void initializeCodonTable() {
    codon_table["TTT"] = "F"; codon_table["TTC"] = "F";
    codon_table["TTA"] = "L"; codon_table["TTG"] = "L";
    // ... (tabla completa)
    codon_table["ATG"] = "M"; // START
    codon_table["TAA"] = "STOP"; codon_table["TAG"] = "STOP"; codon_table["TGA"] = "STOP";
}
```

**Lógica:** Mapeo directo de cada codón de 3 nucleótidos a su aminoácido correspondiente según el código genético estándar.

## Tareas Implementadas

### 1. Búsqueda de Genes (M, S, ORF1AB)

**Objetivo:** Encontrar los índices de aparición de cada gen en el genoma principal.

**Implementación:**
```cpp
vector<GeneInfo> findGenes() {
    vector<GeneInfo> genes;
    
    // Buscar cada gen usando findSubstring
    int m_start = findSubstring(wuhan_genome, gene_M);
    int s_start = findSubstring(wuhan_genome, gene_S);
    int orf1ab_start = findSubstring(wuhan_genome, gene_ORF1AB);
    
    // Almacenar información de cada gen encontrado
}
```

**Resultado:** Para cada gen se muestra:
- Nombre del gen
- Índices de inicio y fin
- Primeros 12 caracteres

### 2. Palíndromos en Genes

**Objetivo:** Encontrar el palíndromo más largo en cada gen.

**Implementación:** Usa la función `findLongestPalindrome()` para cada secuencia genética.

**Resultado:** 
- Longitud del palíndromo más largo
- Posición de inicio
- Secuencia del palíndromo
- Guardado en archivo `palindromes_results.txt`

### 3. Mapeo de Proteínas

**Objetivo:** Encontrar las secciones del genoma que codifican cada proteína.

**Implementación:**
```cpp
void mapProteinsToGenome() {
    for (auto& protein : proteins) {
        // Buscar en diferentes marcos de lectura
        for (int frame = 0; frame < 3; frame++) {
            for (size_t i = frame; i < wuhan_genome.length() - 2; i += 3) {
                // Verificar coincidencia de aminoácidos
                // Traducir codones y comparar
            }
        }
    }
}
```

**Resultado:** Para cada proteína:
- Índices en el genoma
- Primeros 4 aminoácidos
- Codones asociados

### 4. Comparación de Genomas

**Objetivo:** Comparar genomas de Wuhan 2019 vs Texas 2020.

**Implementación:**
```cpp
void compareGenomes() {
    for (size_t i = 0; i < min_length; i++) {
        if (wuhan_genome[i] != texas_genome[i]) {
            // Registrar diferencia
            // Obtener codones afectados
            // Traducir a aminoácidos
            // Verificar si hay cambio de aminoácido
        }
    }
}
```

**Resultado:**
- Posición de cada diferencia
- Nucleótidos diferentes
- Codones afectados
- Aminoácidos resultantes
- Identificación de cambios que afectan la proteína

## Compilación y Ejecución

```bash
# Compilar
make

# Ejecutar
make run

# O compilación manual
g++ -std=c++11 -Wall -Wextra -O2 -o sars_cov2_analyzer sars_cov2_analyzer.cpp

# Ejecutar
./sars_cov2_analyzer
```

## Estructuras de Datos Utilizadas

### GeneInfo
```cpp
struct GeneInfo {
    string name;
    int start_index;
    int end_index;
    string sequence;
    string first_12_chars;
};
```

### PalindromeInfo
```cpp
struct PalindromeInfo {
    string sequence;
    int start_index;
    int length;
};
```

### ProteinInfo
```cpp
struct ProteinInfo {
    string name;
    string sequence;
    int genome_start;
    int genome_end;
    string first_4_amino_acids;
    string first_4_codons;
};
```

### GenomeDifference
```cpp
struct GenomeDifference {
    int position;
    char wuhan_nucleotide;
    char texas_nucleotide;
    string wuhan_codon;
    string texas_codon;
    string wuhan_amino_acid;
    string texas_amino_acid;
};
```

## Complejidad Computacional

- **Búsqueda de substrings:** O(n*m) donde n es la longitud del texto y m la del patrón
- **Búsqueda de palíndromos:** O(n³) para verificar todas las subcadenas
- **Subcadena común más larga:** O(n*m) usando programación dinámica
- **Comparación de genomas:** O(n) donde n es la longitud del genoma más corto

## Aplicaciones en el Problema

### Utilidad de la Subcadena Común Más Larga

En el contexto del análisis genómico, la subcadena común más larga es útil para:

1. **Identificar regiones conservadas** entre diferentes variantes del virus
2. **Detectar secuencias funcionales** que se mantienen a través de la evolución
3. **Encontrar patrones de mutación** al comparar genomas
4. **Validar la calidad** de las secuencias genómicas
5. **Identificar sitios de recombinación** viral

### Marco de Lectura en SARS-CoV-2

El SARS-CoV-2 tiene un genoma de ARN de cadena positiva que se traduce directamente. Los marcos de lectura son cruciales porque:

1. **Determinan la estructura de las proteínas** producidas
2. **Identifican genes superpuestos** (genes que se leen en diferentes marcos)
3. **Permiten la detección de mutaciones** que cambian el marco de lectura
4. **Facilitan la predicción de proteínas** a partir de la secuencia genómica

## Consideraciones Técnicas

1. **Manejo de archivos:** Lectura secuencial de archivos FASTA
2. **Memoria:** Uso eficiente de estructuras de datos para genomas grandes
3. **Precisión:** Comparaciones exactas de nucleótidos y aminoácidos
4. **Robustez:** Manejo de errores en la carga de archivos
5. **Escalabilidad:** Algoritmos que pueden manejar genomas completos

## Resultados Esperados

El programa genera:

1. **Análisis de genes:** Posiciones y secuencias de los genes M, S y ORF1AB
2. **Archivo de palíndromos:** `palindromes_results.txt` con los palíndromos más largos
3. **Mapeo de proteínas:** Correspondencia entre proteínas y regiones genómicas
4. **Comparación de genomas:** Diferencias entre Wuhan 2019 y Texas 2020
5. **Análisis de impacto:** Identificación de mutaciones que cambian aminoácidos

Este análisis proporciona una comprensión profunda de la estructura y variación del genoma SARS-CoV-2, utilizando técnicas de programación avanzadas y conceptos de biología molecular.
