#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>

using namespace std;

// ESTRUCTURAS DE DATOS


struct GeneInfo {
    string name;
    int start_index;
    int end_index;
    string sequence;
    string first_12_chars;
};

struct PalindromeInfo {
    string sequence;
    int start_index;
    int length;
};

struct ProteinInfo {
    string name;
    string sequence;
    int genome_start;
    int genome_end;
    string first_4_amino_acids;
    string first_4_codons;
};

struct GenomeDifference {
    int position;
    char wuhan_nucleotide;
    char texas_nucleotide;
    string wuhan_codon;
    string texas_codon;
    string wuhan_amino_acid;
    string texas_amino_acid;
};


// CLASE BASE PARA UTILIDADES COMUNES

class StringUtils {
public:
    // Comparación de strings sin usar str.find
    static bool stringCompare(const string& str1, const string& str2) {
        if (str1.length() != str2.length()) return false;
        for (size_t i = 0; i < str1.length(); i++) {
            if (str1[i] != str2[i]) return false;
        }
        return true;
    }

    // Búsqueda de substring sin usar str.find
    static int findSubstring(const string& text, const string& pattern, int start_pos = 0) {
        if (pattern.length() > text.length()) return -1;

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

    // Buscar todas las apariciones de un substring
    static vector<int> findAllSubstrings(const string& text, const string& pattern) {
        vector<int> positions;
        if (pattern.length() > text.length()) return positions;

        int start_pos = 0;
        while (start_pos <= (int)text.length() - (int)pattern.length()) {
            int pos = findSubstring(text, pattern, start_pos);
            if (pos == -1) break;

            positions.push_back(pos);
            start_pos = pos + 1; // Continuar desde la siguiente posición
        }

        return positions;
    }

    // Verificar si una cadena es palíndromo
    static bool isPalindrome(const string& str) {
        int left = 0;
        int right = str.length() - 1;
        while (left < right) {
            if (str[left] != str[right]) return false;
            left++;
            right--;
        }
        return true;
    }

    // Encontrar la subcadena común más larga usando programación dinámica
    static string findLongestCommonSubstring(const string& str1, const string& str2) {
        int m = str1.length();
        int n = str2.length();

        vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
        int max_length = 0;
        int end_pos = 0;

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

        if (max_length > 0) {
            return str1.substr(end_pos - max_length + 1, max_length);
        }
        return "";
    }
};


// CLASE PARA MANEJO DE CODONES Y TRADUCCIÓN


class CodonTranslator {
private:
    map<string, string> codon_table;

public:
    CodonTranslator() {
        initializeCodonTable();
    }

    void initializeCodonTable() {
        codon_table["TTT"] = "F"; codon_table["TTC"] = "F";
        codon_table["TTA"] = "L"; codon_table["TTG"] = "L";
        codon_table["TCT"] = "S"; codon_table["TCC"] = "S"; codon_table["TCA"] = "S"; codon_table["TCG"] = "S";
        codon_table["TAT"] = "Y"; codon_table["TAC"] = "Y";
        codon_table["TGT"] = "C"; codon_table["TGC"] = "C";
        codon_table["TGG"] = "W";
        codon_table["CTT"] = "L"; codon_table["CTC"] = "L"; codon_table["CTA"] = "L"; codon_table["CTG"] = "L";
        codon_table["CCT"] = "P"; codon_table["CCC"] = "P"; codon_table["CCA"] = "P"; codon_table["CCG"] = "P";
        codon_table["CAT"] = "H"; codon_table["CAC"] = "H";
        codon_table["CAA"] = "Q"; codon_table["CAG"] = "Q";
        codon_table["CGT"] = "R"; codon_table["CGC"] = "R"; codon_table["CGA"] = "R"; codon_table["CGG"] = "R";
        codon_table["ATT"] = "I"; codon_table["ATC"] = "I"; codon_table["ATA"] = "I";
        codon_table["ATG"] = "M"; // START
        codon_table["ACT"] = "T"; codon_table["ACC"] = "T"; codon_table["ACA"] = "T"; codon_table["ACG"] = "T";
        codon_table["AAT"] = "N"; codon_table["AAC"] = "N";
        codon_table["AAA"] = "K"; codon_table["AAG"] = "K";
        codon_table["AGT"] = "S"; codon_table["AGC"] = "S";
        codon_table["AGA"] = "R"; codon_table["AGG"] = "R";
        codon_table["GTT"] = "V"; codon_table["GTC"] = "V"; codon_table["GTA"] = "V"; codon_table["GTG"] = "V";
        codon_table["GCT"] = "A"; codon_table["GCC"] = "A"; codon_table["GCA"] = "A"; codon_table["GCG"] = "A";
        codon_table["GAT"] = "D"; codon_table["GAC"] = "D";
        codon_table["GAA"] = "E"; codon_table["GAG"] = "E";
        codon_table["GGT"] = "G"; codon_table["GGC"] = "G"; codon_table["GGA"] = "G"; codon_table["GGG"] = "G";
        codon_table["TAA"] = "STOP"; codon_table["TAG"] = "STOP"; codon_table["TGA"] = "STOP";
    }

    string codonToAminoAcid(const string& codon) {
        if (codon_table.find(codon) != codon_table.end()) {
            return codon_table[codon];
        }
        return "?";
    }

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
};


// CLASE PARA MANEJO DE ARCHIVOS


class FileManager {
public:
    static string loadGenomeFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: No se pudo abrir el archivo " << filename << endl;
            return "";
        }

        string genome = "";
        string line;
        getline(file, line); // Saltar la primera línea (header)
        while (getline(file, line)) {
            genome += line;
        }
        file.close();

        return genome;
    }

    static string loadGeneFile(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: No se pudo abrir el archivo " << filename << endl;
            return "";
        }

        string gene = "";
        string line;
        getline(file, line); // Saltar la primera línea (header)
        while (getline(file, line)) {
            gene += line;
        }
        file.close();

        return gene;
    }

    static pair<int, int> extractGenePositions(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: No se pudo abrir el archivo " << filename << endl;
            return {-1, -1};
        }

        string header;
        getline(file, header);
        file.close();

        // Buscar el patrón :start-end en el header
        size_t colon_pos = header.find(':');
        if (colon_pos == string::npos) return {-1, -1};

        size_t dash_pos = header.find('-', colon_pos);
        if (dash_pos == string::npos) return {-1, -1};

        int start = stoi(header.substr(colon_pos + 1, dash_pos - colon_pos - 1));
        int end = stoi(header.substr(dash_pos + 1));

        return {start, end};
    }

    static vector<ProteinInfo> loadProteinsFile(const string& filename) {
        vector<ProteinInfo> proteins;
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: No se pudo abrir el archivo " << filename << endl;
            return proteins;
        }

        string current_protein_name = "";
        string current_protein_seq = "";
        string line;

        while (getline(file, line)) {
            if (line[0] == '>') {
                if (!current_protein_name.empty()) {
                    ProteinInfo protein;
                    protein.name = current_protein_name;
                    protein.sequence = current_protein_seq;
                    proteins.push_back(protein);
                }
                current_protein_name = line.substr(1);
                current_protein_seq = "";
            } else {
                current_protein_seq += line;
            }
        }

        // Agregar la última proteína
        if (!current_protein_name.empty()) {
            ProteinInfo protein;
            protein.name = current_protein_name;
            protein.sequence = current_protein_seq;
            proteins.push_back(protein);
        }

        file.close();
        return proteins;
    }

    static void savePalindromesToFile(const vector<PalindromeInfo>& palindromes, const vector<string>& gene_names) {
        ofstream output("palindromes_results.txt");

        output << "=== PALÍNDROMOS MÁS LARGOS EN CADA GEN ===" << endl;

        for (size_t i = 0; i < palindromes.size(); i++) {
            output << "\nGen " << gene_names[i] << ":" << endl;
            output << "  Longitud del palíndromo más largo: " << palindromes[i].length << endl;
            output << "  Posición de inicio: " << palindromes[i].start_index << endl;
            output << "  Secuencia: " << palindromes[i].sequence << endl;
        }

        output.close();
    }
};


// TAREA 1: BÚSQUEDA DE GENES


class GeneFinder {
private:
    string wuhan_genome;
    string gene_M;
    string gene_S;
    string gene_ORF1AB;

public:
    GeneFinder(const string& genome, const string& m, const string& s, const string& orf1ab) 
        : wuhan_genome(genome), gene_M(m), gene_S(s), gene_ORF1AB(orf1ab) {}

    vector<GeneInfo> findGenes() {
        vector<GeneInfo> genes;

        // Extraer posiciones de los genes desde los archivos
        auto m_pos = FileManager::extractGenePositions("files/gen-M.txt");
        auto s_pos = FileManager::extractGenePositions("files/gen-S.txt");
        auto orf1ab_pos = FileManager::extractGenePositions("files/gen-ORF1AB.txt");

        // Usar las posiciones conocidas para extraer los genes del genoma
        if (m_pos.first != -1 && m_pos.second != -1) {
            // Extraer la secuencia del genoma en la posición conocida
            string genome_m = wuhan_genome.substr(m_pos.first - 1, m_pos.second - m_pos.first + 1);

            // Buscar todas las apariciones de esta secuencia en el genoma
            vector<int> m_positions = StringUtils::findAllSubstrings(wuhan_genome, genome_m);
            for (int pos : m_positions) {
                GeneInfo gene = createGeneInfo("M", pos, genome_m);
                genes.push_back(gene);
            }
        }

        if (s_pos.first != -1 && s_pos.second != -1) {
            // Extraer la secuencia del genoma en la posición conocida
            string genome_s = wuhan_genome.substr(s_pos.first - 1, s_pos.second - s_pos.first + 1);

            // Buscar todas las apariciones de esta secuencia en el genoma
            vector<int> s_positions = StringUtils::findAllSubstrings(wuhan_genome, genome_s);
            for (int pos : s_positions) {
                GeneInfo gene = createGeneInfo("S", pos, genome_s);
                genes.push_back(gene);
            }
        }

        if (orf1ab_pos.first != -1 && orf1ab_pos.second != -1) {
            // Extraer la secuencia del genoma en la posición conocida
            string genome_orf1ab = wuhan_genome.substr(orf1ab_pos.first - 1, orf1ab_pos.second - orf1ab_pos.first + 1);

            // Buscar todas las apariciones de esta secuencia en el genoma
            vector<int> orf1ab_positions = StringUtils::findAllSubstrings(wuhan_genome, genome_orf1ab);
            for (int pos : orf1ab_positions) {
                GeneInfo gene = createGeneInfo("ORF1AB", pos, genome_orf1ab);
                genes.push_back(gene);
            }
        }

        return genes;
    }

    void displayResults(const vector<GeneInfo>& genes) {
        cout << "\n=== BÚSQUEDA DE GENES ===" << endl;

        for (const auto& gene : genes) {
            cout << "\nGen " << gene.name << ":" << endl;
            cout << "  Índices: " << gene.start_index << " - " << gene.end_index << endl;
            cout << "  Primeros 12 caracteres: " << gene.first_12_chars << endl;
        }
    }

private:
    GeneInfo createGeneInfo(const string& name, int start, const string& sequence) {
        GeneInfo gene;
        gene.name = name;
        gene.start_index = start;
        gene.end_index = start + sequence.length() - 1;
        gene.sequence = sequence;
        gene.first_12_chars = sequence.substr(0, min(12, (int)sequence.length()));
        return gene;
    }
};


// TAREA 2: ANÁLISIS DE PALÍNDROMOS


class PalindromeAnalyzer {
public:
    static PalindromeInfo findLongestPalindrome(const string& sequence) {
        PalindromeInfo result;
        result.length = 0;
        result.start_index = -1;
        result.sequence = "";

        // Verificar todos los posibles palíndromos
        for (size_t i = 0; i < sequence.length(); i++) {
            for (size_t j = i + 1; j <= sequence.length(); j++) {
                string substr = sequence.substr(i, j - i);
                if (StringUtils::isPalindrome(substr) && (int)substr.length() > result.length) {
                    result.length = substr.length();
                    result.start_index = i;
                    result.sequence = substr;
                }
            }
        }
        return result;
    }

    static vector<PalindromeInfo> analyzeGenes(const vector<GeneInfo>& genes) {
        vector<PalindromeInfo> palindromes;

        for (const auto& gene : genes) {
            PalindromeInfo palindrome = findLongestPalindrome(gene.sequence);
            palindromes.push_back(palindrome);
        }

        return palindromes;
    }

    static void displayResults(const vector<PalindromeInfo>& palindromes, const vector<GeneInfo>& genes) {
        cout << "\n=== PALÍNDROMOS MÁS LARGOS EN CADA GEN ===" << endl;

        for (size_t i = 0; i < palindromes.size(); i++) {
            cout << "\nGen " << genes[i].name << ":" << endl;
            cout << "  Longitud del palíndromo más largo: " << palindromes[i].length << endl;
            cout << "  Posición de inicio: " << palindromes[i].start_index << endl;
            cout << "  Secuencia: " << palindromes[i].sequence << endl;
        }
    }

    static void saveResults(const vector<PalindromeInfo>& palindromes, const vector<GeneInfo>& genes) {
        vector<string> gene_names;
        for (const auto& gene : genes) {
            gene_names.push_back(gene.name);
        }
        FileManager::savePalindromesToFile(palindromes, gene_names);
    }
};


// TAREA 3: MAPEO DE PROTEÍNAS


class ProteinMapper {
private:
    string wuhan_genome;
    CodonTranslator translator;

public:
    ProteinMapper(const string& genome) : wuhan_genome(genome) {}

    void mapProteinsToGenome(vector<ProteinInfo>& proteins) {
        cout << "\n=== MAPEO DE PROTEÍNAS A SECCIONES DEL GENOMA ===" << endl;

        for (auto& protein : proteins) {
            mapSingleProtein(protein);
        }
    }

private:
    void mapSingleProtein(ProteinInfo& protein) {
        string protein_sequence = protein.sequence;

        // Buscar en diferentes marcos de lectura
        for (int frame = 0; frame < 3; frame++) {
            for (size_t i = frame; i < wuhan_genome.length() - 2; i += 3) {
                string codon = wuhan_genome.substr(i, 3);
                string amino_acid = translator.codonToAminoAcid(codon);

                if (amino_acid == protein_sequence.substr(0, 1)) {
                    if (checkProteinMatch(protein_sequence, i)) {
                        setProteinInfo(protein, i, protein_sequence);
                        displayProteinInfo(protein);
                        return;
                    }
                }
            }
        }
    }

    bool checkProteinMatch(const string& protein_sequence, size_t start_pos) {
        for (size_t j = 0; j < min(protein_sequence.length(), (size_t)10); j++) {
            if (start_pos + j * 3 + 2 >= wuhan_genome.length()) {
                return false;
            }
            string test_codon = wuhan_genome.substr(start_pos + j * 3, 3);
            string test_aa = translator.codonToAminoAcid(test_codon);
            if (test_aa != protein_sequence.substr(j, 1)) {
                return false;
            }
        }
        return true;
    }

    void setProteinInfo(ProteinInfo& protein, size_t start_pos, const string& protein_sequence) {
        protein.genome_start = start_pos;
        protein.genome_end = start_pos + protein_sequence.length() * 3 - 1;

        // Obtener primeros 4 aminoácidos y codones
        protein.first_4_amino_acids = protein_sequence.substr(0, min(4, (int)protein_sequence.length()));
        protein.first_4_codons = "";
        for (int k = 0; k < min(4, (int)protein_sequence.length()); k++) {
            protein.first_4_codons += wuhan_genome.substr(start_pos + k * 3, 3) + " ";
        }
    }

    void displayProteinInfo(const ProteinInfo& protein) {
        cout << "\nProteína: " << protein.name << endl;
        cout << "  Índices en el genoma: " << protein.genome_start << " - " << protein.genome_end << endl;
        cout << "  Primeros 4 aminoácidos: " << protein.first_4_amino_acids << endl;
        cout << "  Codones asociados: " << protein.first_4_codons << endl;
    }
};


// TAREA 4: COMPARACIÓN DE GENOMAS


class GenomeComparator {
private:
    string wuhan_genome;
    string texas_genome;
    CodonTranslator translator;

public:
    GenomeComparator(const string& wuhan, const string& texas) 
        : wuhan_genome(wuhan), texas_genome(texas) {}

    vector<GenomeDifference> compareGenomes() {
        cout << "\n=== COMPARACIÓN DE GENOMAS WUHAN vs TEXAS ===" << endl;

        vector<GenomeDifference> differences;
        size_t min_length = min(wuhan_genome.length(), texas_genome.length());

        // Comparar nucleótido por nucleótido
        for (size_t i = 0; i < min_length; i++) {
            if (wuhan_genome[i] != texas_genome[i]) {
                GenomeDifference diff = createDifference(i);
                differences.push_back(diff);
            }
        }

        displayResults(differences);
        return differences;
    }

private:
    GenomeDifference createDifference(size_t position) {
        GenomeDifference diff;
        diff.position = position;
        diff.wuhan_nucleotide = wuhan_genome[position];
        diff.texas_nucleotide = texas_genome[position];

        // Obtener codones - calcular el codón completo que contiene la diferencia
        int codon_start = (position / 3) * 3;
        if (codon_start + 2 < (int)min(wuhan_genome.length(), texas_genome.length())) {
            // Obtener el codón completo de Wuhan
            diff.wuhan_codon = wuhan_genome.substr(codon_start, 3);

            // Obtener el codón completo de Texas
            diff.texas_codon = texas_genome.substr(codon_start, 3);

            // Traducir ambos codones a aminoácidos
            diff.wuhan_amino_acid = translator.codonToAminoAcid(diff.wuhan_codon);
            diff.texas_amino_acid = translator.codonToAminoAcid(diff.texas_codon);
        }

        return diff;
    }

    void displayResults(const vector<GenomeDifference>& differences) {
        cout << "\nDiferencias encontradas: " << differences.size() << endl;

        for (const auto& diff : differences) {
            cout << "\nPosición: " << diff.position << endl;
            cout << "  Nucleótido Wuhan: " << diff.wuhan_nucleotide << endl;
            cout << "  Nucleótido Texas: " << diff.texas_nucleotide << endl;

            // Calcular la posición dentro del codón (0, 1, o 2)
            int position_in_codon = diff.position % 3;
            cout << "  Posición en el codón: " << (position_in_codon + 1) << "ª posición" << endl;

            cout << "  Codón Wuhan: " << diff.wuhan_codon << endl;
            cout << "  Codón Texas: " << diff.texas_codon << endl;
            cout << "  Aminoácido Wuhan: " << diff.wuhan_amino_acid << endl;
            cout << "  Aminoácido Texas: " << diff.texas_amino_acid << endl;

            if (diff.wuhan_amino_acid != diff.texas_amino_acid) {
                cout << "  *** CAMBIO DE AMINOÁCIDO ***" << endl;
            } else {
                cout << "  (Sin cambio de aminoácido - mutación silenciosa)" << endl;
            }
        }
    }
};


// CLASE PRINCIPAL COORDINADORA


class SARSCOV2Analyzer {
private:
    string wuhan_genome;
    string texas_genome;
    string gene_M;
    string gene_S;
    string gene_ORF1AB;
    vector<ProteinInfo> proteins;

public:
    SARSCOV2Analyzer() {}

    bool loadAllData() {
        cout << "=== CARGANDO DATOS ===" << endl;

        // Cargar genomas
        wuhan_genome = FileManager::loadGenomeFile("files/SARS-COV-2-MN908947.3.txt");
        texas_genome = FileManager::loadGenomeFile("files/SARS-COV-2-MT106054.1.txt");

        if (wuhan_genome.empty() || texas_genome.empty()) {
            return false;
        }

        // Cargar genes
        gene_M = FileManager::loadGeneFile("files/gen-M.txt");
        gene_S = FileManager::loadGeneFile("files/gen-S.txt");
        gene_ORF1AB = FileManager::loadGeneFile("files/gen-ORF1AB.txt");

        // Cargar proteínas
        proteins = FileManager::loadProteinsFile("files/seq-proteins.txt");

        cout << "✓ Archivos cargados exitosamente" << endl;
        cout << "  Longitud del genoma Wuhan: " << wuhan_genome.length() << " nucleótidos" << endl;
        cout << "  Longitud del genoma Texas: " << texas_genome.length() << " nucleótidos" << endl;
        cout << "  Proteínas cargadas: " << proteins.size() << endl;

        return true;
    }

    void runAnalysis() {
        cout << "\n=== ANÁLISIS COMPLETO ===" << endl;

        if (!loadAllData()) {
            cerr << "Error al cargar los archivos" << endl;
            return;
        }

        // TAREA 1: Búsqueda de genes
        GeneFinder geneFinder(wuhan_genome, gene_M, gene_S, gene_ORF1AB);
        vector<GeneInfo> genes = geneFinder.findGenes();
        geneFinder.displayResults(genes);

        // TAREA 2: Análisis de palíndromos
        vector<PalindromeInfo> palindromes = PalindromeAnalyzer::analyzeGenes(genes);
        PalindromeAnalyzer::displayResults(palindromes, genes);
        PalindromeAnalyzer::saveResults(palindromes, genes);

        // TAREA 3: Mapeo de proteínas
        ProteinMapper proteinMapper(wuhan_genome);
        proteinMapper.mapProteinsToGenome(proteins);

        // TAREA 4: Comparación de genomas
        GenomeComparator genomeComparator(wuhan_genome, texas_genome);
        genomeComparator.compareGenomes();

    }
};


// FUNCIÓN PRINCIPAL

int main() {
    SARSCOV2Analyzer analyzer;
    analyzer.runAnalysis();
    return 0;
}
