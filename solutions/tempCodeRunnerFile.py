def count_nucleotides():
    """Count nucleotides in DNA string from rosalind_dna.txt"""
    path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\rosalind_dna.txt"
    with open(path, 'r') as file:
        text = file.read().strip()
    
    nucleotide = {'A': 0, "C": 0, 'G': 0, 'T': 0}
    for i in text:
        if i != '\n':
            nucleotide[i] += 1
    
    # Print counts in order A C G T
    result = " ".join(str(nucleotide[base]) for base in ['A', 'C', 'G', 'T'])
    print(result)
    return result

def dna_to_rna():
    """Convert DNA string to RNA by replacing T with U from rosalind_rna.txt"""
    path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\rosalind_rna.txt"
    with open(path, 'r') as file:
        text = file.read().strip()
    
    out = ""
    for i in range(len(text)):
        if text[i] == 'T':
            out += 'U'
        else:
            out += text[i]
    
    print(out)
    return out

def reverse_complement():
    """Generate reverse complement of DNA string from rosalind_revc.txt"""
    path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\rosalind_revc.txt"
    with open(path, 'r') as file:
        text = file.read().strip()
    
    out = ""
    complements = {'A': 'T', 'G': 'C', "T": "A", "C": "G"}
    for i in range(len(text) - 1, -1, -1):
        if text[i] in complements.keys():
            out += complements[text[i]]
        else:
            out += text[i]
    
    print(out)
    return out

def fibonacci_rabbits(n=None, k=None):
    """Calculate rabbit population with given parameters (wascally wabbits problem)"""
    if n is None or k is None:
        # If no parameters provided, use default values from original code
        n, k = 33, 5
    
    def fib(n, k):
        if n == 1:
            return 1
        if n <= 0:
            return 0
        return fib(n - 1, k) + k * fib(n - 2, k)
    
    result = fib(n, k)
    print(result)
    return result

def highest_gc_content():
    """Find sequence with highest GC content from rosalind_gc.txt"""
    path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\rosalind_gc.txt"
    with open(path, 'r') as file:
        full_text = file.read()
    
    temp = full_text.split(">")
    
    maxCG = -1
    name = ""
    for i in temp:
        if i == "\n" or i == "":
            continue
        else:
            lines = i.strip().split('\n')
            head = lines[0]
            body = ''.join(lines[1:])
            
            count = 0
            filtered = ""
            
            for j in body:
                if j == "G" or j == "C":
                    count += 1
                if j != "\n":
                    filtered += j
            
            if len(filtered) > 0:
                cur = round(count / len(filtered) * 100, 6)
                if cur > maxCG:
                    maxCG = cur
                    name = head
    
    print(name)
    print(maxCG)
    return name, maxCG
def hamming_distance():
    """Calculate Hamming distance between two DNA strings"""
    # Using the strings from your original code - you can modify to read from file later
    string1 = "ATAGTGAAACACGTGTACGCCGAGCCTAATGATGCACGCCCCAGGACCCTACAACAACCGCGTAGTTAGCGGGTAGCCGTTCTCCCTCTATACAGGACAGGAGGTCCCACCGTGACGTAATACTATATGTGGTTATTACCTTATCTCGTAAGATGAGTAGCTTGCATCGATCTCAAGGTCTTTTATTAGATGAAGGGAACAAAGGGTCTTTCACGAACGCCGTTCCACTCCTGTATATGTGGAGACGGGGCACGATAACATGTATCCCGATCGGTCCGTTTCATGTCTGATCGAAGCTAGTTCCGAGGTAGCTAAAATATCTGAGTCATTGCCCTAAACGCACGCCCTTTTGATCGCCCTGAGATAACTCCGAATGTGAGAGTTTCAGTGTCTTCCCAGAATTCGCGAGACAAATCGATCATCTTAACTACGTATCAATGGCATTCCATGCCCCTCAGATTTTGGGTAACCCTCATCTATCTCCTTCTTCAATGTTTATACAAGACTACTGTGTGACGTCAGATGCAAGTGGATTAGTTAACTCAGGTCTTATCATGCCCGGTCTATCGGCTCGGATGCCTAACAATCGTAGGGCTTAGTGCTCGGGACGGACCCTGAGACAGGCTGATAAAACCGTCGCGCCTGTGTCCGTCAACAGTAAATCCTATTGGCTGCAGAGTACGTCGAGCGCCGTGGTTGGGTACATGCAACAATCGACGTGGTGAATTTTGCGGGCCTATGCGTCATACCGAACAGCATCAAAGTACTAACCCCTCCCAATAGCGATTCGACCTGTGTCGGTGCCCGACGGCGACATCTCTTGGGTGCGGGAAACACGAGAGTTCCCGTTAGTAGATGTCGTTGCATCCTTAAGGTTGGATACCGACAAAGGTCCTGTACTTCACCACAGGCAGACCATGGTGTCCATACTTAGGGCACTGTGTAGTGTAAGGGGTCGCCGGCTCGCATCGGGATTGGAATTCGCGGAAATT"
    string2 = "CTAGTTACGGGCTAGCTCCGTGGGACCGCTCAGCTCCGCGCCTACCCCCTGCATCGCTTACATAGGTCGAGGAGAACCGCACCGTTTATACCTCGGAAAGGCAACCAAAACGGGATGTAATGCTGAATGGAGTCGCGCGTTTATCTCTAAGGAAGCTTCGCTAGCATGCCTCTCCCTGACGTTAATCCGTAGAAGGGAATCCAGGATACTTCAGGATTGGGCTTCGCTTCATTTATCCGCATAGGGGGGGTCGAACTGCATCCACTAACGCCGACACTTTTCATGTAAAATTTAACCAGATTCCGATGCACATACGATGCATGGTAAGTCGAACATAACAAAGGTAATTCTGTTGCCCCTCAGATATCGCTAAAAGGTAGAGGTTAAATACCCTTTGGGAATTCCTCGCACAAAGGGTGCAGCTAGGCATCTGATGAATTTCTCTCCTCCCTCATCAAATATGGGGGAACACTTACCCTTCAGTGGCGTCCATTGTAATTCCTCACAAGTCAGTATTGAGAGATCCACTACGTTAAGTTAATCACGCGATGATCTTGCAAGCTCTCCCATCGCGAATGGGCCGCATGCATAGTGAGTCGCTATCGGGACGAACACGGAACGCGGCACCTATAACAACCCCGCCTTGTTATTACATCACTGACTGCGTTTTTCCCCAAGACACATGGCACGCAAAATTTCGTATAATGCTACATACGCCGAACGCAATTAAGCTTAGGTAGGCGTCACACAGGACAGGATAGAAGTCAGAGTATTCATAATAAATGTACCGACTCGAGAACCAGACGGAGACCGACCTCTACTCGGGCAGGCACAATGAAAAGGCCAGCTTTCTTGGGGGACTAGCGCGAATGTGCTTGCATGCCGAAGACCGGGTTGTAGTGCGCTACAAGTCCGGCCGTGTGTGATCTCCGACGTAAGTAACTATGGTACGGGCTAACTAGCTAAGCACGGAAGTCTACTTAGCGTCCCTT"
    
    counter = len(string1)
    c = 0
    for i in range(0, counter):
        if string1[i] != string2[i]:
            c += 1
    
    print(c)
    return c

def hamming_distance_from_file():
    """Calculate Hamming distance between two DNA strings from file"""
    # You can create this function to read from a hamming distance dataset file
    path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\rosalind_hamm.txt"
    with open(path, 'r') as file:
        lines = file.read().strip().split('\n')
        string1 = lines[0]
        string2 = lines[1]
    
    c = 0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            c += 1
    
    print(c)
    return 

def mendel_dominant_probability(k=None, m=None, n=None):
    """Calculate probability of dominant phenotype in Mendelian genetics"""
    if k is None or m is None or n is None:
        # If no parameters provided, get from user input
        k, m, n = map(int, input("Enter k m n: ").split())
    
    total = k + m + n
    total_possible_matings = ((k + m + n) * (k + m + n - 1)) / 2
    population_of_total_individuals = ((k * (k - 1) / 2) + k * m + k * n + 
                                     m * (m - 1) / 2 * 0.75 + m * n * 0.5 + 
                                     n * (n - 1) / 2 * 0.0)
    probability_of_dominant_phenotype = population_of_total_individuals / total_possible_matings
    
    print(probability_of_dominant_phenotype)
    return probability_of_dominant_phenotype

def mendel_dominant_probability_from_file():
    """Calculate probability of dominant phenotype from file"""
    # You can create this function to read from a mendel dataset file
    path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\rosalind_iprb.txt" 
    with open(path, 'r') as file:
        k, m, n = map(int, file.read().strip().split())
    return mendel_dominant_probability(k, m, n)
    pass

def translate_rna_to_protein():
    """Translate RNA string to protein using codon table"""
    # Load codon table
    codon_table_path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\codon table.txt"
    with open(codon_table_path, 'r') as file:
        codon_table_text = file.read()
    
    # Parse codon table
    pattern = r'([AUGC]{3})\s+([A-Z][a-z]*|Stop)'
    matches = re.findall(pattern, codon_table_text)
    codon_table = dict(matches)
    
    # Load RNA string
    rna_path = r"C:\Users\Sadnan\Documents\University\Rosalind_Coding\datasets\rosalind_prot.txt"
    with open(rna_path, 'r') as file:
        dna_string = file.read().strip()
    
    c = ""
    out = ""
    for i in range(len(dna_string)):
        c += dna_string[i]
        
        if len(c) == 3:
            if c in codon_table:
                amino_acid = codon_table[c]
                if amino_acid == "Stop":
                    break
                out += amino_acid
            c = ""
    
    print(out)
    return out
