/**
 * @file dna_to_2bit.hpp
 * @brief DNA to 2-bit encoding (from MEGAnE)
 * @author cTEA Development Team
 * 
 * Optimized from MEGAnE's dna_to_2bit.hpp
 * Used for k-mer calculation in repeat filtering
 */

#ifndef DNA_TO_2BIT_HPP
#define DNA_TO_2BIT_HPP

#include <string>
#include <unordered_map>
#include <cstdint>

namespace dna_to_2bit_hpp {

/**
 * @brief Convert a single DNA character to 2-bit representation
 * A->00, T->11, C->01, G->10
 */
inline uint32_t dna_to_2bitf_32(char nt) {
    switch (nt) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;  // N or other -> treat as A
    }
}

/**
 * @brief Convert a string to 2-bit encoded k-mer
 * Used for rolling k-mer calculation
 */
inline uint32_t string_to_2bit(const std::string& seq, size_t start, size_t k) {
    uint32_t bit2 = 0;
    for (size_t i = 0; i < k && (start + i) < seq.length(); i++) {
        bit2 <<= 2;
        bit2 |= dna_to_2bitf_32(seq[start + i]);
    }
    return bit2;
}

/**
 * @brief Update 2-bit k-mer for rolling window
 * @param old_kmer: previous k-mer value
 * @param new_nt: new nucleotide character
 * @param k: k-mer size
 */
inline uint32_t update_2bit_kmer(uint32_t old_kmer, char new_nt, size_t k) {
    uint32_t new_kmer = old_kmer << 2;
    new_kmer |= dna_to_2bitf_32(new_nt);
    // Clear the high bits beyond k*2 bits
    new_kmer &= (1U << (k * 2)) - 1;
    return new_kmer;
}

/**
 * @brief Build a lookup table for DNA to 2-bit conversion
 * Faster than switch statement for large sequences.
 * 
 * Fixed: Removed const from array to allow initialization in constructor.
 */
static uint32_t dna_to_2bit_table[256];  // Non-const, initialized in constructor

// Initialize the table properly
struct DNATo2BitInitializer {
    DNATo2BitInitializer() {
        // Initialize all to 0 first
        for (int i = 0; i < 256; i++) {
            dna_to_2bit_table[i] = 0;
        }
        // A, a
        dna_to_2bit_table[(unsigned char)'A'] = 0;
        dna_to_2bit_table[(unsigned char)'a'] = 0;
        // C, c
        dna_to_2bit_table[(unsigned char)'C'] = 1;
        dna_to_2bit_table[(unsigned char)'c'] = 1;
        // G, g
        dna_to_2bit_table[(unsigned char)'G'] = 2;
        dna_to_2bit_table[(unsigned char)'g'] = 2;
        // T, t
        dna_to_2bit_table[(unsigned char)'T'] = 3;
        dna_to_2bit_table[(unsigned char)'t'] = 3;
    }
};

static DNATo2BitInitializer dna_to_2bit_initializer;

} // namespace dna_to_2bit_hpp

#endif // DNA_TO_2BIT_HPP