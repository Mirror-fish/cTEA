/**
 * @file complementary_seq.hpp
 * @brief Generate complementary strand sequences
 * @author cTEA Development Team
 * 
 * From MEGAnE's complementary_seq.hpp
 * Used for reverse complement in read processing
 */

#ifndef COMPLEMENTARY_SEQ_HPP
#define COMPLEMENTARY_SEQ_HPP

#include <string>
#include <unordered_map>

namespace complementary_seq_hpp {

/**
 * @brief Get complementary base
 */
inline char complement(char base) {
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'N': case 'n': return 'N';
        default: return base;  // Unknown, return as-is
    }
}

/**
 * @brief Reverse complement a sequence
 * (Simplified from MEGAnE's implementation)
 */
inline std::string reverse_complement(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.length());
    
    for (int i = seq.length() - 1; i >= 0; i--) {
        rc += complement(seq[i]);
    }
    
    return rc;
}

/**
 * @brief Complement a sequence (without reversing)
 * MEGAnE uses this for processing soft-clipped reads
 */
inline std::string complement_string(const std::string& seq) {
    std::string c;
    c.reserve(seq.length());
    
    for (char base : seq) {
        c += complement(base);
    }
    
    return c;
}

/**
 * @brief Lookup table for faster complement
 * Fixed: Removed const from array to allow initialization in constructor.
 */
static char complement_table[256];

// Initialize the table properly
struct ComplementInitializer {
    ComplementInitializer() {
        // Initialize all to 0 first
        for (int i = 0; i < 256; i++) {
            complement_table[i] = 0;
        }
        // A, a -> T
        complement_table[(unsigned char)'A'] = 'T';
        complement_table[(unsigned char)'a'] = 'T';
        // C, c -> G
        complement_table[(unsigned char)'C'] = 'G';
        complement_table[(unsigned char)'c'] = 'G';
        // G, g -> C
        complement_table[(unsigned char)'G'] = 'C';
        complement_table[(unsigned char)'g'] = 'C';
        // T, t -> A
        complement_table[(unsigned char)'T'] = 'A';
        complement_table[(unsigned char)'t'] = 'A';
        // N, n -> N
        complement_table[(unsigned char)'N'] = 'N';
        complement_table[(unsigned char)'n'] = 'N';
    }
};

static ComplementInitializer complement_initializer;

} // namespace complementary_seq_hpp

#endif // COMPLEMENTARY_SEQ_HPP