/*
 * bit-string.hpp
 *
 *  Created on: Jun 23, 2026
 *      Author: connor
 */

#ifndef GENERAL_CC_INCLUDE_BIT_STRING_HPP_
#define GENERAL_CC_INCLUDE_BIT_STRING_HPP_

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"

#include <vector>
#include <string>

#ifndef PSI_API
#define PSI_API
#endif

namespace psi {
namespace general_cc {

class PSI_API OrbitString final {
public:
    enum Spin : unsigned short {
        ALPHA = 0, BETA = 1
    };

    // We use unsigned short for length since we don't expect to have more than 65536 orbitals, but could reasonably have more than 256.
    // We want to pack things well, so we need to cut out unnecessary bytes.
    OrbitString(std::vector<unsigned short> orbitals, unsigned short max_oribtal, Spin spin);

    OrbitString(unsigned short length, unsigned short max_orbital, Spin spin);

    OrbitString(size_t index, unsigned short length, unsigned short max_orbital, Spin spin);

    OrbitString() = default;

    OrbitString(OrbitString const &other);

    OrbitString(OrbitString &&other);

    OrbitString& operator=(OrbitString const &other);

    OrbitString& operator=(OrbitString &&other);

    ~OrbitString();

    [[nodiscard]] inline unsigned short const* get_indices() const {
        return indices_;
    }

    [[nodiscard]] inline unsigned short* get_indices() {
        return indices_;
    }

    [[nodiscard]] unsigned short get_index(unsigned short i) const;

    [[nodiscard]] unsigned short& get_index(unsigned short i);

    [[nodiscard]] inline unsigned short& operator[](unsigned short i) {
        return indices_[i];
    }

    [[nodiscard]] inline unsigned short operator[](unsigned short i) const {
        return indices_[i];
    }

    [[nodiscard]] inline unsigned short get_num_indices() const {
        return length_;
    }

    [[nodiscard]] inline unsigned short get_max_orbital() const {
        return max_orbital_;
    }

    [[nodiscard]] inline Spin get_spin() const {
        return spin_;
    }

    OrbitString& operator++();

    OrbitString operator++(int);

    [[nodiscard]] bool has_next() const;

    [[nodiscard]] static size_t number_of_strings(unsigned short length, unsigned short max_orbital);

    [[nodiscard]] inline size_t number_of_strings() const {
        return number_of_strings(length_, max_orbital_);
    }

    [[nodiscard]] size_t to_index() const;

    [[nodiscard]] bool sanity_check() const;

    [[nodiscard]] std::string to_string() const;

    friend bool operator<(OrbitString const &left, OrbitString const &right);
    friend bool operator==(OrbitString const &left, OrbitString const &right);

private:
    // C-style array because I know how to do C-style memory management, and we may need the memory efficiency that
    // it provides.
    unsigned short *indices_ { nullptr };

    unsigned short length_ { 0 }; // number of excitations
    unsigned short max_orbital_ { 0 };

    Spin spin_ { ALPHA }; // Whether this is an alpha string or beta string.
    // I won't specify padding. The compiler can figure that out if it needs it.
};

bool operator<(OrbitString const &left, OrbitString const &right);

inline bool operator>(OrbitString const &left, OrbitString const &right) {
    return right < left;
}

inline bool operator<=(OrbitString const &left, OrbitString const &right) {
    return !(right < left);
}

inline bool operator>=(OrbitString const &left, OrbitString const &right) {
    return !(left < right);
}

bool operator==(OrbitString const &left, OrbitString const &right);

inline bool operator!=(OrbitString const &left, OrbitString const &right) {
    return !(left == right);
}

}
}

#endif /* GENERAL_CC_INCLUDE_BIT_STRING_HPP_ */
