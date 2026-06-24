/*
 * orbit-string.cpp
 *
 *  Created on: Jun 23, 2026
 *      Author: connor
 */

#include <general_cc/orbit-string.hpp>
#include <general_cc/utils.hpp>
#include <vector>
#include <map>
#include <sstream>
#include <string>
#include <memory>
#include <limits>

namespace psi::general_cc {

OrbitString::OrbitString(std::vector<unsigned short> orbitals, unsigned short max_orbital, Spin spin) {
    // First, make sure there aren't too many orbitals.

    if (orbitals.size() > std::numeric_limits<unsigned short>::max()) {
        throw std::runtime_error("Too many orbitals passed to orbit string! Keep the number under 65536.");
    }

    if (orbitals.size() > 0 && orbitals.data() == nullptr) {
        throw std::runtime_error("Input list of orbitals is in an undefined state! Size is non-zero, but data is null.");
    }

    // Next, ensure that all orbitals are within the specified range.

    for (unsigned short p : orbitals) {
        if (p >= max_orbital) {
            throw std::runtime_error("One of the orbitals passed into the orbit string was too large!");
        }
    }

    // Next, make sure that the orbitals are ascending.
    for (int i = 1; i < orbitals.size(); i++) {
        if (orbitals[i - 1] >= orbitals[i]) {
            throw std::runtime_error("The orbitals are not in ascending order!");
        }
    }

    // Also, we can't have more indices than the max orbital due to the restriction that the
    // indices are in strictly ascending order.

    // Sanity checks over. Now we make the object.

    length_ = static_cast<unsigned short>(orbitals.size());

    indices_ = new unsigned short[length_];

    spin_ = spin;

    max_orbital_ = max_orbital;

    if (indices_ == nullptr) {
        throw std::runtime_error("Could not allocate orbit string!");
    }

    std::memcpy(indices_, orbitals.data(), length_ * sizeof(unsigned short));
}

OrbitString::OrbitString(unsigned short length, unsigned short max_orbital, Spin spin) : length_ { length }, max_orbital_ { max_orbital }, spin_ {
        spin } {

    if (length_ >= max_orbital_) {
        throw std::runtime_error("Determinants with more indices than orbitals will be zero!");
    }

    indices_ = new unsigned short[length_];

    if (indices_ == nullptr) {
        throw std::runtime_error("Could not allocate orbit string.");
    }

    for (unsigned short i = 0; i < length_; i++) {
        indices_[i] = i;
    }
}

OrbitString::OrbitString(size_t index, unsigned short length, unsigned short max_orbital, Spin spin) : spin_ { spin }, max_orbital_ {
        max_orbital }, length_ { length } {

    if (index >= number_of_strings(length, max_orbital)) {
        throw std::runtime_error("The index is out of range for the given length and maximum orbital index.");
    }

    detail::binomial_storage_type dynamic_list;

    indices_ = new unsigned short[length_];

    if (indices_ == nullptr) {
        throw std::runtime_error("Could not allocate orbit string.");
    }

    size_t curr_index = index;

    std::exception_ptr eptr;

    try {
        for (unsigned short k = length_; k > 0; k--) {
            // Loop through and try to find the needed element.
            if (curr_index == 0) {
                indices_[k - 1] = k - 1;
                continue;
            }
            signed short ck = max_orbital_ - 1;

            if (k < length_) {
                ck = indices_[k] - 1;
            }
            indices_[k - 1] = k - 1;
            for (; ck >= k - 1; ck--) {
                size_t binomial = binomial_coefficient_driver(ck, k, &dynamic_list);

                if (binomial <= curr_index) {
                    curr_index -= binomial;
                    indices_[k - 1] = ck;
                    break;
                }
            }
        }
    } catch (...) {
        delete[] indices_;

        eptr = std::current_exception();
    }

    if (eptr) {
        std::rethrow_exception(eptr);
    }
}

OrbitString::OrbitString(OrbitString const &other) : length_ { other.get_num_indices() }, max_orbital_ { other.get_max_orbital() }, spin_ {
        other.get_spin() } {
    if (!other.sanity_check()) {
        throw std::runtime_error("Copying a malformed orbit string!");
    }

    indices_ = new unsigned short[length_];

    std::memcpy(indices_, other.get_indices(), length_ * sizeof(unsigned short));
}

OrbitString::OrbitString(OrbitString &&other) : length_ { other.get_num_indices() }, max_orbital_ { other.get_max_orbital() }, spin_ {
        other.get_spin() } {
    if (!other.sanity_check()) {
        throw std::runtime_error("Moving a malformed orbit string!");
    }

    indices_ = other.indices_;

    other.length_ = 0;
    other.max_orbital_ = 0;
    other.spin_ = ALPHA;
    other.indices_ = nullptr;
}

OrbitString& OrbitString::operator=(OrbitString const &other) {

    if (!other.sanity_check()) {
        throw std::runtime_error("Copying a malformed orbit string!");
    }
    length_ = other.get_num_indices();
    max_orbital_ = other.get_max_orbital();
    spin_ = other.get_spin();

    if (other.get_indices() == nullptr && length_ != 0) {
        throw std::runtime_error("Copy value is in an undefined state where its pointer is null but its length is non-zero.");
    }

    indices_ = new unsigned short[length_];

    std::memcpy(indices_, other.get_indices(), length_ * sizeof(unsigned short));

    return *this;
}

OrbitString& OrbitString::operator=(OrbitString &&other) {
    if (!other.sanity_check()) {
        throw std::runtime_error("Moving a malformed orbit string!");
    }

    length_ = other.get_num_indices();
    max_orbital_ = other.get_max_orbital();
    spin_ = other.get_spin();
    indices_ = other.get_indices();

    other.length_ = 0;
    other.max_orbital_ = 0;
    other.spin_ = ALPHA;
    other.indices_ = nullptr;

    return *this;
}

OrbitString::~OrbitString() {
    length_ = 0;
    max_orbital_ = 0;
    spin_ = ALPHA;

    if (indices_ != nullptr) {
        delete[] indices_;
    }
}

unsigned short OrbitString::get_index(unsigned short i) const {
    if (i < length_) {
        return indices_[i];
    } else {
        throw std::length_error("Index out of range");
    }
}

unsigned short& OrbitString::get_index(unsigned short i) {
    if (i < length_) {
        return indices_[i];
    } else {
        throw std::length_error("Index out of range");
    }
}

OrbitString& OrbitString::operator++() {
    unsigned short carry = 1;

    for (unsigned short i = 0; i < length_ && carry > 0; i++) {
        indices_[i] += carry;
        carry = 0;

        if (length_ > 0 && i < length_ - 1 && indices_[i] >= indices_[i + 1]) {
            if (i > 0) {
                indices_[i] = indices_[i - 1] + 1;
            } else {
                indices_[i] = 0;
            }
            carry = 1;
        }
    }

    if (length_ > 0 && indices_[length_ - 1] >= max_orbital_) {
        throw std::runtime_error("Went too far! The highest index went past the highest orbital.");
    }

    return *this;
}

OrbitString OrbitString::operator++(int) {
    OrbitString out = *this;
    ++*this;

    return out;
}

bool OrbitString::has_next() const {
    if (!sanity_check()) {
        return false;
    }

    if (indices_ != nullptr) {
        return indices_[0] < max_orbital_ - length_;
    } else {
        return false;
    }
}

size_t OrbitString::number_of_strings(unsigned short length, unsigned short max_orbital) {
// There are max_orbital balls and length bins.
    return binomial_coefficient(max_orbital, length);
}

size_t OrbitString::to_index() const {
// Apparently, this is relatively easy. You just put the value in the top of the binomial expansion and
// the position in the bottom. The top is zero-indexed, but the bottom is 1-indexed.
    detail::binomial_storage_type dynamic_list;

    size_t index = 0;

    for (unsigned short i = 0; i < length_; i++) {
        index += detail::binomial_coefficient_driver(indices_[i], i + 1, &dynamic_list);
    }

    return index;
}

bool OrbitString::sanity_check() const {
    if (indices_ == nullptr && length_ != 0) {
        return false;
    }

    if (length_ >= max_orbital_) {
        return false;
    }

    if (spin_ != ALPHA && spin_ != BETA) {
        return false;
    }

    for (unsigned short i = 0; i < length_; i++) {
        if (indices_[i] >= max_orbital_) {
            return false;
        }
    }

    if (length_ > 1) {
        for (unsigned short i = 0; i < length_ - 1; i++) {
            if (indices_[i] >= indices_[i + 1]) {
                return false;
            }
        }
    }
    return true;
}

std::string OrbitString::to_string() const {
    std::ostringstream out;

    out << "(Spin: ";

    switch (spin_) {
    case ALPHA:
        out << "ALPHA";
        break;
    case BETA:
        out << "BETA";
        break;
    default:
        out << "UNKNOWN";
    }

    out << ", Length: " << length_ << ", Max value: " << max_orbital_ << ") {";

    for (unsigned short i = 0; i < length_; i++) {
        out << indices_[i];
        if (i != length_ - 1) {
            out << ", ";
        }
    }
    out << '}';

    return out.str();
}

bool operator<(OrbitString const &left, OrbitString const &right) {
    if (left.get_num_indices() < right.get_num_indices()) {
        return true;
    }

    if (left.get_num_indices() > right.get_num_indices()) {
        return false;
    }

    if (left.get_indices() == nullptr) {
        if (right.get_indices() == nullptr) {
            return false;
        }
        return true;
    }

    if (right.get_indices() == nullptr) {
        return false;
    }

    if (left.get_spin() < right.get_spin()) {
        return true;
    }

    if (left.get_spin() > right.get_spin()) {
        return false;
    }

    if (left.get_max_orbital() < right.get_max_orbital()) {
        return true;
    }

    if (left.get_max_orbital() > right.get_max_orbital()) {
        return false;
    }

    for (unsigned short i = 0; i < left.get_num_indices(); i++) {
        if (left.get_index(i) < right.get_index(i)) {
            return true;
        }
        if (left.get_index(i) > right.get_index(i)) {
            return false;
        }
    }
    return false;
}

bool operator==(OrbitString const &left, OrbitString const &right) {
    if (left.get_num_indices() != right.get_num_indices()) {
        return false;
    }

    if (left.get_spin() != right.get_spin()) {
        return false;
    }

    if (left.get_max_orbital() != right.get_max_orbital()) {
        return false;
    }

    // Shortcut. If the pointers are the same, the indices are the same.
    if (left.get_indices() == right.get_indices()) {
        return true;
    }

    for (unsigned short i = 0; i < left.get_num_indices(); i++) {
        if (left.get_index(i) != right.get_index(i)) {
            return false;
        }
    }

    return true;
}

}

