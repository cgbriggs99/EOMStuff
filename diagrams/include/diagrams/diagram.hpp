/*
 * diagram.hpp
 *
 *  Created on: Jun 3, 2026
 *      Author: connor
 */

#ifndef DIAGRAM_HPP_
#define DIAGRAM_HPP_

#include <list>
#include <string>
#include <tuple>
#include <limits>

namespace diagram {

using triplet = std::tuple<unsigned short, unsigned char, unsigned char>;

std::string triple_string(triplet const &triple, int start_i, int start_a, int start_m, int start_e);

enum Interaction : unsigned char {
    Fab = 0, // f_ab
    Fji = 1, // f_ji
    Fia = 2, // f_ia
    Wadbc = 3, // <ad||bc>
    Wjkil = 4, // <jk||il>
    Wjaib = 5, // <ja||ib>
    Wicab = 6, // <ic||ab>
    Wijak = 7, // <ij||ak>
    Wacib = 8, // <ac||ib>
    Wajik = 9, // <aj||ik>
    Wijab = 10, // <ij||ab>
    Fai = 11, // f_ai
    Wabij = 12 // <ab||ij>
};

// Based on the work of Kallay and Surjan.
class Diagram {
public:
    constexpr Diagram(Interaction interaction, triplet const &t1, triplet const &t2 = { 0, 0, 0 }, triplet const &t3 = { 0, 0, 0 },
            triplet const &t4 = { 0, 0, 0 }) : interaction { interaction }, triple_1 { t1 }, triple_2 { t2 }, triple_3 { t3 }, triple_4 { t4 } {
    }

    constexpr Interaction get_interaction() const {
        return interaction;
    }

    constexpr triplet get_t1() const {
        return triple_1;
    }
    constexpr triplet get_t2() const {
        return triple_2;
    }
    constexpr triplet get_t3() const {
        return triple_3;
    }
    constexpr triplet get_t4() const {
        return triple_4;
    }

    std::string to_string() const;

private:
    triplet triple_1, triple_2, triple_3, triple_4;
    Interaction interaction;
};

class RightDiagram {
public:
    constexpr RightDiagram(Interaction interaction, triplet const &right, signed short particle_excess = 0, triplet const &t1 = { 0, 0, 0 },
            triplet const &t2 = { 0, 0, 0 }, triplet const &t3 = { 0, 0, 0 }) : interaction { interaction }, triple_1 { t1 }, triple_2 { t2 }, triple_3 {
            t3 }, right { right }, particle_excess { particle_excess } {
    }

    constexpr Interaction get_interaction() const {
        return interaction;
    }

    constexpr signed short get_particle_excess() const {
        return particle_excess;
    }

    constexpr triplet get_t1() const {
        return triple_1;
    }
    constexpr triplet get_t2() const {
        return triple_2;
    }
    constexpr triplet get_t3() const {
        return triple_3;
    }
    constexpr triplet get_right() const {
        return right;
    }

private:
    triplet triple_1, triple_2, triple_3, right;
    signed short particle_excess;
    Interaction interaction;
};

class LambdaDiagram {
public:

    constexpr LambdaDiagram(Interaction interaction, unsigned short lambda_level, triplet const &t1, triplet const &t2 = { 0, 0, 0 },
            triplet const &t3 = { 0, 0, 0 }) : interaction { interaction }, lambda { lambda_level }, triple_1 { t1 }, triple_2 { t2 }, triple_3 {
            t3 } {
    }

    constexpr Interaction get_interaction() const {
        return interaction;
    }

    constexpr unsigned short get_lambda() const {
        return lambda;
    }

    constexpr triplet get_t1() const {
        return triple_1;
    }

    constexpr triplet get_t2() const {
        return triple_2;
    }

    constexpr triplet get_t3() const {
        return triple_3;
    }

    std::string to_string() const;

private:
    triplet triple_1, triple_2, triple_3;
    unsigned short lambda;
    Interaction interaction;
};

namespace detail {
void assign_Fab_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Fij_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Fia_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Fia_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Fia_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wabcd_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wabcd_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wabcd_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijkl_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijkl_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijkl_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Waijb_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Waijb_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Waijb_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wiabc_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wiabc_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wiabc_t3_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wiabc_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijak_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijak_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijak_t3_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijak_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wabic_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Waijk_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijab_t_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijab_t2_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijab_t3_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijab_t4_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1,
        int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Wijab_terms(int residual, int max_excitation, std::list<Diagram> *out, bool skip_t1, int total_order =
        std::numeric_limits<int>::max(), int t1_order = 2);

void assign_Fai_terms(int residual, std::list<Diagram> *out);

void assign_Wabij_terms(int residual, std::list<Diagram> *out);
}

// skip_t1 = true and canonical = true gives the T1 transformed form.
std::list<Diagram> compute_cc_residual(int residual_exc, int max_exc, bool skip_t1, bool canonical);

std::list<Diagram> compute_ccn_residual(int residual, int max_excitation, bool skip_t1, bool canonical);

std::list<Diagram> compute_ccp_residual(int residual_exc, int max_exc, bool skip_t1, bool canonical);

std::list<Diagram> compute_lcc_residual(int residual_exc, int max_exc, bool skip_t1, bool canonical);

} // namespace diagram

#endif /* DIAGRAM_HPP_ */
