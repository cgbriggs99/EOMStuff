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
#include <cstdint>

namespace diagram {

using triplet = std::tuple<signed short, unsigned char, unsigned char>;

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

    constexpr triplet get(int index) const {
        switch (index) {
        case 0:
            return triple_1;
        case 1:
            return triple_2;
        case 2:
            return triple_3;
        case 3:
            return triple_4;
        default:
            return triple_1;
        }
    }

    std::string to_string() const;

private:
    triplet triple_1, triple_2, triple_3, triple_4;
    Interaction interaction;
};

constexpr bool operator==(triplet const &left, triplet const &right) {
    return std::get < 0 > (left) == std::get < 0 > (right) && std::get < 1 > (left) == std::get < 1 > (right)
            && std::get < 2 > (left) == std::get < 2 > (right);
}

constexpr bool operator==(Diagram const &left, Diagram const &right) {
    return left.get_t1() == right.get_t1() && left.get_t2() == right.get_t2() && left.get_t3() == right.get_t3()
            && left.get_t4() == right.get_t4() && left.get_interaction() == right.get_interaction();
}

class FactoredDiagram {
public:
    FactoredDiagram() = default;

    FactoredDiagram(int level);

    FactoredDiagram(std::list<Diagram> const &diagrams);

    [[nodiscard]] inline std::list<FactoredDiagram> const& get_branches() const {
        return branches_;
    }

    [[nodiscard]] inline std::list<FactoredDiagram>& get_branches() {
        return branches_;
    }

    [[nodiscard]] inline std::list<Diagram> const& get_leaves() const {
        return leaves_;
    }

    [[nodiscard]] inline std::list<Diagram>& get_leaves() {
        return leaves_;
    }

    [[nodiscard]] inline triplet get_triplet() const {
        return curr_triple_;
    }

    [[nodiscard]] inline triplet& get_triplet() {
        return curr_triple_;
    }

    [[nodiscard]] inline int get_level() const {
        return level_;
    }

    void add_diagram(Diagram const &diagram);

private:
    std::list<FactoredDiagram> branches_;

    std::list<Diagram> leaves_;

    triplet curr_triple_ { 0, 0, 0 };
    int level_ { 0 };
};

//class LambdaDiagram {
//public:
//
//    constexpr LambdaDiagram(Interaction interaction, unsigned short lambda_level, triplet const &t1, triplet const &t2 = { 0, 0, 0 },
//            triplet const &t3 = { 0, 0, 0 }) : interaction { interaction }, lambda { lambda_level }, triple_1 { t1 }, triple_2 { t2 }, triple_3 {
//            t3 } {
//    }
//
//    constexpr Interaction get_interaction() const {
//        return interaction;
//    }
//
//    constexpr unsigned short get_lambda() const {
//        return lambda;
//    }
//
//    constexpr triplet get_t1() const {
//        return triple_1;
//    }
//
//    constexpr triplet get_t2() const {
//        return triple_2;
//    }
//
//    constexpr triplet get_t3() const {
//        return triple_3;
//    }
//
//    std::string to_string() const;
//
//private:
//    triplet triple_1, triple_2, triple_3;
//    unsigned short lambda;
//    Interaction interaction;
//};
//
//class DensityDiagram {
//public:
//    DensityDiagram(uint32_t lambda, std::list<triplet> const &triples);
//    ~DensityDiagram();
//
//    uint32_t get_lambda() const;
//
//    triplet const* begin() const;
//    triplet* begin();
//    triplet const* cbegin() const;
//
//    triplet const* end() const;
//    triplet* end();
//    triplet const* cend() const;
//
//    std::string to_string() const;
//
//private:
//    uint32_t lambda_order { 0 }, num_t_amps { 0 };
//    triplet *triples { nullptr };
//};

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

// Lambda assignments
/*
 void assign_lFab_terms(int residual, int level, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lFij_terms(int residual, int level, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lFia_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lFia_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lFia_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1,
 int total_order = std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWabcd_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWabcd_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijkl_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijkl_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWaijb_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWaijb_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWiabc_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWiabc_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWiabc_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijak_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijak_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijak_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWabic_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWaijk_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijab_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijab_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijab_t3_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);

 void assign_lWijab_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order =
 std::numeric_limits<int>::max(), int t1_order = 2);
 */
}

// skip_t1 = true and canonical = true gives the T1 transformed form.
std::list<Diagram> compute_cc_residual(int residual_exc, int max_exc, bool skip_t1, bool canonical);

std::list<Diagram> compute_ccn_residual(int residual, int max_excitation, bool skip_t1, bool canonical);

std::list<Diagram> compute_ccp_residual(int residual_exc, int max_exc, bool skip_t1, bool canonical);

std::list<Diagram> compute_lcc_residual(int residual_exc, int max_exc, bool skip_t1, bool canonical);

/*
 std::list<LambdaDiagram> compute_lambda_cc_residual(int residual, int level, bool skip_t1, bool canonical);

 std::list<LambdaDiagram> compute_lambda_ccn_residual(int residual, int level, bool skip_t1, bool canonical);

 std::list<LambdaDiagram> compute_lambda_ccp_residual(int residual, int level, bool skip_t1, bool canonical);

 std::list<DensityDiagram> compute_denstity_block(unsigned int level, unsigned int hole_creators, unsigned int particle_annihilators,
 unsigned int particle_creators, unsigned int hole_annihilators);
 */

} // namespace diagram

#endif /* DIAGRAM_HPP_ */
