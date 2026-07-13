/*
 * lambda_diagrams.cpp
 *
 *  Created on: Jul 8, 2026
 *      Author: connor
 */

#include <diagrams/diagram.hpp>
#include <string>
#include <list>

namespace diagram {

std::string LambdaDiagram::to_string() const {
    std::ostringstream out;

    if (interaction < 0 || interaction > 12) {
        throw std::runtime_error("Interaction id is outside of the expected range!");
    }

    int curr_i_ind = 1;
    int curr_a_ind = 1;
    int curr_m_ind = 1;
    int curr_e_ind = 1;

    out << triple_string(triple_1, curr_i_ind, curr_a_ind, curr_m_ind, curr_e_ind);

    int internal_lines = std::get < 1 > (triple_1);

    int internal_particles = std::get < 2 > (triple_1);

    int internal_holes = internal_lines - internal_particles;

    int external_particles = std::get < 0 > (triple_1) - internal_particles;

    int external_holes = std::get < 0 > (triple_1) - internal_holes;

    curr_i_ind += external_holes;
    curr_a_ind += external_particles;
    curr_m_ind += internal_holes;
    curr_e_ind += internal_particles;

    out << triple_string(triple_2, curr_i_ind, curr_a_ind, curr_m_ind, curr_e_ind);

    internal_lines = std::get < 1 > (triple_2);

    internal_particles = std::get < 2 > (triple_2);

    internal_holes = internal_lines - internal_particles;

    external_particles = std::get < 0 > (triple_2) - internal_particles;

    external_holes = std::get < 0 > (triple_2) - internal_holes;

    curr_i_ind += external_holes;
    curr_a_ind += external_particles;
    curr_m_ind += internal_holes;
    curr_e_ind += internal_particles;

    out << triple_string(triple_3, curr_i_ind, curr_a_ind, curr_m_ind, curr_e_ind);

    internal_lines = std::get < 1 > (triple_3);

    internal_particles = std::get < 2 > (triple_3);

    internal_holes = internal_lines - internal_particles;

    external_particles = std::get < 0 > (triple_3) - internal_particles;

    external_holes = std::get < 0 > (triple_3) - internal_holes;

    curr_i_ind += external_holes;
    curr_a_ind += external_particles;
    curr_m_ind += internal_holes;
    curr_e_ind += internal_particles;

    switch (interaction) {
    case Fab:
        out << "f(a" << curr_a_ind << ", e1)";
        break;
    case Fji:
        out << "f(m1, i" << curr_i_ind << ")";
        break;
    case Fia:
        out << "f(m1, e1)";
        break;
    case Wadbc:
        out << "<a" << curr_a_ind << " a" << curr_a_ind + 1 << "||e1 e2>";
        break;
    case Wjkil:
        out << "<m1 m2||i" << curr_i_ind << " i" << curr_i_ind + 1 << ">";
        break;
    case Wjaib:
        out << "<m1 a" << curr_a_ind << "||i" << curr_i_ind << " e1>";
        break;
    case Wicab:
        out << "<m1 a" << curr_a_ind << "||e1 e2>";
        break;
    case Wijak:
        out << "<m1 m2||e1 i" << curr_i_ind << ">";
        break;
    case Wacib:
        out << "<a" << curr_a_ind << " a" << curr_a_ind + 1 << "||i" << curr_i_ind << " m1>";
        break;
    case Wajik:
        out << "<a" << curr_a_ind << " m1||i" << curr_i_ind << " i" << curr_i_ind + 1 << ">";
        break;
    case Wijab:
        out << "<m1 m2||e1 e2>";
        break;
    case Fai:
        out << "f(a" << curr_a_ind << ", i" << curr_i_ind << ")";
        break;
    case Wabij:
        out << "<a" << curr_a_ind << " a" << curr_a_ind + 1 << "||i" << curr_i_ind << " i" << curr_i_ind + 1 << ">";
        break;
    default:
        throw std::runtime_error("Interaction id is outside of the expected range!");
    }

    return out.str();
}

namespace detail {
void assign_lFab_terms(int residual, int level, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order) {
    if(residual == 0) {
        return;
    }


}

void assign_lFij_terms(int residual, int level, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lFia_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lFia_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lFia_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWabcd_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWabcd_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijkl_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijkl_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWaijb_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWaijb_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWiabc_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWiabc_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWiabc_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijak_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijak_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijak_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWabic_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWaijk_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijab_t_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijab_t2_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijab_t3_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);

void assign_lWijab_terms(int residual, int max_excitation, std::list<LambdaDiagram> *out, bool skip_t1, int total_order, int t1_order);
}

}
