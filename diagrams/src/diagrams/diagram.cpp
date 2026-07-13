/*
 * diagram.cpp
 *
 *  Created on: Jun 3, 2026
 *      Author: connor
 */

#include <diagrams/diagram.hpp>
#include <sstream>
#include <limits>

namespace diagram {

std::string triple_string(triplet const &triple, int start_i, int start_a, int start_m, int start_e) {

    std::ostringstream out;
    if (std::get < 0 > (triple) != 0) {
        if (std::get < 0 > (triple) > 0) {
            out << "t(";
        } else {
            out << "r(";
        }

        int internal_lines = std::get < 1 > (triple);

        int internal_particles = std::get < 2 > (triple);

        int internal_holes = internal_lines - internal_particles;

        int external_particles = std::get < 0 > (triple) - internal_particles;

        int external_holes = std::get < 0 > (triple) - internal_holes;

        bool begin = true;

        for (int i = 0; i < external_holes; i++) {
            if (begin) {
                out << "i" << i + start_i;
                begin = false;
            } else {
                out << ", i" << i + start_i;
            }
        }

        for (int i = 0; i < internal_holes; i++) {
            if (begin) {
                out << "m" << i + start_m;
                begin = false;
            } else {
                out << ", m" << i + start_m;
            }
        }

        for (int a = 0; a < external_particles; a++) {
            out << ", a" << a + start_a;
        }

        for (int a = 0; a < internal_particles; a++) {
            out << ", e" << a + start_e;
        }

        out << ")";
    }

    return out.str();
}

std::string Diagram::to_string() const {

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

    out << triple_string(triple_4, curr_i_ind, curr_a_ind, curr_m_ind, curr_e_ind);
    internal_lines = std::get < 1 > (triple_4);

    internal_particles = std::get < 2 > (triple_4);

    internal_holes = internal_lines - internal_particles;

    external_particles = std::get < 0 > (triple_4) - internal_particles;

    external_holes = std::get < 0 > (triple_4) - internal_holes;

    curr_i_ind += external_holes;
    curr_a_ind += external_particles;
    curr_m_ind = 1;
    curr_e_ind = 1;

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
    curr_m_ind = 1;
    curr_e_ind = 1;

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

FactoredDiagram::FactoredDiagram(std::list<Diagram> const &diagrams) : curr_triple_ { 0, 0, 0 } {
    for (auto diagram : diagrams) {
        add_diagram(diagram);
    }
}

FactoredDiagram::FactoredDiagram(int level) : level_ { level } {
}

void FactoredDiagram::add_diagram(Diagram const &diagram) {
    bool found = false;

    if (level_ >= 3) {
        leaves_.push_back(diagram);
        return;
    }

    if (diagram.get(level_) == triplet { 0, 0, 0 }) {
        leaves_.push_back(diagram);
        return;
    }

    for (auto &branch : branches_) {
        if (branch.get_triplet() == diagram.get(level_)) {
            branch.add_diagram(diagram);
            found = true;
            break;
        }
    }

    if (!found) {
        auto &new_factor = branches_.emplace_back(level_ + 1);

        new_factor.add_diagram(diagram);
    }
}

} // namespace diagram
