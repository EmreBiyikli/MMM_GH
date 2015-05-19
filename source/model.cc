// Copyright (C) 2011-2014 University of Pittsburgh. All rights reserved.
// See COPYING.txt for details.
// Author: Emre Biyikli (biyikli.emre@gmail.com)

// model.cc
// *****************************************************************************

#include "model.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "atom_group.h"
#include "matrix.h"
#include "mesh.h"
#include "mmm.h"
#include "output.h"
#include "potential.h"

using std::abs;
using std::string;
using std::swap;
using std::vector;

// Independent functions
// *****************************************************************************

namespace {

// determinant of 4x4
// *****************************************************************************
// (independent)
// Returns determinant of a 4-by-4 matrix.
double determinant44(double a[4][4]) {
  double det = 0.0;
  det += a[0][0] * a[1][1] * a[2][2] * a[3][3];
  det -= a[0][0] * a[1][1] * a[2][3] * a[3][2];
  det += a[0][0] * a[1][2] * a[2][3] * a[3][1];
  det -= a[0][0] * a[1][2] * a[2][1] * a[3][3];
  det += a[0][0] * a[1][3] * a[2][1] * a[3][2];
  det -= a[0][0] * a[1][3] * a[2][2] * a[3][1];
  det -= a[0][1] * a[1][2] * a[2][3] * a[3][0];
  det += a[0][1] * a[1][2] * a[2][0] * a[3][3];
  det -= a[0][1] * a[1][3] * a[2][0] * a[3][2];
  det += a[0][1] * a[1][3] * a[2][2] * a[3][0];
  det -= a[0][1] * a[1][0] * a[2][2] * a[3][3];
  det += a[0][1] * a[1][0] * a[2][3] * a[3][2];
  det += a[0][2] * a[1][3] * a[2][0] * a[3][1];
  det -= a[0][2] * a[1][3] * a[2][1] * a[3][0];
  det += a[0][2] * a[1][0] * a[2][1] * a[3][3];
  det -= a[0][2] * a[1][0] * a[2][3] * a[3][1];
  det += a[0][2] * a[1][1] * a[2][3] * a[3][0];
  det -= a[0][2] * a[1][1] * a[2][0] * a[3][3];
  det -= a[0][3] * a[1][0] * a[2][1] * a[3][2];
  det += a[0][3] * a[1][0] * a[2][2] * a[3][1];
  det -= a[0][3] * a[1][1] * a[2][2] * a[3][0];
  det += a[0][3] * a[1][1] * a[2][0] * a[3][2];
  det -= a[0][3] * a[1][2] * a[2][0] * a[3][1];
  det += a[0][3] * a[1][2] * a[2][1] * a[3][0];
  return det;
}

// gaussj
// *****************************************************************************
// (independent)
// Returns x of the linear matrix equation a * x = b for the input a and b
// solved by Gauss-Jordan elimination. a must be n-by-n, b must be n-by-m. m > 1
// for solving multiple cases. The returned version of a consists of inverse of
// input a whereas return version of b consists of solutions (i.e., x).
void gaussj(int n, int m, vector< vector<double> >& a,
            vector< vector<double> >& b) {
  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv;
  int *indxc, *indxr, *ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];
  irow = icol = 0;  // to fool the compiler not to get a warning
  for (j = 0; j < n; j++) ipiv[j] = 0;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
    if (ipiv[j] != 1)
    for (k = 0; k < n; k++) {
      if (ipiv[k] == 0) {
        if (abs(a[j][k]) >= big) {
          big = abs(a[j][k]);
          irow = j;
          icol = k;
        }
      }
    }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l = 0; l < n; l++) swap(a[irow][l], a[icol][l]);
      for (l = 0; l < m; l++) swap(b[irow][l], b[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0) throw("gaussj: Singular Matrix");
    pivinv = 1.0 / a[icol][icol];
    a[icol][icol] = 1.0;
    for (l = 0; l< n; l++) a[icol][l] *= pivinv;
    for (l = 0; l< m; l++) b[icol][l] *= pivinv;
    for (ll = 0; ll< n; ll++)
    if (ll != icol) {
      dum = a[ll][icol];
      a[ll][icol] = 0.0;
      for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
      for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
    }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
    for (k = 0; k < n; k++)
      swap(a[k][indxr[l]], a[k][indxc[l]]);
  }
  delete[] indxc;
  delete[] indxr;
  delete[] ipiv;
}

}  // namespace

// Public
// *****************************************************************************

// BuildModel
// *****************************************************************************
void Model::BuildModel() {
  if (style().compare("mmm") == 0) {
    BuildAtomElement();
    BuildElementAtom();
    BuildElementCenterAtom();
    BuildShapeFunction();
    BuildScheme();    
    BuildForceWeight();
    BuildMass();
  } else if (style().compare("full_atomistic") == 0) {
    atom_mass_.set_all_to(mmm()->mass());
    atom_type_.set_all_to(2);
    is_rep_.set_all_to(true);
  }
  CountAtomTypeNum();
}

// InitMMM
// *****************************************************************************
void Model::InitMMM(string scheme) {
  atom_element_.set_size(atom_group()->atom_num(), kMaxAtomElementNum);
  atom_element_num_.set_1d_size(atom_group()->atom_num());
  atom_mass_.set_1d_size(atom_group()->atom_num());
  atom_type_.set_1d_size(atom_group()->atom_num());
  atom_type_num_.set_1d_size(5);
  element_atom_.set_size(mesh()->element_num(), kMaxElementAtomNum);
  element_atom_num_.set_1d_size(mesh()->element_num());
  element_center_atom_.set_1d_size(mesh()->element_num());
  force_weight_.set_1d_size(atom_group()->atom_num());
  is_rep_.set_1d_size(atom_group()->atom_num());
  shape_function_.set_size(atom_group()->atom_num(),
                           mesh()->element_node_num());
  set_scheme(scheme);
  set_style("mmm");
  set_is_full_atomistic(false);
}

// InitFullAtomistic
// *****************************************************************************
void Model::InitFullAtomistic() {
  atom_mass_.set_1d_size(atom_group()->atom_num());
  atom_type_.set_1d_size(atom_group()->atom_num());
  atom_type_num_.set_1d_size(5);
  is_rep_.set_1d_size(atom_group()->atom_num());
  set_style("full_atomistic");
  set_is_full_atomistic(true);
}

// SetAtomListType
// *****************************************************************************
// Divides tasks of BuildModel into 2 parts and performs setting of types of
// atom list between the two parts.
void Model::SetAtomListType(const vector<int>& atom_list, int input_type) {
  // build model - part 1
  BuildAtomElement();
  BuildElementAtom();
  BuildElementCenterAtom();
  BuildShapeFunction();
  BuildScheme();
  // set atom type
  for (int i = 0; i < atom_list.size(); i++) {
    atom_type_.set_1d_element(atom_list[i], input_type);
  }
  // type constraints
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    // ssmp be inside an element
    if (atom_type_.get_1d_element(i) == 4 &&
        atom_element_num_.get_1d_element(i) == 0) {
      output()->ThrowError("atom" + output()->ToString(i) +
                           " is set to type ssmp, but it is idle");
    }
    // nirep inside an element turn others into nirep
    if (atom_type_.get_1d_element(i) == 2) {
      if (atom_element_num_.get_1d_element(i) > 0) {
        int element = atom_element_.get_element(i, 0);
        for (int j = 0; j < element_atom_num_.get_1d_element(element); j++) {
          if (atom_type_.get_1d_element(element_atom_.get_element(element,
                                                                  j)) == 1) {
            continue;
          }
          atom_type_.set_1d_element(element_atom_.get_element(element, j), 2);
        }
      }
    }
  }
  // is rep
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (atom_type_.get_1d_element(i) <= 2) is_rep_.set_1d_element(i, true);
  }
  // build model - part 2
  BuildForceWeight();
  BuildMass();
  CountAtomTypeNum();
}

// Private
// *****************************************************************************

// BuildAtomElement
// *****************************************************************************
void Model::BuildAtomElement() {
  if (mmm()->dimension() == 1) {
    BuildAtomElementLine();
  } else if (mmm()->dimension() == 2) {
    BuildAtomElementTriangle();
  } else if (mmm()->dimension() == 3) {
    BuildAtomElementTetrahedron();
  }
}

// BuildAtomElementLine
// *****************************************************************************
// Assigns an element to an atom if the atom is inside the element (including
// boundary). An atom is inside if it is simply between the x-coordinates of the
// nodes of the element.
void Model::BuildAtomElementLine() {
  vector<int> atom_element;
  vector<int> atom_element_num;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    atom_element_num.push_back(0);
    for (int j = 0; j < mesh()->element_num(); j++) {
      int element_node_1 = mesh()->element_node_at(j, 0);
      int element_node_2 = mesh()->element_node_at(j, 1);
      if (atom_group()->position_.get_element(i, 0) >=
              atom_group()->position_.get_element(element_node_1, 0) &&
              atom_group()->position_.get_element(i, 0) <=
              atom_group()->position_.get_element(element_node_2, 0)) {
        atom_element.push_back(j);
        atom_element_num[i]++;
      }
    }
    atom_element.insert(atom_element.end(),
                        kMaxAtomElementNum - atom_element_num[i], -1);
  }
  atom_element_.set(atom_element);
  atom_element_num_.set(atom_element_num);
}

// BuildAtomElementTetrahedron
// *****************************************************************************
// Assigns an element to an atom if the atom is inside the element (including
// boundary). For further details, see
// http://steve.hollasch.net/cgindex/geometry/ptintet.html.
void Model::BuildAtomElementTetrahedron() {
  vector<int> atom_element;
  vector<int> atom_element_num;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    atom_element_num.push_back(0);
    for (int j = 0; j < mesh()->element_num(); j++) {
      int sign[5] = {0, 0, 0, 0, 0};
      double matrix_1[4][4];
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 3; l++) {
          matrix_1[k][l] =
              atom_group()->position_.get_element(mesh()->element_node_at(j, k),
                                                  l);
        }
        matrix_1[k][3] = 1.0;
      }
      double matrix_1_determinant = determinant44(matrix_1);
      if (abs(matrix_1_determinant) < 1.0e-6) matrix_1_determinant = 0.0;
      sign[0] = matrix_1_determinant > 0 ? 1 : -1;
      for (int k = 0; k < 4; k++) {
        double matrix_2[4][4];
        for (int l = 0; l < 4; l++) {
          for (int m = 0; m < 4; m++) {
            matrix_2[l][m] = matrix_1[l][m];
          }
        }
        for (int l = 0; l < 3; l++) {
          matrix_2[k][l] = atom_group()->position_.get_element(i, l);
        }
        double matrix_2_determinant = determinant44(matrix_2);
        // Include few exterior points
        if (abs(matrix_2_determinant) < 1.0e-2) matrix_2_determinant = 0.0;
        sign[k + 1] = matrix_2_determinant > 0 ? 1 : -1;
      }
      bool is_positive = true;
      bool is_negative = true;
      for (int k = 0; k < 5; k++) {
        if (sign[k] < 0) {
          is_positive = false;
        } else {  // Case = 0 is included here
          is_negative = false;
        }
      }
      if (is_positive || is_negative) {
        atom_element.push_back(j);
        atom_element_num[i]++;
        break;
      }
    }
    if (atom_element_num[i] > kMaxAtomElementNum) {
      output()->ThrowError("atom element number exceededs allowed maximum");
    }
    atom_element.insert(atom_element.end(),
                        kMaxAtomElementNum - atom_element_num[i], -1);    
  }
  atom_element_.set(atom_element);
  atom_element_num_.set(atom_element_num);
}

// BuildAtomElementTriangle
// *****************************************************************************
// Assigns an element to an atom if the atom is inside the element (including
// boundary). An atom is inside if the angles it makes with every pair of
// element nodes sums up to 2PI.
void Model::BuildAtomElementTriangle() {
  vector<int> atom_element;
  vector<int> atom_element_num;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    atom_element_num.push_back(0);
    for (int j = 0; j < mesh()->element_num(); j++) {
      double triangle_coordinate[3][2];
      for (int k = 0; k < mesh()->element_node_num(); k++) {
        triangle_coordinate[k][0] =
            atom_group()->position_.get_element(mesh()->element_node_at(j, k),
                                                0);
        triangle_coordinate[k][1] =
            atom_group()->position_.get_element(mesh()->element_node_at(j, k),
                                                1);
      }
      bool is_magnitude_0 = false;
      double magnitude[3];  // magnitude of the vector below
      double vector[3][2];  // vector between atom and nodes of the element
      for (int k = 0; k < mesh()->element_node_num(); k++) {
        vector[k][0] = triangle_coordinate[k][0] -
                           atom_group()->position_.get_element(i, 0);
        vector[k][1] = triangle_coordinate[k][1] -
                           atom_group()->position_.get_element(i, 1);
        magnitude[k] = sqrt(vector[k][0] * vector[k][0] + vector[k][1] *
            vector[k][1]);
        if (magnitude[k] < 0.1) is_magnitude_0 = true;
      }
      double angle[3];
      int angle_count = 0;
      for (int k = 0; k < mesh()->element_node_num() - 1; k++) {
        for (int l = k + 1; l < mesh()->element_node_num(); l++) {          
          double cos_value = (vector[k][0] * vector[l][0] + vector[k][1] *
              vector[l][1]) / (magnitude[k] * magnitude[l]);
          if (cos_value < -1 && (cos_value + 1) < 1e-3) cos_value = -1.0;
          if (cos_value > 1 && (cos_value - 1) < 1e-3) cos_value = 1.0;
          angle[angle_count] = acos(cos_value);
          angle_count++;
        }
      }
      double angle_sum = angle[0] + angle[1] + angle[2];
      if (abs((angle_sum - 2 * kPi)) < 0.1 || is_magnitude_0) {
        atom_element.push_back(j);
        atom_element_num[i]++;
        if (!mesh()->is_node_ID(i)) break;
      }
    }
    atom_element.insert(atom_element.end(),
                        kMaxAtomElementNum - atom_element_num[i], -1);
  }
  atom_element_.set(atom_element);
  atom_element_num_.set(atom_element_num);
}

// BuildElementAtom
// *****************************************************************************
// Constructs from atom_element_.
void Model::BuildElementAtom() {
  vector<int> element_atom;
  vector<int> element_atom_num(mesh()->element_num());
  for (int i = 0; i < mesh()->element_num(); i++) {
    element_atom_num[i] = 0;
    for (int j = 0; j < atom_group()->atom_num(); j++) {
      for (int k = 0; k < atom_element_num_.get_1d_element(j); k++) {
        if (i == atom_element_.get_element(j, k)) {
          element_atom.push_back(j);
          element_atom_num[i]++;
        }
      }
    }
    if (element_atom_num[i] > kMaxElementAtomNum) {
      output()->ThrowError("element atom number exceeds allowed maximum");
    }
    element_atom.insert(element_atom.end(),
                        kMaxElementAtomNum - element_atom_num[i], -1);
  }
  element_atom_.set(element_atom);
  element_atom_num_.set(element_atom_num);
}

// BuildElementCenterAtom
// *****************************************************************************
// Selects atom closest to the element center. Distances are measure in square.
void Model::BuildElementCenterAtom() {
  vector<int> element_center_atom;
  for (int i = 0; i < mesh()->element_num(); i++) {
    double element_center_position[3] = {0, 0, 0};
    for (int j = 0; j < mesh()->element_node_num(); j++) {
      for (int k = 0; k < mmm()->dimension(); k++) {
        element_center_position[k] +=
            atom_group()->position_.get_element(mesh()->element_node_at(i, j),
            k) / static_cast<double>(mesh()->element_node_num());
      }
    }
    double minimum_distance = 1.0e6;
    int closest_atom = -1;    
    for (int j = 0; j < element_atom_num_.get_1d_element(i); j++) {
      int atom = element_atom_.get_element(i, j);
      double distance = 0;
      for (int k = 0; k < mmm()->dimension(); k++) {
        distance +=
            pow(atom_group()->position_.get_element(atom, k) -
            element_center_position[k], 2);
      }
      if (distance < minimum_distance) {
        minimum_distance = distance;
        closest_atom = atom;
      }
    }
    element_center_atom.push_back(closest_atom);
  }
  element_center_atom_.set(element_center_atom);
}

// BuildForceWeight
// *****************************************************************************
void Model::BuildForceWeight() {
  vector<double> force_weight;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    double weight = 1.0;  // if irep, nirep, or ssmp
    if (atom_type_.get_1d_element(i) == 3) {  // if psmp
      weight = 1.0;
      int element = atom_element_.get_element(i, 0);
      for (int j = 0; j < element_atom_num_.get_1d_element(element); j++) {
        if (atom_type_.get_1d_element(element_atom_.get_element(element, j)) ==
            5) {
          weight++;
        }
      }
    }
    if (atom_type_.get_1d_element(i) == 5) weight = 0.0;  // if nsmp
    force_weight.push_back(weight);
  }
  force_weight_.set(force_weight);
}

// BuildMass
// *****************************************************************************
void Model::BuildMass() {
  double temp_mass;
  vector<double> atom_mass;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (atom_type_.get_1d_element(i) == 1) {  // if irep
      temp_mass = mmm()->mass();
      for (int j = 0; j < atom_element_num_.get_1d_element(i); j++) {
        int element = atom_element_.get_element(i, j);
        int ghost_atom_num = 0;
        for (int k = 0; k < element_atom_num_.get_1d_element(element); k++) {
          if (!is_rep_.get_1d_element(element_atom_.get_element(element, k))) {
            ghost_atom_num++;
          }
        }
        double ghost_mass = static_cast<double>(ghost_atom_num) * mmm()->mass();
        temp_mass += ghost_mass /
            static_cast<double>(mesh()->element_node_num());
      }
    } else if (atom_type_.get_1d_element(i) == 2) {  // if nirep
      temp_mass = mmm()->mass();
    } else {  // if psmp, ssmp, or nsmp
      temp_mass = 0.0;
    }
    atom_mass.push_back(temp_mass);
  }
  atom_mass_.set(atom_mass);
}

// BuildScheme
// *****************************************************************************
void Model::BuildScheme() {
  int type;
  vector<int> atom_type;
  if (scheme_.compare("no_ssmp") == 0) {
    for (int i = 0; i < atom_group()->atom_num(); i++) {
      type = 5;  // nsmp if nothing else
      for (int j = 0; j < mesh()->element_num(); j++) {
        for (int k = 0; k < mesh()->element_node_num(); k++) {
          if (i == mesh()->element_node_at(j, k)) type = 1;  // irep if node
        }
      }
      for (int j = 0; j < mesh()->element_num(); j++) {
        if (i == element_center_atom_.get_1d_element(j)) {
          type = 3;  // psmp if center atom
        }
      }
      if (atom_element_num_.get_1d_element(i) == 0) type = 2;  // nirep if idle
      atom_type.push_back(type);
    }
  } else if (scheme_.compare("all_ssmp") == 0) {
    for (int i = 0; i < atom_group()->atom_num(); i++) {
      type = 4;  // ssmp if nothing else
      for (int j = 0; j < mesh()->element_num(); j++) {
        for (int k = 0; k < mesh()->element_node_num(); k++) {
          if (i == mesh()->element_node_at(j, k)) type = 1;  // irep if node
        }
      }
      for (int j = 0; j < mesh()->element_num(); j++) {
        if (i == element_center_atom_.get_1d_element(j)) {
          type = 3;  // psmp if center atom
        }
      }
      if (atom_element_num_.get_1d_element(i) == 0) type = 2;  // nirep if idle
      atom_type.push_back(type);
    }
  } else if (scheme_.compare("ssmp_around_irep") == 0) {
    for (int i = 0; i < atom_group()->atom_num(); i++) {
      type = 5;  // nsmp if nothing else
      for (int j = 0; j < mesh()->element_num(); j++) {
        for (int k = 0; k < mesh()->element_node_num(); k++) {
          if (i == mesh()->element_node_at(j, k)) type = 1;  // irep if node
        }
      }
      for (int j = 0; j < mesh()->element_num(); j++) {
        if (i == element_center_atom_.get_1d_element(j)) {
          type = 3;  // psmp if center atom
        }
      }
      if (atom_element_num_.get_1d_element(i) == 0) type = 2;  // nirep if idle
      if (type == 5) {
        for (int j = 0; j < mesh()->element_num(); j++) {
          for (int k = 0; k < mesh()->element_node_num(); k++) {
            int node = mesh()->element_node_at(j, k);
            double distance_square = 0;
            for (int l = 0; l < mmm()->dimension(); l++) {
              double difference = atom_group()->position_.get_element(i, l) -
                                  atom_group()->position_.get_element(node, l);
              distance_square += difference * difference;
            }
            if (distance_square <= potential()->cutoff_radius_square()) {
              type = 4;  // ssmp if around irep
            }
          }
        }
      }
      atom_type.push_back(type);
    }
  } else if (scheme_.compare("ssmp_around_psmp") == 0) {
    for (int i = 0; i < atom_group()->atom_num(); i++) {
      type = 5;  // nsmp if nothing else
      for (int j = 0; j < mesh()->element_num(); j++) {
        for (int k = 0; k < mesh()->element_node_num(); k++) {
          if (i == mesh()->element_node_at(j, k)) type = 1;  // irep if node
        }
      }
      for (int j = 0; j < mesh()->element_num(); j++) {
        if (i == element_center_atom_.get_1d_element(j)) {
          type = 3;  // psmp if center atom
        }
      }
      if (atom_element_num_.get_1d_element(i) == 0) type = 2;  // nirep if idle
      if (type == 5) {
        for (int j = 0; j < mesh()->element_num(); j++) {
          int psmp = element_center_atom_.get_1d_element(j);
          double distance_square = 0;
          for (int k = 0; k < mmm()->dimension(); k++) {
            double difference = atom_group()->position_.get_element(i, k) -
                                atom_group()->position_.get_element(psmp, k);
            distance_square += difference * difference;
          }
          if (distance_square <= potential()->cutoff_radius_square()) {
            type = 4;  // ssmp if around psmp
          }          
        }
      }
      atom_type.push_back(type);
    }
  }
  atom_type_.set(atom_type);
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (atom_type_.get_1d_element(i) <= 2) is_rep_.set_1d_element(i, true);
  }
}

// BuildShapeFunction
// *****************************************************************************
// Constructs linear shape functions. For further details, see:
// http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch16.d/IFEM.Ch16.pdf.

void Model::BuildShapeFunction() {
  vector<double> shape_function;
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    if (atom_element_num_.get_1d_element(i) == 0) {
      shape_function.insert(shape_function.end(), mesh()->element_node_num(), 
                            0.0);
      continue;
    }
    int element = atom_element_.get_element(i, 0);
    vector< vector<double> > A(mesh()->element_node_num(),
                               vector<double>(mesh()->element_node_num()));
    for (int j = 0; j < mesh()->element_node_num(); j++) {
      int node = mesh()->element_node_at(element, j);      
      A[0][j] = 1.0;
      for (int k = 1; k < mmm()->dimension() + 1; k++) {
        A[k][j] = atom_group()->position_.get_element(node, k - 1);
      }
    }
    vector< vector<double> > b(mesh()->element_node_num(), vector<double>(1));
    b[0][0] = 1.0;
    for (int j = 1; j < mesh()->element_node_num(); j++) {
      b[j][0] = atom_group()->position_.get_element(i, j - 1);
    }
    gaussj(mesh()->element_node_num(), 1, A, b);  // solve for A x = b
    for (int j = 0; j < mesh()->element_node_num(); j++) {
      if (abs(b[j][0]) < 1e-6) b[j][0] = 0.0;
      shape_function.push_back(b[j][0]);
    }
  }
  shape_function_.set(shape_function);
}

// CountAtomTypeNum
// *****************************************************************************
void Model::CountAtomTypeNum() {
  atom_type_num_.reset();
  for (int i = 0; i < atom_group()->atom_num(); i++) {
    int type = model()->atom_type_.get_1d_element(i);
    atom_type_num_.add_1d_element(type - 1, 1);
  }
}

