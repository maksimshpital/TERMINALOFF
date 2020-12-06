#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "fsr.h"

int main() {
  size_t h_size, v_size;
  Matrix* matrix = NULL;

  scanf_s("%zu %zu", &v_size, &h_size);
  matrix = create_matrix(v_size, h_size);

  double temp = 0.0;
  for (size_t i = 0; i < v_size; ++i) {
    for (size_t j = 0; j < h_size; ++j) {
      scanf_s("%lf", &temp);
      insert_value_matrix(matrix, i, j, temp);
    }
  }
  solve_matrix(matrix);
  print_fsr(matrix);

  free_matrix(matrix);
  return 0;
}