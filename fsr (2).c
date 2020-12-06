#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include "fsr.h"

struct _matrix {
  size_t v_size;
  size_t h_size;
  double** arr;
  int basis_size;
};

void print_matrix(Matrix *matrix) {
  for (size_t i = 0; i < matrix->v_size; ++i) {
    for (size_t j = 0; j < matrix->h_size; ++j) {
      printf("%10.3llf ", matrix->arr[i][j]);
    }
    printf("\n");
  }

  printf("\n");
}

Matrix* create_matrix(size_t v_size, size_t h_size) {
  Matrix* new_matrix = NULL;

  new_matrix = (Matrix*) malloc(sizeof(Matrix));
  new_matrix->h_size = h_size;
  new_matrix->v_size = v_size;
  new_matrix->basis_size = h_size;
  new_matrix->arr = (double**) malloc(v_size * sizeof(double*));

  for (size_t i = 0; i < v_size; ++i) {
    new_matrix->arr[i] = (double*) malloc(h_size * sizeof(double));
  }

  return new_matrix;
}

void insert_value_matrix(Matrix* matrix, size_t row, size_t col, double val) {
  if (row >= matrix->v_size || col >= matrix->h_size) {
    fprintf(stderr, 
            "Range error in func: insert_value_matrix\n Row: %d of %d. Column: %d of %d\n",
            row, matrix->v_size - 1, col, matrix->h_size - 1);
    free_matrix(matrix);
    exit(139);
  }

  matrix->arr[row][col] = val;
}

double show_value_matrix(Matrix *matrix, size_t row, size_t col) {
  if (row >= matrix->v_size || col >= matrix->h_size) {
    fprintf(stderr,
            "Range error in func: show_value_matrix\n Row: %d of %d. Column: %d of %d\n",
            row, matrix->v_size - 1, col, matrix->h_size - 1);
    free_matrix(matrix);
    exit(139);
  }

  return matrix->arr[row][col];
}

void swap_row(double* left, double* right, size_t length) {
  double temp;

  for (size_t i = 0; i < length; ++i) {
    temp = left[i];
    left[i] = right[i];
    right[i] = temp;
  }
}

void summ_row(double* left, double* right, double coef, size_t length) {
  for (size_t i = 0; i < length; ++i) {
    left[i] += right[i] * coef;
  }
}

int is_zero(double *arr, size_t length) {
  for (size_t i = 0; i < length; ++i) {
    if (arr[i] != 0.0) {
      return 0;
    }
  }

  return 1;
}

size_t count_free_size(Matrix *matrix) {
  double **arr = matrix->arr;
  size_t res = matrix->v_size;

  for (size_t i = 0; i < matrix->v_size; ++i) {
    if (is_zero(arr[i], matrix->h_size)) {
      --res;
    }
  }

  return res;
}

void calc_fsr(Matrix* matrix) {
  double** arr = matrix->arr;
  size_t free_size = count_free_size(matrix);
  matrix->basis_size = matrix->h_size - free_size;

  // move all zero string to end
  for (size_t i = 0; i < matrix->v_size; ++i) {
    if (is_zero(arr[i], matrix->h_size)) {
      for (size_t j = i; j < matrix->v_size - 1; ++j) {
        swap_row(arr[j], arr[j + 1], matrix->h_size);
      }
    }
  }

  // solve for free variables
  size_t basis_size = matrix->basis_size;
  for (size_t i = 0; i < free_size; ++i) {
    if (arr[i][i + basis_size] == 0.0) {
      for (size_t j = i + 1; j < matrix->v_size; ++j) {
        if (arr[j][i + basis_size] != 0.0) {
          swap_row(arr[i], arr[j], matrix->h_size);
          break;
        }
      }
      if (arr[i][i + basis_size] == 0.0) continue;
    }

    double temp = arr[i][i + basis_size];
    for (size_t j = 0; j < matrix->h_size; ++j) {
      arr[i][j] /= temp;
    }

    for (size_t j = 0; j < free_size; ++j) {
      if (j == i)
        continue;
      if (arr[j][i + basis_size] == 0.0)
        continue;

      double coef = (-1.0) * arr[j][i + basis_size];
      summ_row(arr[j], arr[i], coef, matrix->h_size);
    }
  }
}

void solve_matrix(Matrix* matrix) {
  double coef = 0.0;
  double** arr = matrix->arr;
  size_t max_size = matrix->v_size < matrix->h_size ? matrix->v_size : 
                                                      matrix->h_size;

  for (size_t i = 0; i < max_size; ++i) {
    if (arr[i][i] == 0.0) {
      // try to find with non-zero element and swap with current
      for (size_t j = i + 1; j < matrix->v_size; ++j) {
        if (arr[j][i] != 0.0) {
          swap_row(arr[i], arr[j], matrix->h_size);
          break;
        }
      }
      // no so such row was found, then just carry on
      if (arr[i][i] == 0.0) continue;
    }

    // normalize row
    double temp = arr[i][i];
    for (size_t j = 0; j < matrix->h_size; ++j) {
      arr[i][j] /= temp;
    }

    // there we make arr[i][{1..m}\{i}] equal to zero
    for (size_t j = 0; j < max_size; ++j) {
      if (j == i) continue;
      if (arr[j][i] == 0.0) continue;
      
      coef = (-1.0) * arr[j][i];
      summ_row(arr[j], arr[i], coef, matrix->h_size);
    }
  }

  calc_fsr(matrix);
}

void print_fsr(Matrix* matrix) {
  double** arr = matrix->arr;
  for (size_t i = 0; i < matrix->basis_size; ++i) {
    for (size_t j = 0; j < matrix->basis_size; ++j) {
      if (i == j) {
        printf("%10.3lf ", 1.0);
      } else {
        printf("%10.3lf ", 0.0);
      }
    }
    for (size_t j = 0; j < matrix->h_size - matrix->basis_size; ++j) {
      printf("%10.3lf ", -arr[j][i]);
    }
    printf("\n");
  }
}

void free_matrix(Matrix* matrix) { 
  for (size_t i = 0; i < matrix->v_size; ++i) {
    free(matrix->arr[i]);
  }

  free(matrix->arr);
  free(matrix);
}