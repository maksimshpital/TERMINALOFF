#ifndef _FSR_H
#define _FSR_H
#include <stdint.h>

typedef struct _matrix Matrix;  // Hide implementation=

//                h_size
//               (columns)
//          | a11 a12 ... a1n |
//          | a21 a22 ... a2n |
//  v_size  |     . . .       |
//  (rows)  |     . . .       |
//          |     . . .       |
//          | am1 am2 ... amn |
//
//  indexation in matrix: a[row][column]

Matrix* create_matrix(size_t v_size, size_t h_size);  

void    insert_value_matrix(Matrix*, size_t row, size_t col, double val); // Sorry son, it's not a C++
double  show_value_matrix(Matrix*, size_t row, size_t col); // so enjoy no operator overloading
void    solve_matrix(Matrix* matrix);
void    print_fsr(Matrix* matrix);
void    free_matrix(Matrix* matrix);

#endif