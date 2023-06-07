/*
 *  Main header file for percolation exercise.
 */

 /*
  *  System size L
  */

#define L 768

/*
 *  Cell density rho (between 0 and 1)
 */

#define RHO 0.49

#define MASTER_RANK 0

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void cellwrite(char* cellfile, int cell[L][L]);
void cellwritedynamic(char* cellfile, int** cell, int l);

/*
  *  Random numbers
 */  

void rinit(int ijkl);
float uni(void);
