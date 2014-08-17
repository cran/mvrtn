#include<stdio.h>
#include<math.h>
#include"MVTN.h"
#include"R.h"

void MVTN_Left (double *zmu, double *zsig, double *c, double * mvtn_left)
{

        // Computes mean of left truncated normal distribution
        mvtn_left[0] = (sqrt(2/M_PI)**zsig + exp((pow((*c - *zmu),2))/(2*(pow(*zsig,2))))**zmu*erfc(((*c - *zmu)/(sqrt(2)**zsig))))/ (exp( pow((*c - *zmu),2) / (2*(pow(*zsig,2))) )*(1 + erf(((-*c + *zmu)/(sqrt(2)**zsig)) )  ));


        // Computes variance of left truncated normal distribution
        mvtn_left[1] = (0.22507907903927651**zsig*erfc((0.7071067811865475*(*c - 1.**zmu))/ (*zsig))*(3.5449077018110318*pow(2.718281828459045,(0.5*pow(*c + *zmu,2))/pow(*zsig,2))*(*c - 1.**zmu)*erfc((0.7071067811865475*(*c - 1.**zmu))/(*zsig)) - 1.4142135623730951**zsig*(2.*pow(2.718281828459045,(2.**c**zmu)/pow(*zsig,2)) - 3.141592653589793*pow(2.718281828459045,(pow(*c,2) + pow(*zmu,2))/pow(*zsig,2))*pow(erfc((0.7071067811865475*(*c - 1.**zmu))/(*zsig)),2))))/(pow(2.718281828459045,(1.*(pow(*c,2) + pow(*zmu,2)))/pow(*zsig,2))*pow(1. + erf((0.7071067811865475*(-1.**c + *zmu))/(*zsig)),3));
}

