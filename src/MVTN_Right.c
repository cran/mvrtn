#include<stdio.h>
#include<math.h>
#include"MVTN.h"
#include"R.h"


void MVTN_Right (double *zmu, double *zsig, double *c, double * mvtn_right)
{
        // Computes mean right truncation of normal distribution
        mvtn_right[0] = ((*zmu - (sqrt(2/M_PI)**zsig)/exp(pow((*c - *zmu),2)/(2*(pow(*zsig,2)))) + (*zmu)*erf(((*c - *zmu)/(sqrt(2)**zsig))))/erfc(((-*c + *zmu)/(sqrt(2)**zsig)) ) );

        // Calculates right truncated variance of normal distribution
        mvtn_right[1] = (0.22507907903927651**zsig*(1. + erf((0.7071067811865475*(*c - 1.**zmu))/(*zsig)))*(1.4142135623730951**zsig*(-2.*pow(2.718281828459045,(2.**c**zmu)/pow(*zsig,2)) + 3.141592653589793*pow(2.718281828459045,(pow(*c,2) + pow(*zmu,2))/pow(*zsig,2))*pow(1. + erf((0.7071067811865475*(*c - 1.**zmu))/(*zsig)),2)) + 3.5449077018110318*pow(2.718281828459045,(0.5*pow(*c + *zmu,2))/pow(*zsig,2))*(*c - 1.**zmu)*(-2. + erfc((0.7071067811865475*(*c - 1.**zmu))/(*zsig)))))/(pow(2.718281828459045,(1.*(pow(*c,2) + pow(*zmu,2)))/pow(*zsig,2))*pow(erfc((0.7071067811865475*(-1.**c + *zmu))/(*zsig)),3));
}


