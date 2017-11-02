#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <complex.h> 
#include <stdbool.h>
#include <math.h>

void main()
{
  double GasDiskScaleLength =4.899684e-5;
  double ColdGas =1.508247e-5;
  
  double new_rad=2.905167;
  double new_mass =1.110189e-11;
  
  double StellarDiskScaleLength =1;
  double StellarMass =1;
  
  int gas = 1;

  if ((gas==0)&((StellarDiskScaleLength==0)|(StellarMass==0)))
  { 
    printf("%e - 1\n",new_rad);
  }
  else if ((gas==1)&((GasDiskScaleLength==0)|(ColdGas==0)))
  { 
    printf("%e -2 \n",new_rad);
  }
  else if ((gas==0)&(new_rad<1e-6))
    printf("%e -3 \n",StellarDiskScaleLength);
  else if ((gas==1)&(new_rad<1e-6))
    printf("%e -4 \n",GasDiskScaleLength);
  else if ((gas==0)&(StellarDiskScaleLength==new_rad))
    printf("%e -3 \n",StellarDiskScaleLength);
  else if ((gas==1)&(GasDiskScaleLength==new_rad))
    printf("%e -4 \n",GasDiskScaleLength);
  else
  {
    double aa=0;
    double rad=0;
    if (gas==0)
    {
      aa = 0.159155*StellarMass/pow(StellarDiskScaleLength,2);
      rad = StellarDiskScaleLength;
    }
    else
    {
      aa = 0.159155*ColdGas/pow(GasDiskScaleLength,2);
      rad = GasDiskScaleLength;
    }
    double cc = 0.159155*new_mass/pow(new_rad,2);
    if (aa==0)
    {
      printf("%e -5 \n",new_rad);
    }
    else if (cc==0)
    {
      printf("%e -6 \n",rad);
    }
    double yy[30];
    double rr[30]; 
    double sumrr = 0;
    double sumyy = 0;
    int i;
    double maxr=fmax(rad,new_rad);
    rr[0]=log(maxr)-1;
    double step = ((log(maxr)+1)/rr[0])/30;
    for (i = 0; i < 30; i++)
    {
      yy[i]=aa*exp(-exp(rr[i])/rad)+cc*exp(-exp(rr[i])/new_rad);
      sumrr+=exp(rr[i]);
      sumyy+=log(yy[i]);
      if (i<29)
        rr[i+1]=rr[i]+step;
    }
    double meanr=sumrr/30;
    double meany=sumyy/30;
    double rr_meanr_yy_meany=0;
    double rr_meanr_sq=0;
    
    for (i = 0; i < 30; i++)
    {
      rr_meanr_yy_meany+=(exp(rr[i])-meanr)*(log(yy[i])-meany);
      rr_meanr_sq+=pow(exp(rr[i])-meanr,2);
    }
    printf("SCALE RADIUS = %e, -7\n",-1/(rr_meanr_yy_meany/rr_meanr_sq));
  }
}

