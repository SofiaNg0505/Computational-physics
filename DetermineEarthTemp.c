#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "lapacke.h"
#include <math.h>
#include <time.h>

int main(){
    int no_Pillars = 10000;
    double heights_pillars[no_Pillars];
    double density[no_Pillars];
    double preassure[no_Pillars];
    double h = 12000/no_Pillars;
    double M = 0.02896; //moon: 0.044([kg/mol]//([kg/mol]
    double g = 9.807;//3.721 
    double T_0= 288.15; //210 Kelvin 
    double R = 8.3143; // J/mol/K
    //mars
    //double T_0 = 210;
    //double  M =0.044;
   // double g = 3.721 ;
   //double R = 8.3143; // J/mol/K
    
    double absorption[no_Pillars];
    double temperature[no_Pillars]; 

    double T_inV[no_Pillars];
    double T_inR[no_Pillars];
    double T_out[no_Pillars];

    double sigma = 5.67*pow(10,-8); 

    FILE *albed=fopen("albedo", "w");
    FILE *t_earth=fopen("t_of_earth", "w");
    double albedo[20];

    for (int al = 0;al<20; al++) //mars = 0.250
    {
        albedo[al] = 0.25+al*0.01;
        printf("%lf\n",albedo[al]);
        fprintf(albed,"%lf\n",albedo[al]);



    double flux_earth = 344;
    double flux_mars = 344*0.431; 
    double flux_0 = flux_earth*(1-albedo[al]); //[]
    printf("flux %lf\n",flux_0);

    double sigma_vis  = 3*pow(10,-5);
    double sigma_ir = 3*pow(10,-4);

    for (int i =0 ; i<no_Pillars; i++){
        heights_pillars[i] = i*h;
        //Transmission and emission is 0 in the beginning
     //   printf("Height for pillar %d: %lf\n", i,heights_pillars[i]);
    }   
    preassure[0] = 101325; //mars = 655.00194 [Pascal]/[N/m][kg/s^2m]
  
    for (int i =0 ; i<no_Pillars; i++){ 
        preassure[i] = preassure[0]*exp((-M*g/(T_0*R))*heights_pillars[i]); //BAROMETRIC FORMULA
  
    }
    
    double H_n = R*T_0/(M*g); // Height scale. 
    density[0] = preassure[0] * M/(R*T_0);
    FILE *den=fopen("density", "w");
    FILE *pre=fopen("preassure", "w");
    FILE *temp=fopen("temperature", "w");
    FILE *height=fopen("height", "w");
    FILE *test=fopen("test", "w");
    for (int i =0 ; i<no_Pillars; i++){
        density[i] = density[0] * exp(-heights_pillars[i]/(H_n));
        temperature[i] = preassure[i]* M /(R*density[i]);
        //printf("Density for pillar %d: %lf [g/L]\n",i, density[i]);
        fprintf(test,"%lf\n",density[i]);
        fprintf(den,"%lf\n", density[i]);
        fprintf(pre,"%lf\n", preassure[i]);
        fprintf(temp,"%lf\n", temperature[i]);
        fprintf(height,"%lf\n", heights_pillars[i]*h);

    }
    fclose (den);
    fclose (test);
    fclose (pre);
    fclose (temp);
    fclose (height);
 
    T_inV[no_Pillars-1] = flux_0;
    absorption[no_Pillars-1] = flux_0*(1-exp(-sigma_vis*density[no_Pillars-1]*h));
    for (int pillar = no_Pillars-2; pillar>=0;pillar --){
        T_inV[pillar] = T_inV[pillar+1]*exp(-density[pillar]*h*sigma_vis);
        absorption[pillar] = T_inV[pillar+1]*(1-exp(-density[pillar]*h*sigma_vis));

    }
    absorption[0] =  T_inV[0];

    for (int i =0 ; i<no_Pillars; i++){
            T_inR[i] =0;
        }
  
    int iterations = 10000;
    for (int a = 0; a <iterations; a++){
        //earth surface-pillar
        T_inV[0] = T_inV[1]*exp(-density[0]*h*sigma_vis); //incoming visible light transmitted from cell above
        T_inR[0] =(T_inR[1]+absorption[1]*0.5) * exp(-density[0]*sigma_ir*h); //incoming infrared transmitted from cell above which also gets half of the emitted long wave
        absorption[0] = T_inV[1] + (0.5*absorption[1] + T_inR[1]); // absorb visible light and infrared from above and half of the emitted long wave. 
        T_out[0] = 0 ; //no incoming long wave from earth

        T_inV[1] = T_inV[2]*exp(-density[1]*h*sigma_vis); 
        T_inR[1] =(T_inR[2]+absorption[2]*0.5)*exp(-density[1]* h*sigma_ir);
        absorption[1] = (1-exp(density[1]*sigma_vis*h))*T_inV[2] + (1-exp(-sigma_ir*density[1]*h)) * (0.5*(absorption[2])+absorption[0] + T_inR[2]);
        T_out[1] = absorption[0]; // 100% reflection

    for (int pillar = no_Pillars-2; pillar>1;pillar --){
        T_inV[pillar] = T_inV[pillar+1]*exp(-density[pillar]*h*sigma_vis);
        T_inR[pillar] = (T_inR[pillar+1] + absorption[pillar+1]*0.5)*exp(-sigma_ir*density[pillar]*h);
        T_out[pillar] = (T_out[pillar-1] + absorption[pillar-1]*0.5)*exp(-sigma_ir*density[pillar]*h);
        absorption[pillar] =(1-exp(-sigma_ir*h*density[pillar]))* (T_out[pillar-1] +T_inR[pillar + 1] + 0.5*(absorption[pillar+1] + absorption[pillar-1]) )  + (1-exp(-sigma_vis*h*density[pillar]))*T_inV[pillar+1] ; 
        }   
        
        T_inV[no_Pillars-1] = flux_0; //at the surface we get 100% of the flux from sun
        T_inR[no_Pillars-1] = 0; //No IR coming in from above the atmosphere
        T_out[no_Pillars-1] = (T_out[no_Pillars-2] + 0.5*absorption[no_Pillars-2])*exp(-density[no_Pillars-1]*h*sigma_ir);
        absorption[no_Pillars-1] = flux_0*(1-exp(-sigma_vis*density[no_Pillars-1]*h)) + (0.5*absorption[no_Pillars-2] + T_out[no_Pillars-2])*(1-exp(-sigma_ir*density[no_Pillars-1]*h));

        }

        
        double F,Temp;
        F = absorption[0];
        Temp = sqrt(sqrt(F/sigma));
        printf(" %lf\n", Temp-274.15);
        fprintf(t_earth," %lf\n", Temp-274.15);

        FILE *absorb=fopen("absorption", "w");
        FILE *vis=fopen("Tvisible", "w");
        FILE *ir=fopen("Tir", "w");
        FILE *out=fopen("Tout", "w");
        for (int i = 0; i<no_Pillars; i ++){
            fprintf(absorb,"%lf\n", absorption[i]);
            fprintf(vis,"%lf\n", T_inV[i]);
            fprintf(ir,"%lf\n", T_inR[i]);
            fprintf(out,"%lf\n", T_out[i]);
        }
        
        fclose (absorb);
        fclose (vis);
        fclose (ir);
        fclose (out);
        flux_0 = 0;
        F=0;
        Temp = 0;
        for (int ok = 0;ok<no_Pillars; ok++){
            absorption[ok] = 0;
            density[ok] = 0;
            T_inR[ok] = 0;
            T_inV[ok] = 0;
            T_out[ok] = 0;
         }

    }
    fclose (albed);
    fclose (t_earth);
    
    }
    



   



