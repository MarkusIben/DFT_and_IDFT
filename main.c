#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIG_LENGTH 320

extern double InputSignal_f32_1kHz_15kHz[SIG_LENGTH];

//global variables
double Output_REX[SIG_LENGTH / 2];
double Output_IMX[SIG_LENGTH / 2];
double Output_MAG[SIG_LENGTH / 2];
double Output_IDFT[SIG_LENGTH];

double PI = 3.14159265359;

//prototypes
void calc_sig_dft(double *sig_src_arr, double *sig_dest_rex_arr, double *sig_dest_imx_arr, int sig_length);
void get_dft_output_mag(double *sig_dest_mag_arr);
void calc_idft(double *idft_out_arr, double *sig_src_rex_arr, double *sig_src_imx_arr, int idft_length);


//main-----------------------------------------------------------------------------------------------------
int main()
{
   calc_sig_dft((double *) &InputSignal_f32_1kHz_15kHz[0],
                (double *) &Output_REX[0],
                (double *) &Output_IMX[0],
                (int) SIG_LENGTH);

   get_dft_output_mag((double*) &Output_MAG[0]);

   calc_idft((double *) &Output_IDFT, (double *) &Output_REX, (double *) &Output_IMX, (int) SIG_LENGTH);


   //Write waveform and spectra into files --------------------------------------------------------------------------------
   FILE *input_sig_fptr, *output_rex_fptr, *output_imx_fptr, *output_mag_fptr, *output_idft_fptr;      //create 5 file pointers


   //open files -------------------------------------
   input_sig_fptr    = fopen("input_signal.dat", "w");   //open file for input signal
   output_rex_fptr   = fopen("output_rex.dat", "w");     //open file for real part of spectrum
   output_imx_fptr   = fopen("output_imx.dat", "w");     //open file for imaginary part of spectrum
   output_mag_fptr   = fopen("output_mag.dat", "w");     //open file for magnitude of spectrum
   output_idft_fptr  = fopen("output_idft.dat", "w");    //open file for IDFT signal


   //write array content into files -------------------
   for (int i = 0; i < SIG_LENGTH; i++)                  //arrays with length "SIG_LENGTH"
   {
      fprintf(input_sig_fptr, "\n%f", InputSignal_f32_1kHz_15kHz[i]);   //write input signal
      fprintf(output_idft_fptr, "\n%f", Output_IDFT[i]);                //write idft output
   }
   for (int i = 0; i < (SIG_LENGTH / 2); i++)            //arrays with length "SIG_LENGTH" / 2
   {
      fprintf(output_rex_fptr, "\n%f", Output_REX[i]);   //write real part of spectrum (rex)
      fprintf(output_imx_fptr, "\n%f", Output_IMX[i]);   //write imaginary part of spectrum (imx)
      fprintf(output_mag_fptr, "\n%f", Output_MAG[i]);   //write magnitude of spectrum (mag)
   }


   //close files -------------------------------------
   fclose(input_sig_fptr);    //close file
   fclose(output_rex_fptr);   //close file
   fclose(output_imx_fptr);   //close file
   fclose(output_mag_fptr);   //close file
   fclose(output_idft_fptr);   //close file



    //printf("Hello world!\n");
    return 0;
}

//functions --------------------------------------------------------------------------------------------------
void calc_sig_dft(double *sig_src_arr, double *sig_dest_rex_arr, double *sig_dest_imx_arr, int sig_length)
{
   for (int j = 0; j < (sig_length / 2); j++)
   {
      sig_dest_rex_arr[j] = 0;
      sig_dest_imx_arr[j] = 0;
   }

   for (int k = 0; k < (sig_length / 2); k++)
   {
      for (int i = 0; i < sig_length; i++)
      {
         sig_dest_rex_arr[k] = sig_dest_rex_arr[k] + sig_src_arr[i] * cos(2 * PI * k * i / sig_length); //real part of spectrum
         sig_dest_imx_arr[k] = sig_dest_imx_arr[k] - sig_src_arr[i] * sin(2 * PI * k * i / sig_length); //imaginary part of spectrum
      }
   }
}


void get_dft_output_mag(double *sig_dest_mag_arr)
{
   for(int i = 0; i < (SIG_LENGTH / 2); i++)
   {
      sig_dest_mag_arr[i] = sqrt(pow(Output_REX[i], 2) + pow(Output_IMX[i], 2));
   }
}


void calc_idft(double *idft_out_arr, double *sig_src_rex_arr, double *sig_src_imx_arr, int idft_length)
{
   for (int k = 0; k < (idft_length / 2); k++)
   {
      sig_src_rex_arr[k] = sig_src_rex_arr[k] / (idft_length / 2);
      sig_src_imx_arr[k] = -sig_src_imx_arr[k] / (idft_length / 2);
   }
   sig_src_rex_arr[0] = sig_src_rex_arr[0] / 2;
   sig_src_imx_arr[0] = -sig_src_imx_arr[0] / 2;

   for(int i = 0; i < idft_length; i++)
   {
      idft_out_arr[i] = 0;
   }

   for(int k = 0; k < (idft_length / 2); k++)
   {
      for(int i = 0; i < idft_length; i++)
      {
         idft_out_arr[i] = idft_out_arr[i] + sig_src_rex_arr[k] * cos(2 * PI * k * i / idft_length);
         idft_out_arr[i] = idft_out_arr[i] + sig_src_imx_arr[k] * sin(2 * PI * k * i / idft_length);
      }
   }
}


