#include <stdio.h>
#include "mcce.h"

void welcome();

int main(int argc, char *argv[])
{
   welcome();

   /* STEP 0, initialization */
   db_open();
   printf("Step 0. Initialize enviroment\n"); fflush(stdout);
   if (init()) {
      db_close();
      printf("Help message: double check file \"run.prm\" in current directory.\n");
      return USERERR;
   }
   else printf("Step 0 Done.\n");

   /* STEP 1, premcce */
   if (env.do_premcce) {
      printf("Step 1. Test and format structral file\n"); fflush(stdout);
      if (premcce()) {db_close(); return USERERR;}
      else printf("Step 1 Done.\n");
   }
   else printf("Not doing \"Step 1. Test and format structral file\"\n");

   /* STEP 2. rotamers */
   if (env.do_rotamers) {
      printf("Step 2. Make multi side chain conformers\n"); fflush(stdout);
   if (env.sidechain_opt==1) {//genetic algorithm sidechain packing optimization
        printf("***Using Pascal Comte's (Brock University-Computer Science,2010) [Paper Ref.: ] Sidechain Packing GA & Bi-directional Evolutionary Sampling***\n");
        rotamers_GA(argc,argv);
      }
      else if (rotamers()) {
		db_close(); return USERERR;
      }
      else printf("Step 2 Done.\n");
   }
   else printf("Not doing \"Step 2. Make multi side chain conformers\"\n");

   /* STEP 3. energies */
   if (env.do_energies) {
      printf("Step 3. Compute energy lookup table\n"); fflush(stdout);
      if (energies()) {
         db_close();
         return USERERR;}
      else printf("Step 3 Done.\n");
   }
   else printf("Not doing \"Step 3. Compute energy lookup table\"\n");

   /* STEP 4. Monte Carlo */
   if (env.do_monte) {
      printf("Step 4. Monte Carlo Sampling\n"); fflush(stdout);
      if (!env.monte_adv_opt) {
         if (monte()) {
            db_close(); 
            return USERERR;}
         else printf("Step 4 Done.\n");
       }
       else {
           if (monte2()) {db_close(); return USERERR;}
           else printf("Step 4 Done.\n");
       }
   }
   else printf("Step 4: Not doing\n");

   db_close();
   return 0;
}

void welcome()
{  printf("<<<    MCCE Multi-Conformation Continuum Electrostatics    >>> \n");
   printf("Vers.: 2.5.1 with 'extra energy' titr in monte & monte2 2014/4  \n");
   printf("===============================================================\n");
   fflush(stdout);
   return;
}
