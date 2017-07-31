#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "mcce.hpp"

using namespace std;

/* check for tabs
if (line.find("\t") != string::npos) {
    cout << "\n   Tab charater is not recognized by MCCE\n";
    cout << "   Please replace tabs with spaces in file " << param_file.pszFileName << endl;
    return USERERR;
}
*/

#define LEN_KEY1 9
#define LEN_KEY2 6
#define LEN_KEY3 5
#define LEN_KEY  LEN_KEY1+LEN_KEY2+LEN_KEY3
Param& Param::load_param(ifstream& param_file)
{
    string line;
    string key1;                 /* key1 of the parameter entry */
    string key2;                 /* key2 of the parameter entry */
    string key3;                 /* key3 of the parameter entry */
    string key;
    //    int Counter;                           /* atom counter */
    //    int i;

    while(getline(param_file, line)) {
      /* strip off the comment of the line */
      line = strip(rm_comment(line));
      
      /* if the line is shorter than the total length of 3 keys,
       * then it has no value, we proceed to the next line */
      if (line.size() < (LEN_KEY)) continue;

      /* create 3 keys */
      key1 = strip(line.substr(0,                 LEN_KEY1));
      key2 = strip(line.substr(LEN_KEY1,          LEN_KEY2));
      key3 = strip(line.substr(LEN_KEY1+LEN_KEY2, LEN_KEY3));
      
      /* Then we convert the value string to appropiate type */

      /* integer type */
      if ( key1 == "NATOM" || key1 == "IATOM" ) {
         int* value = new int (atoi(line.substr(LEN_KEY).c_str()));
         save(key1, key2, key3, value, sizeof(int));
      }
      /* float */
      else if (key1 = "PROTON") ||
               key1 = "ELECTRON")
          !strcmp(key1, "PKA      ") ||
               !strcmp(key1, "CHARGE   ") ||
               !strcmp(key1, "RADIUS   ") ||
               !strcmp(key1, "RXN      ") ||
               !strcmp(key1, "VDWAMBER ") ||
               !strcmp(key1, "RADCOVAL ") ||
               !strcmp(key1, "EM       ") ||
               !strcmp(key1, "EXTRA    ") ||
               !strcmp(key1, "SCALING  ")
              ) {    /* if falls into one of them, convert string to float */
         float value;
         value = atof(line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         param_sav(key1, key2, key3, &value, sizeof(float));
      }
      /* bool */
      else if (!strcmp(key1, "CAL_VDW  ") ||
               !strcmp(key1, "RELAX    ")
              ) {
          int value;
          if (strchr(line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 't') ||
              strchr(line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 'T') ) value = 1;
          else value = 0;
         param_sav(key1, key2, key3, &value, sizeof(int));
      }
      /* array of strings */
      if (!strcmp(key1, "CONFLIST ")) {
         STRINGS value;
         strip(sbuff, line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         value.n = strlen(sbuff)/6+1;  /* 5 char for each conformer name, 1 separation space, # of conformers */

         /* allocate string array from the number of conf names */
         value.strings = (char **) malloc(value.n*sizeof(char *));

         /* loop over all conformer names in the value string */
         for (Counter=0; Counter<value.n; Counter++) {
            /* allocate the string for 5 char in each conf entry */
            value.strings[Counter] = (char *) malloc(6*sizeof(char));
            /* copy the conformer name */
            strncpy(value.strings[Counter], sbuff+Counter*6, 5);
            value.strings[Counter][5] = '\0';
         }
         /* save this STRINGS structure */
         param_sav(key1, key2, key3, &value, sizeof(value));
      }
      /* connectivity */
      else if (!strcmp(key1, "CONNECT  ")) {
         CONNECT connect;              /* temporary CONNECT record */
         CONNECTED_ATOM connected;     /* temporary CONNECTED_ATOM record */

         /* get the orbital name from column LEN_KEY1+LEN_KEY2+LEN_KEY3, 10 characters */
         strncpy(sbuff, line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 10); sbuff[10] = '\0';
         strip(connect.orbital, sbuff);

         /* set the line pointer to the first atom */
         if (strlen(line) <= LEN_KEY1+LEN_KEY2+LEN_KEY3 + 10)
            ptr = line+strlen(line);                    /* no atom, step to the end */
         else
            ptr = line+LEN_KEY1+LEN_KEY2+LEN_KEY3 + 10; /* 1st atom, step to the 1st atom */

         /* set atom counter */
         Counter = 0;

         /* start looping over of atoms */
         while (strlen(ptr)) {
            strncpy(sbuff, ptr, 5); sbuff[5] = '\0';  /* cut the first 5 char */
            connected.res_offset = atoi(sbuff);       /* convert to integer */
            if (strstr(sbuff, "LIG")) connected.ligand = 1;/* ligand flag is on */
            else connected.ligand = 0;                     /* ligand flag is off */

            strncpy(connected.name, ptr+5, 4); /* cut the next 4 char as the atom name */
            connected.name[4] = '\0';
            while (strlen(connected.name)<4)   /* If atom name is too short, complete with space */
                strcat(connected.name," ");

            connect.atom[Counter] = connected; /* save this connected atom */
            Counter ++;                        /* Counter steps forward by 1 */
            if (strlen(ptr) > 10) ptr += 10;   /* string pointer steps forward by 1 atom */
            else ptr += strlen(ptr);           /* string pointer steps forward to the end */
         }

         /* store the value */
         connect.n = Counter;
         param_sav(key1, key2, key3, &connect, sizeof(CONNECT));

      }
      else if (!strcmp(key1, "ROTAMER  ") ||
               !strcmp(key1, "HYDROXYL ")) { /* convert to ROTAMER */
         ROTAMER value;
         strncpy(value.atom1, line+LEN_KEY1+LEN_KEY2+LEN_KEY3, 4); value.atom1[4] = '\0';
         strncpy(value.atom2, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+5, 4); value.atom2[4] = '\0';
         strip(value.affected, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+10);
         /* save this STRINGS structure */
         param_sav(key1, key2, key3, &value, sizeof(value));
      }
      else if (!strcmp(key1, "TORSION  ")) {
          TORS value;
          char  sbuffer[10];
          int i_term;
          strncpy(value.atom1, line+LEN_KEY1+LEN_KEY2+LEN_KEY3,    4); value.atom1[4] = '\0';
          strncpy(value.atom2, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+5,  4); value.atom2[4] = '\0';
          strncpy(value.atom3, line+LEN_KEY1+LEN_KEY2+LEN_KEY3+10, 4); value.atom3[4] = '\0';
          strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+15, 4);     sbuffer[4] = '\0';
          if (strchr(sbuffer, 't') ||
              strchr(sbuffer, 'T') ) value.opt_hyd = 1;
          else value.opt_hyd = 0;
          for(i=strlen(line)-1;i>=0;i--) {
              if (line[i]==' ') line[i]='\0';
              else if (line[i]=='\t') line[i]='\0';
              else break;
          }
          value.n_term = ((float)(strlen(line)-(LEN_KEY1+LEN_KEY2+LEN_KEY3+20)))/30.+1;
          if (value.n_term < 1) value.n_term = 1;
          for (i_term=0;i_term<value.n_term; i_term++) {
              strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+20+i_term*30, 9);     sbuffer[9] = '\0'; value.V2[i_term]    = atof(sbuffer);
              strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+30+i_term*30, 9);     sbuffer[9] = '\0'; value.n_fold[i_term]  = atof(sbuffer);
              strncpy(sbuffer,     line+LEN_KEY1+LEN_KEY2+LEN_KEY3+40+i_term*30, 9);     sbuffer[9] = '\0'; value.gamma[i_term] = env.d2r * atof(sbuffer);
          }
         /* save */
         param_sav(key1, key2, key3, &value, sizeof(TORS));
      }
      /* the rest converted to string */
      else {
         strcpy(sbuff, line+LEN_KEY1+LEN_KEY2+LEN_KEY3);
         param_sav(key1, key2, key3, sbuff, sizeof(sbuff));
         /* save with the string terminating character */
      }
   }


   fclose(fp);
   return 0;
}
