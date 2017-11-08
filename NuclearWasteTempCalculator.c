/* Computational Physics and Modelling Assignment 4- Solving A Partial Differential Equation **
** A program to solve for the temperature around a rod of nuclear waste **
** 10/12/2016 **
** Jamie Bull */

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h> /*access*/
#include <string.h>
#include <limits.h>

#include <gsl/gsl_linalg.h>



#define PROGRAM_NAME "NuclearWastetempCalculator"    
#define OUT_EXT_TXT ".txt"
#define F_NAME_MAX 128    //Max file name of 127 characters.
#define F_EXT 4    //Number of characters in the file extension.
#define F_NUM 2    //Max number of integers allowed to be added onto a file name.
#define BUFFER_SIZE 256    //Size of buffer used in scanf_buff.
#define LINE_MAX_SIZE 16384    //Max size of a single line of data which can be plotted.

#define GNUPLOT "gnuplot"
#define PLOT_SCRIPT "PlotNuclearWaste.script"

/*Default values*/
#define DELTA_R_DEFAULT 0.01    
#define R_MAX_DEFAULT 1.0    //Note: the max radius and time calculated for will be R_MAX-DELTA_R or T_MAX-DELTA_T.
#define DELTA_T_DEFAULT 3.15576e8     //10 Years.
#define T_MAX_DEFAULT 3.15576e9    //100 Years.
#define A_DEFAULT 0.25    //Radius of rod in (m).
#define KAPPA_DEFAULT 2e5
#define HALF_LIFE_DEFAULT 3.15576e9   //100 Years.
#define TEMP_ROD_DEFAULT 1.0    
#define TEMP_INF_LIMIT 300.0    //Temperature at r->infinity.

#define YEAR 3.15576e7    //Seconds in a year
#define LEGEND_LIMIT 14    //Max number of entries before the legend is disabled to ensure the graph can still be read.


/*Names of output file(s)*/
typedef struct outFiles{
	
	char fName[F_NAME_MAX-F_EXT];
	char outFile[F_NAME_MAX];    
	
}outFiles;


/*values for the calculation*/
typedef struct values{
	
	double deltaR;
	double rMax;
	double deltaT;
	double tMax;
	double a;
	double kappa;
	double halfLife;
	double tempRod;
	
}values;


/*Initialises values to to default values/ takes user input to set them*/
static void set_values(values *v){
	
	v->deltaR = DELTA_R_DEFAULT;
	v->rMax = R_MAX_DEFAULT;
	v->deltaT = DELTA_T_DEFAULT;
	v->tMax = T_MAX_DEFAULT;
	v->a = A_DEFAULT;
	v->kappa = KAPPA_DEFAULT;
	v->halfLife = HALF_LIFE_DEFAULT;
	v->tempRod = TEMP_ROD_DEFAULT;
	
	/*Uncomment to allow user to set distance steps, max distance, time steps and max time. Otherwise defaults will be used. **
	**printf("Enter the size of distance steps, maximum distance, time steps and max time to be used (in meters and seconds).\n"); **
	**scanf("%lg %lg %lg %lg", deltaR, rMax, deltaT, tMax); */
}


/*Builds output file name using outFiles structure*/
static void build_f_name(outFiles *file){
	
	snprintf( file->outFile, sizeof(file->outFile), "%s%s", file->fName, OUT_EXT_TXT );
}


/*Safely takes input for a string of size F_NAME_MAX-F_EXT using a buffer*/
static void scanf_buff(char string[F_NAME_MAX-F_EXT]){
	
	char buffer[BUFFER_SIZE];
	scanf("%255s", buffer);
	snprintf(string, F_NAME_MAX-F_EXT, "%s", buffer);
}


/*Checks if a file allready exists, gives user option to overwrite or create new file*/
static void check_file(outFiles *file){
	
	build_f_name(file);
	
	if( access(file->outFile, F_OK)==-1 ){
		return;
	}
	else{
		int validC=0;
		char choice;
		while( validC==0 ) {
			
		    printf("\n%s all-ready exists would you like to:\n"
		           "1. Overwrite %s\n"
		    	   "2. Create a new file\n", file->outFile, file->outFile);
		    scanf(" %c", &choice);
		    if( choice=='1' ){
				validC=1;
		    	FILE *f=fopen(file->outFile, "w");    /*Clearing file*/
		    	if( f==NULL ){
		    		fprintf(stderr, "Error opening %s\n", file->outFile);
		    		exit(errno);
		    	}
		    	else{
					fclose(f);
		    		return; 
		    	}
		    }
		    else if( choice=='2' ){
				validC=1;
		        printf("Enter a new name for your file. (Max %i characters)\n", F_NAME_MAX-F_EXT);
				scanf_buff(file->fName);
				build_f_name(file);
				if( access(file->outFile, F_OK)!=-1 ){
					validC=0;
				}
	        }
			else{
				printf("Invalid choice, please try again.\n");
			}
	    }
	}
	return;
}


/*Writes temperature to a file*/
static void write_temp_file(FILE *f, int arraySizeDist, int arraySizeTime, double **temp, values *v){
	int i,
	    j;
	
	fprintf(f, "#Data produced by NuclearWastetempCalculator.c#\n"
	           "#First row = distance from origin, the rest = temperature with time separation of %lgs#\n\n", v->deltaT);
	for( i=0; i<arraySizeDist; i++ ){
		fprintf(f, "%lg	", ((double)i)*v->deltaR);
		for( j=0; j<arraySizeTime; j++ ){
			fprintf(f, "%lg	", temp[j][i]);
		}
		fprintf(f, "\n");
	}
}


/*Solves a tridiagonal matrix to find temperature as a function of r using GSL*/
static void calc_temp(gsl_vector *t, int arraySizeDist, double prevTemp[arraySizeDist], double sourceTemp[arraySizeDist], values *v, double s){
	
	double diagMidData[arraySizeDist],
	       diagHighData[arraySizeDist-1],
		   diagLowData[arraySizeDist-1],
		   bData[arraySizeDist];
	int i;
	
	for( i=0; i<arraySizeDist; i++ ){
		diagMidData[i]=1.0+2.0*s;
	}
	diagMidData[0]+=-s+(s/2.0);    //Boundary condition at r=0
	gsl_vector_view DiagMid = gsl_vector_view_array(diagMidData, arraySizeDist);
	
	for( i=0; i<arraySizeDist-1; i++ ){
		diagHighData[i]=-s-(s/(2.0*(i+1.0)));
		diagLowData[i]=-s+(s/(2.0*(i+2.0)));
	}
	gsl_vector_view DiagHigh = gsl_vector_view_array(diagHighData, arraySizeDist-1);
	gsl_vector_view DiagLow = gsl_vector_view_array(diagLowData, arraySizeDist-1);
	
	for( i=0; i<arraySizeDist; i++ ){
		bData[i] = prevTemp[i] + (sourceTemp[i] * v->kappa * v->deltaT);
	}
	bData[arraySizeDist-1]=bData[arraySizeDist-1]-(-s-(s/(2.0*(arraySizeDist-1.0))))*TEMP_INF_LIMIT;    //Boundary condition for r->infinity
	gsl_vector_view b = gsl_vector_view_array(bData, arraySizeDist);
	
	gsl_linalg_solve_tridiag( &DiagMid.vector, &DiagHigh.vector, &DiagLow.vector, &b.vector, t );  
}


/*Loop to perform the calculation at different times*/
static void calc_temp_loop(FILE *f, values *v){
	
	int arraySizeDist = v->rMax / v->deltaR,
	    arraySizeTime = v->tMax / v->deltaT,
	    i,
		j;
	
	/*Allocating array for temperature*/
	double **temp=malloc(arraySizeTime*sizeof(double*));
	if( temp==NULL ){
		fprintf(stderr, "Error allocating memory.\n");
		exit(errno);
	}
	for( i=0; i<arraySizeTime; i++ ){
		temp[i]=malloc(arraySizeDist*sizeof(double));
		if( temp[i]==NULL ){
			fprintf(stderr, "Error allocating memory.\n");
			exit(errno);
		}
	}
	double *sourceTemp=malloc(arraySizeDist*sizeof(double));
	if( sourceTemp==NULL ){
		fprintf(stderr, "Error allocating memory.\n");
		exit(errno);
	}
	
	double prevTemp[arraySizeDist],
	       s = v->kappa * v->deltaT / (v->deltaR * v->deltaR),
		   Time=0;
    
	/*Initialising previous temperature (before Time=0)*/
	for( i=0; i<arraySizeDist; i++ ){
	prevTemp[i]=300;
	}
	
	gsl_vector *t = gsl_vector_alloc(arraySizeDist);
	printf("Performing calculation...\n");
	for( i=0; i<arraySizeTime; i++ ){
		Time = i * v->deltaT;
		/*Calculating sourceTemp*/
	    for( j=0; j<arraySizeDist; j++ ){
		    if( j*v->deltaR <= v->a ){    //Inside the rod
		        sourceTemp[j]=(v->tempRod*exp(-Time/v->halfLife))/(v->a*v->a);
		    }
		    else{    //Outside the rod
			    sourceTemp[i]=0;
		    }
	    }
		calc_temp(t, arraySizeDist, prevTemp, sourceTemp, v, s);
		for( j=0; j<arraySizeDist; j++ ){
	        temp[i][j]=gsl_vector_get(t, j);
			prevTemp[j]=temp[i][j];
		}
	}
	gsl_vector_free(t);
	printf("Calculation complete, writing to file...\n");
	write_temp_file(f, arraySizeDist, arraySizeTime, temp, v);
	for( i=0; i<arraySizeTime; i++ ){
		free(temp[i]);
	}
	free(temp);
	free(sourceTemp);
}


/*Writes the script to be used to plot temperature*/
static int write_plot_temp_script(outFiles *file, values *v){
	
	char plotFileName[F_NAME_MAX],
	     plotScript1[LINE_MAX_SIZE],
		 plotScript2[LINE_MAX_SIZE],
		 plotScript3[LINE_MAX_SIZE];
	int i;
		 
	FILE *script=fopen(PLOT_SCRIPT, "w");
	if( script==NULL ){
		fprintf(stderr, "Error unable to open %s.\n", PLOT_SCRIPT);
		exit(errno);
	}
	
	printf("Enter the filename for your graph (max %i characters, without extension)\n", F_NAME_MAX-F_EXT-1);
	scanf_buff(plotFileName);
	strcat(plotFileName, ".pdf");
	
	sprintf(plotScript1, "#Script for plotting data from NuclearWaste.c#\n"
	                "set autoscale\n"
					"set title \"Temperature Surrounding a Rod of Nuclear Waste\"\n"
					"set xlabel \"r(m)\"\n"
				    "set ylabel \"T(K)\"\n"
					"plot ");
	
    if( v->tMax/v->deltaT > LEGEND_LIMIT ){
        printf("Too many entries, disabling legend.\n");		
	    for( i=0; (v->deltaT*i) < v->tMax; i++ ){
		    sprintf(plotScript2, "\"%s\" using 1:%i with lines notitle, ", file->outFile, i+2);
		    if( strlen(plotScript1) + strlen(plotScript2) < LINE_MAX_SIZE ){
		        strcat(plotScript1, plotScript2);
		    }
		    else{
			    printf("Not enough memory to write plot script, consider using a smaller number of time steps.\n");
				fclose(script);
			    return 0;
		    }
	    }
	}
	else{
		for( i=0; (v->deltaT*i) < v->tMax; i++ ){
		    sprintf(plotScript2, "\"%s\" using 1:%i with lines title \"%lg 	Years\", ", file->outFile, i+2, (i*v->deltaT)/YEAR);
		    if( strlen(plotScript1) + strlen(plotScript2) < LINE_MAX_SIZE ){
		        strcat(plotScript1, plotScript2);
		    }
		    else{
			    printf("Not enough memory to write plot script, consider using a smaller number of time steps.\n");
				fclose(script);
			    return 0;
		    }
	    }
	}

	sprintf(plotScript3, "\nset term pdf\n"
                         "set output \"%s\"\n"
		                 "replot\n"
					     "exit", plotFileName);
	
	if( strlen(plotScript1) + strlen(plotScript3) < LINE_MAX_SIZE ){    /*Need to guarantee that script is not cut off at end and arrays do not overflow*/
		strcat(plotScript1, plotScript3);
	}
	else{
		printf("Not enough memory to write plot script.\n");
		fclose(script);
		return 0;
	}
	
	fprintf(script, "%s", plotScript1);
		
	fclose(script);
	return 1;    //Script sucessfully written.
}


/*Plots temperature as a function of r*/
static void plot_temp(FILE *f){
	
	char plotCommand[PATH_MAX];
	
	snprintf(plotCommand, sizeof(plotCommand), "%s --persist %s", GNUPLOT, PLOT_SCRIPT);
	system(plotCommand);
}



int main(){
	
	outFiles *file = malloc(sizeof(outFiles));
	if( file==NULL ){
		fprintf( stderr, "Error allocating memory for structure.\n" );
		exit(errno);
	}
	
	printf("%s:\nCalculates temperature of region surrounding a rod of nuclear waste.\n", PROGRAM_NAME);
	printf("\nPlease enter the desired name for the output file. (max %i characters, without extension)\n", F_NAME_MAX-F_EXT-1);
	scanf_buff(file->fName);
	check_file(file);
	
	values *v = malloc(sizeof(values));
	if( v==NULL ){
		fprintf(stderr, "Error unable to allocate memory for structure.\n");
		exit(errno);
	}
	set_values(v);
	
	FILE *fWrite=fopen(file->outFile, "w");
	if( fWrite==NULL ){
		fprintf( stderr, "Error opening %s.\n", file->outFile );
	}
	calc_temp_loop(fWrite, v);
	fclose(fWrite);
	
	int scriptSuccess=write_plot_temp_script(file, v);
	
	if( scriptSuccess==1 ){
	    FILE *fRead=fopen(file->outFile, "r");
	        if( fRead==NULL ){
		        fprintf(stderr, "Error unable to open %s.\n", file->outFile);
		        exit(errno);
	        }
	
	    if( scriptSuccess==1 ){
	        plot_temp(fRead);
	    }
	
	    fclose(fRead);
	}
	
	free(v);
	free(file);
		
	return 0;
}