#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "von_neumann.h"
#include "aes.h"
#include "mersenne_twister.h"

#define ARRAY_MAX_SIZE 1000
#define OLDRAND_MAX 2147483647
#define NB_IT 1024


 
struct file_attente
{
	double * arrivee; 
	int taille_arrivee; 
	double *depart;
	int taille_depart;
}; 
typedef struct file_attente file_attente;

struct evolution
{
	double *temps; 
	unsigned int *nombre; 
};
typedef struct evolution evolution; 



static unsigned int next;
static u32 Kex[NB*NR];   //Pour AES
static u32 Px[NB];   //Pour AES


int rdtsc()
{
	// cette fonction suivante cree un warning : c'est normal.
	__asm__ __volatile__("rdtsc");
}

void oldinit_rand(unsigned int seed)
{
	next = seed;
}

unsigned int poidsfort(unsigned int nb)
{
	return (nb >> 27) & 0x0F ;
}

unsigned int poidsfaible(unsigned int nb)
{
	return nb & 0x0F ;
}

unsigned int oldrand()
{
	next = next * 1103515245 + 12345;
	return (unsigned int)(next % OLDRAND_MAX);
}

double Frequency ( int n, unsigned int * tab )
{
	int s = 0, i, j;
	for ( i = 0; i < NB_IT ; i++ )
	{
		for ( j = 0 ; j < n ; j++ )
		{
			if ( ((tab[i] >> j) & 0x01) == 1 )
			{
				s++;
			}
			else
			{
				s--;
			}
		}
	}
	double Sobs = abs(s)/sqrt(NB_IT*n);
	//double Pvaleur = erfc(Sobs/sqrt(2));
	return erfc(Sobs/sqrt(2));
}

double Run(int n, unsigned int * tab){
	int s = 0;
	int i, j;
	int Vn = 0;
	
	for ( i = 0; i < NB_IT ; i++ )
	{
		for ( j = 0 ; j < n ; j++ )
		{
			if ( ((tab[i] >> j) & 0x01) == 1 )
			{
				s++;
			}
		}
	}

	
	double propUn = (double)s/(double)(NB_IT*n);

	if( fabs(propUn - (0.5)) >= (2/sqrt(NB_IT*n)))
	{
		return 0.0;
	}
	
	
	for ( i = 0; i < NB_IT ; i++ )
	{
		for ( j = 0 ; j < n-1 ; j++ )
		{
			if ( ((tab[i] >> j) & 0x01) != ((tab[i] >> (j+1)) & 0x01) )
			{
				Vn++;
			}
		}
		
		if(i < (NB_IT-1)){
			if ( ((tab[i+1] >> j) & 0x01) != ((tab[i]) & 0x01)  )
			{
				Vn++;
			}
		}
	}

	Vn++;
	
	
	double Pvaleur = erfc( abs( Vn - 2*(NB_IT)*n*propUn*(1-propUn)) / (2*sqrt(2*(NB_IT)*n)*propUn*(1-propUn)));
	
	/*
	printf("Pi %f \n",propUn); 
	printf("taux %f\n", (2/sqrt(NB_IT*n)));
	printf("abs pi - 1/2 %f\n", fabs(propUn - (0.5)));
	printf("Vn %d\n",Vn);
	printf("Pvaleur %f\n",Pvaleur);
	*/
	return Pvaleur;
}


double Alea(){
	return (double) AES(Px, Kex) / UINT32_MAX ;
}

double Exponentielle(double lambda)
{
	return -(log(1-Alea())/lambda);
}

file_attente FileMM1 (double lambda, double mu, double D)
{
		//Construction de la file d'arrivee
		double *arrivee = (double*)calloc(ARRAY_MAX_SIZE, sizeof(double)); 
		int taille_arrivee = 0; 
		double ti = 0.0f; 

		for(;;)
		{
			ti += Exponentielle(lambda); 
			if(ti>D)
				break;
			arrivee[taille_arrivee] = ti ; 
			taille_arrivee ++; 
		}
		
		//Construction de la file de départ
		double *depart = (double*)calloc(ARRAY_MAX_SIZE,sizeof(double));
		int taille_depart = 0;
		double tprimei = 0;
		
		
		//On remplit pour le premier client séparement
		tprimei = arrivee[0] + Exponentielle(mu) ; 
		depart[0] = tprimei; 
		taille_depart++; 
		
		
		int i; 
		for(i = 1 ; i<taille_arrivee; i++)
		{
			double tiplusUn = arrivee[i];
			tprimei = depart[i-1]; 
			
			if(tiplusUn > tprimei ) 
			{
				depart[i] = arrivee[i] + Exponentielle(mu); 
				if(depart[i] < D){
					taille_depart++; 
				}
			}else
			{
				depart[i] = depart[i-1] + Exponentielle(mu); 
				if(depart[i] < D){
					taille_depart++; 
				}			
			}
		}

		
		
		file_attente a = {arrivee,taille_arrivee,depart, taille_depart};
		return a; 
}

evolution Calcul_evolution(file_attente a)
{

	double *temps = (double*)calloc(a.taille_arrivee + a.taille_depart,sizeof(double));
	int *nombre = (int*)calloc(a.taille_arrivee + a.taille_depart,sizeof(double));
	
	int indiceArrivee = 0; 
	int indiceDepart = 0; 
	int indiceEvolution = 0; 
	

	//Premier tour de boucle 
	temps[indiceEvolution] = a.arrivee[indiceArrivee]; 
	nombre[indiceEvolution] = 1; 
	indiceEvolution++;
	indiceArrivee++; 

	
	//Autres
	while( (indiceArrivee+indiceDepart) < (a.taille_arrivee + a.taille_depart) )
	{
		//Si on a atteint la fin du tableau d'arrivee, on remplit qu'avec les depart
		//Rq: On ne peut jamais atteindre la fin du tableau départ avant celui d'arrivee
		if(indiceArrivee == a.taille_arrivee){
			temps[indiceEvolution] = a.depart[indiceDepart]; 
			nombre[indiceEvolution] = nombre[indiceEvolution-1] - 1; 
			indiceDepart++; 
		}else{
	
			if(a.arrivee[indiceArrivee] < a.depart[indiceDepart] )
			{
				temps[indiceEvolution] = a.arrivee[indiceArrivee]; 
				nombre[indiceEvolution] = nombre[indiceEvolution-1] + 1; 
				indiceArrivee++; 
			}else
			{
				temps[indiceEvolution] = a.depart[indiceDepart]; 
				nombre[indiceEvolution] = nombre[indiceEvolution-1] - 1; 
				indiceDepart++; 
			}
			
		}
		
		indiceEvolution++;
	}
	
	evolution e = {temps,nombre};
	return e; 
}


double Calcul_nombre_moyen(evolution e, int taille)
{
	double nbMoyen = 0; 
	double denominateur = 0; 
	
	int i; 
	
	for(i =0; i< taille-1 ; i++)
	{
		nbMoyen += ((e.temps[i+1] - e.temps[i])* (double) (e.nombre[i]));
		denominateur += (e.temps[i+1] - e.temps[i]) ; 
	}

	return nbMoyen/denominateur; 
}

double Calcul_temps_moyen(file_attente a)
{
	double tpsMoyen = 0; 
	
	int i ; 
	for ( i = 0 ; i < a.taille_depart ; i++)
	{
		tpsMoyen += a.depart[i] - a.arrivee[i]; 
	}	
	
	return tpsMoyen/(double)a.taille_depart ; 
}

int main()
{
	word16 x=1111; // nombre entre 1000 et 9999 pour Von Neumann
	struct mt19937p mt; // Pour Mersenne-Twister
	int tmp = rand(), seed, i; // Pour Mersenne-Twister
	u32 Kx[NK]; // pour l'AES, le reste en haut en static
	unsigned int Test[NB_IT],tabFort[NB_IT], tabFaible[NB_IT], tabVN[NB_IT], tabMT[NB_IT], tabAES[NB_IT];

	unsigned int output_rand; // sortie du rand du C	
	unsigned int output_rand_fort;	//sortie du rand poids fort
	unsigned int output_rand_faible; //sortie du rand poids faible
	word32 output_AES; // sortie pour l'AES
	word16 output_VN; // sortie pour pour Von Neumann
	word32 output_MT; // sortie pour Mersenne-Twister
	

               
	// initialisation des graines des generateurs

	srand(rdtsc());  // rand du C 
	seed = rand();
	oldinit_rand(seed);
	sgenrand(time(NULL)+(tmp), &mt); // Mersenne-Twister
	// Initialisation de la clé et du plaintext pour l'AES 
	// 45 est un paramètre qui doit changer à chaque initialisation
	init_rand(Kx, Px, NK, NB, 45);
	KeyExpansion(Kex,Kx); // AES : sous-clefs
	
	FILE* fichier = NULL;


	
	//---GENERATION des nombres
	
		fichier = fopen("rand_c.txt", "w");
		if(fichier != NULL)
		{
			// sorties des generateurs	
			for (i=0; i < NB_IT; i++)
			{
				output_rand= oldrand(); // Von Neumann
				fprintf(fichier,"%u\n",output_rand);
			}
		}
		close(fichier);
		
		
		fichier = fopen("VN.txt", "w");
		if(fichier != NULL)
		{
			for (i=0; i < NB_IT; i++)
			{
				output_VN = Von_Neumann(&x); // Von Neumann
				fprintf(fichier,"%u\n",output_VN);
			}
		}
		fclose(fichier);
		
		fichier = fopen("MT.txt", "w");
		if(fichier != NULL)
		{
			for (i=0; i < NB_IT; i++)
			{
				output_MT = genrand(&mt); // Mersenne-Twister
				fprintf(fichier,"%u\n",output_MT);
			}
		}
		fclose(fichier);
		
		fichier = fopen("AES.txt", "w");
		if(fichier != NULL)
		{
			for (i=0; i < NB_IT; i++)
			{
				output_AES = AES(Px, Kex); // AES
				fprintf(fichier,"%u\n",output_AES);
			}
		}
		fclose(fichier);
		
		fichier = fopen("rand_fort.txt", "w");
		if(fichier != NULL)
		{
			for (i=0; i < NB_IT; i++)
			{
				output_rand_fort = poidsfort(oldrand());
				fprintf(fichier,"%u\n",output_rand_fort);
			}
		}
		fclose(fichier);
		
		fichier = fopen("rand_faible.txt", "w");
		if(fichier != NULL)
		{
			for (i=0; i < NB_IT; i++)
			{
				output_rand_faible = poidsfaible(oldrand());
				fprintf(fichier,"%u\n",output_rand_faible);
			}
		}
		fclose(fichier);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
			
		
		
		//-------------------Generation pour les rands
	fichier = fopen("tests.txt","w");
	if(fichier != NULL)
	{
	int ij; 
	for(ij = 0; ij < 20; ij++)
	{
		int i; 
		for (i=0; i < NB_IT; i++)
		{
			output_rand = oldrand(); // rand du C
			tabFort[i] = poidsfort(output_rand);
			tabFaible[i] = poidsfaible(output_rand);
			tabVN[i] = Von_Neumann(&x); // Von Neumann
			tabMT[i] = genrand(&mt); // Mersenne-Twister
			tabAES[i] = AES(Px, Kex); // AES
		}
	
		
		
		//---------------------Test Frequence monobit
		fprintf(fichier, "\n\n-----Test %d----- \n\n",ij);
		fprintf(fichier,"\n\n Frequences \n\n");
		fprintf(fichier,"Frequence Von Neumann : %f\n",Frequency(13,tabVN));
		fprintf(fichier,"Frequence Mersenne Twister : %f\n",Frequency(32,tabMT));
		fprintf(fichier,"Frequence AES : %f\n",Frequency(32,tabAES));
		fprintf(fichier,"Frequence 4 bits de poids fort : %f\n", Frequency(4,tabFort));
		fprintf(fichier,"Frequence 4 bits de poids faible : %f\n", Frequency(4,tabFaible));
	
	
	
	

		//---------------------Test RUN
		fprintf(fichier,"\n\n Run \n\n");
		fprintf(fichier,"Run Von Neumann : %f\n",Run(13,tabVN));
		fprintf(fichier,"Run Mersenne Twister : %f\n",Run(32,tabMT));
		fprintf(fichier,"Run AES : %f\n",Run(32,tabAES));
		fprintf(fichier,"Run 4 bits de poids fort : %f\n", Run(4,tabFort));
		fprintf(fichier,"Run 4 bits de poids faible : %f\n", Run(4,tabFaible));
	}
	}
	fclose(fichier);

	
	// affichage
	/*
	printf("- Generation de nombres aleatoires -\n");
	printf("rand du C : %u \n",output_rand); 
	printf("Von Neumann : %u\n",output_VN);
	printf("Mersenne Twister : %u\n",output_MT);
	printf("AES : %u\n",output_AES);
	printf("4 bits de poids fort : %u\n", output_rand_fort);
	printf("4 bits de poids faible : %u\n", output_rand_faible);
	*/
	


	//------------------------PARTIE3
	
	
	//printf("%f\n",Exponentielle(1));
	//printf("%f\n",Alea());

	
	
	file_attente a = FileMM1(0.2f,0.33f,180.0f);
	file_attente b = FileMM1(0.2f,2*0.33f,180.0f);
	evolution e = Calcul_evolution(a);
	double nbMoyen = Calcul_nombre_moyen(e, a.taille_arrivee+a.taille_depart);
	double tpsMoyen = Calcul_temps_moyen(a); 
	
		
	fichier = fopen("fileMM1.txt","w");
	if(fichier != NULL)
	{
	int j; 
	fprintf(fichier,"\n\n\ Arrivee               Departs  \n"); 
	for(j = 0; j<a.taille_depart; j++)
	{
		fprintf(fichier,"%f ",a.arrivee[j]);
		fprintf(fichier," \t\t%f \n",a.depart[j]); 	
	}
	for(j = a.taille_depart; j<a.taille_arrivee; j++)
	{
		fprintf(fichier,"%f \t\t \n ",a.arrivee[j]);
	}
	
	
	fprintf(fichier,"\n\n\n Evolution \n\n");
	 
	fprintf(fichier,"Temps    \t\t\t Nombre\n"); 
	for(j = 0; j<(a.taille_arrivee+a.taille_depart) ; j++)
	{
		fprintf(fichier,"%f \t\t\t %d\n",e.temps[j],e.nombre[j]); 
	}
		
	
	fprintf(fichier,"\n\n\n\n\nNb Moyen : %f \n\n",nbMoyen); 
	fprintf(fichier,"\nTps Moyen: %f \n\n",tpsMoyen);
	}
	fclose(fichier);
	
	
	
	
	
	free(a.arrivee); 
	free(a.depart); 
	free(e.temps); 
	free(e.nombre); 
	return 0;
}
