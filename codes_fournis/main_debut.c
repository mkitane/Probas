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

static unsigned int next;

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
    struct mt19937p mt; // Pour Mersenne-Twister
	int tmp = rand();
    sgenrand(time(NULL)+(tmp), &mt); // Mersenne-Twister


	unsigned int i =genrand(&mt); // Mersenne-Twister
	
	double f = (float)i ; // sur nb max possible en mersenne
	printf("Le nombre est %d \n", i);
    return 0.0f;
}

int main()
{
	word16 x=1111; // nombre entre 1000 et 9999 pour Von Neumann
	struct mt19937p mt; // Pour Mersenne-Twister
	int tmp = rand(), seed, i; // Pour Mersenne-Twister
	u32 Kx[NK], Kex[NB*NR], Px[NB]; // pour l'AES
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
	
	
	/*FILE* fichier = NULL;
	fichier = fopen("rand_c.txt", "w");
	if(fichier != NULL)
	{*/
		// sorties des generateurs	
		
		
		
		
		
		  //---------------------Test Frequence monobit
		for (i=0; i < NB_IT; i++)
		{
			Test[i] = 619;
			output_rand = oldrand(); // rand du C
			//fprintf(fichier,"%u\n",output_rand);
			tabFort[i] = poidsfort(output_rand);
			tabFaible[i] = poidsfaible(output_rand);
			tabVN[i] = Von_Neumann(&x); // Von Neumann
			tabMT[i] = genrand(&mt); // Mersenne-Twister
			tabAES[i] = AES(Px, Kex); // AES
		}
		
		
		/*
	printf("Von Neumann : %f\n",Frequency(13,tabVN));
	printf("Mersenne Twister : %f\n",Frequency(32,tabMT));
	printf("AES : %f\n",Frequency(32,tabAES));
	printf("4 bits de poids fort : %f\n", Frequency(4,tabFort));
	printf("4 bits de poids faible : %f\n", Frequency(4,tabFaible));
	*/
	
	
		//---------------------Test RUN
	printf("Test : %f\n",Run(10,Test));
	printf("Von Neumann : %f\n",Run(13,tabVN));
	printf("Mersenne Twister : %f\n",Run(32,tabMT));
	printf("AES : %f\n",Run(32,tabAES));
	printf("4 bits de poids fort : %f\n", Run(4,tabFort));
	printf("4 bits de poids faible : %f\n", Run(4,tabFaible));
	
	
	
	
	
	
	
	
	//}
	/*close(fichier);
	
	fichier = fopen("VN.txt", "w");
	if(fichier != NULL)
	{
		for (i=0; i < NBIT; i++)
		{
			output_VN = Von_Neumann(&x); // Von Neumann
			fprintf(fichier,"%u\n",output_VN);
		}
	}
	fclose(fichier);
	
	fichier = fopen("MT.txt", "w");
	if(fichier != NULL)
	{
		for (i=0; i < NBIT; i++)
		{
			output_MT = genrand(&mt); // Mersenne-Twister
			fprintf(fichier,"%u\n",output_MT);
		}
	}
	fclose(fichier);
	
	fichier = fopen("AES.txt", "w");
	if(fichier != NULL)
	{
		for (i=0; i < NBIT; i++)
		{
			output_AES = AES(Px, Kex); // AES
			fprintf(fichier,"%u\n",output_AES);
		}
	}
	fclose(fichier);
	
	fichier = fopen("rand_fort.txt", "w");
	if(fichier != NULL)
	{
		for (i=0; i < NBIT; i++)
		{
			output_rand_fort = poidsfort(oldrand());
			fprintf(fichier,"%u\n",output_rand_fort);
		}
	}
	fclose(fichier);
	
	fichier = fopen("rand_faible.txt", "w");
	if(fichier != NULL)
	{
		for (i=0; i < NBIT; i++)
		{
			output_rand_faible = poidsfaible(oldrand());
			fprintf(fichier,"%u\n",output_rand_faible);
		}
	}*/

	// affichage
	/*printf("- Generation de nombres aleatoires -\n");
	printf("rand du C : %u \n",output_rand); 
	printf("Von Neumann : %u\n",output_VN);
	printf("Mersenne Twister : %u\n",output_MT);
	printf("AES : %u\n",output_AES);
	printf("4 bits de poids fort : %u\n", output_rand_fort);
	printf("4 bits de poids faible : %u\n", output_rand_faible);*/

	
	Alea();
	Alea();
	return 0;
}
