#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "von_neumann.h"
#include "aes.h"
#include "mersenne_twister.h"

#define ARRAY_MAX_SIZE 1000
#define OLDRAND_MAX 2147483647

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

unsigned int oldrand()
{
	next = next * 1103515245 + 12345;
	return (unsigned int)(next % OLDRAND_MAX);
}


int main()
{
	word16 x=1111; // nombre entre 1000 et 9999 pour Von Neumann
	struct mt19937p mt; // Pour Mersenne-Twister
	int tmp = rand(), seed; // Pour Mersenne-Twister
	u32 Kx[NK], Kex[NB*NR], Px[NB]; // pour l'AES

	unsigned int output_rand; // sortie du rand du C	
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


	// sorties des generateurs	
	output_rand = oldrand(); // rand du C
	output_VN = Von_Neumann(&x); // Von Neumann
	output_MT = genrand(&mt); // Mersenne-Twister
	output_AES = AES(Px, Kex); // AES

	// affichage
	printf("- Generation de nombres aleatoires -\n");
	printf("rand du C : %u \n",output_rand); 
	printf("Von Neumann : %u\n",output_VN);
	printf("Mersenne Twister : %u\n",output_MT);
	printf("AES : %u\n",output_AES);


	return 0;
}
