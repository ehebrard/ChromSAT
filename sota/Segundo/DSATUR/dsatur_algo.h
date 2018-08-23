#ifndef DSATUR_algo_H
#define DSATUR_algo_H

#include <iostream>
#include <algorithm>
#include <set>
#include <string>
#include <cmath>

#include <unistd.h>

#include "graphe.h"
#include "bornes.h"
#include "cliquer/graph.h"
#include "cliquer/cliquer.h"

using namespace std;

class DSATUR_{
	public:
		/************************************************************/
		/********************STRUCTURE DE DSATUR*********************/
		/************************************************************/
		C_Graphe G;
		int LB;
		int UB;
		double time_limit;
		double time_spent;
		int regle_sel;

		clock_t start;

		int quit;
		unsigned long long int nombre_noeuds;

		int *meilleure_coloration;
		int *coloration_courante;
		int *dsat_courant;
		int solution_courante;

		/************************************************************/
		/************************CONSTRUCTEUR************************/
		/************************************************************/
		//Constructeur
		DSATUR_(){
			quit = -1;
			time_limit = 0;
			time_spent = 0;
			LB = -1;
			UB = -1;
			solution_courante = 0;
			nombre_noeuds = 1;
		};
		//Destructeur
		virtual ~DSATUR_();

		/************************************************************/
		/*************************FONCTIONS**************************/
		/************************************************************/
		//Algorithme DSATUR
		void DSATUR_preprocessing();
		int DSATUR_algo(C_Graphe &G_param, double time_param, int regle_input, int LB_input, int UB_input);
		void DSATUR_algo_rec(int profondeur);

		int choisir_sommet(int *size_candidats);
		int choisir_sommet_PASS(int size_candidats);
		int same(int v1, int v2);

		void afficher_solution_courante(int profondeur);

		void choisir_couleurs(int choix_sommet, int *choix_couleurs, int *size_couleurs);
		void mise_a_jour_dsat(int choix_sommet, int choix_couleur, int *liste_changements_F, int *size_changements_F);
		void mise_a_jour_dsat_noupdate(int choix_sommet, int choix_couleur);

		int DSATUR_h(C_Graphe &GC);
};
#endif
