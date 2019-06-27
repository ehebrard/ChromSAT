#include "node.hpp"
#include "ca_graph.hpp"
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

ILOSTLBEGIN

/******************************************************************************/
/* void Node::load(const string fileName)                                     *
 * - Lit le fichier "name" et charge son contenu dans Node.m_graph.           *
 * - Extrait le radical du nom de fichier et le stock dans Node.m_name.       */
/******************************************************************************/
void Node::load(const string fileName) {
	// Déclaration 
	int            V(0);     // Vertices
	string         line;     // Buffer de lecture
	vector<string> tks;      // Vecteur de token
	int            a,    b;  // Index des nodes lues
	ifstream       f;        // Flux de lecture entrant
	stringstream   sepline;  // Flux de conversion

	// Ouvrir fichier
	f.open(fileName);
	
	// Vérifier si fichier ouvert
	if (f.is_open()) {
		cout << "Impossible d'ouvrir le fichier en entré." << endl;
		exit(EXIT_FAILURE);
	}

	// Lire les lignes correspondant aux nodes seules
	while (getline(f, line)) {
		if (line == "#") {
			break;
		} else {
			V++;
		}
	}

	// Instancier le graphe
	m_graph = gc::ca_graph(V);

	// Lire les lignes correspondant aux arcs
	while (getline(f, line)) {
		sepline << line;
		sepline >> a >> b;
		m_graph.add_edge(a,b);
	}

	// Fermer fichier
	f.close();
	
	// Séparer et sauvegarder le nom du fichier
	sepline << fileName;
	tks = vector<string>();
	while (getline(sepline, line, '/')) {
		tks.push_back(line);
	} 
	if (!tks.empty()){
		sepline << tks.back();
		getline(sepline, m_name, '.');
	} else {
		m_name = "unknown";
	}
}

/******************************************************************************/
/* void Node::begin()                                                         *
 * - Crée la variable d'environnement Node.m_env.                             */
/******************************************************************************/
void Node::begin() {
	m_env = IloEnv();
}

/******************************************************************************/
/* void Node::end()                                                           *
 * - Détruit la variable d'environnement.                                     */
/******************************************************************************/
void Node::end() {
	m_env.end();
}
/******************************************************************************/
/* void Node::init(const int nbColor)                                         *
 * - Initialise le problème principale.	                                      */
/******************************************************************************/
void Node::init(const int nbColor) {
	// Sauvegarder nbColor
	IloNumArray line;
	m_nbColor = nbColor;

	// Instancier le modèle maître
	m_master            = IloModel(m_env, "Master Problem");

	IloObjective  cover = IloAdd(m_master, IloMinimize(m_env)); //min( S[lambda j] )
	IloRangeArray range = IloAdd(m_master, IloRangeArray(m_env,m_nbUser,0,1));

	IloNumVarArray choice(m_env);
	

	// Instancier le sous-modèle
	m_sub = IloModel(m_env, "Subproblem");

	IloObjective   reducedCost = IloAdd(m_sub, IloMinimize(m_env, 1));
	IloNumVarArray use(m_env, m_nbUser, 0.0, 1.0, ILOINT);

	for (int u = 0 ; u < m_nbUser ; u++) {
		for (int v = u ; v < m_nbUser ; v++) {
			if (m_graph.matrix[u][v] == 1) {
				// Pas de couleurs juxtaposées
				m_sub.add(use[u] + use[v] <= 1.0);
			}
		}
		// Toujours soi-même ou un voisin de coloré
		line = IloNumArray(m_env, m_nbUser);
		for(IloInt i = 0 ; i < m_nbUser ; i++) {
			line[i] = m_graph.matrix[u][i];
		}
		m_sub.add(use[u] + IloScalProd(use, line) >= 1);
	}

	// Crée les solvers :
	IloCplex masterSolver(m_master);
	IloCplex subSolver(m_sub);

	// Initialisation :
 	IloNumArray newPatt(m_env, m_nbUser);
	IloNumArray price(m_env, m_nbUser);

	for(;;) {
		// Optimiser en l'état
		masterSolver.solve();
		
		// Récuperer les valeurs d'évaluation de prix
		for(IloInt i = 0 ; i < m_nbUser ; i++) {
			price[i] = -masterSolver.getDual(range[i]);
		}
		reducedCost.setLinearCoefs(use, price);
	
		// Recherche de nouveaux motifs :
		subSolver.solve();
	
		// Verification si on améliore les résultats :
		if (subSolver.getValue(reducedCost) > -RC_EPS) break;

		// Obtention du nouveau motif
		cout << "New pattern!" << endl;
		subSolver.getValues(newPatt, use);

		// Ajout du nouveau motif
		m_master.add(IloNumVar(cover(1) + range(newPatt)));
	}

	m_master.add(IloConversion(m_env, choice, ILOINT));
	masterSolver.solve();
}

int main() {
	return 0;
}

/******************************************************************************/
/* void Node::retrieveColor()                                                 *
 * - Récupère le résultat de coloration de notre problème.                    *
 * - Ne récupère que les résultat en nombre entier (defaut:  pas colorier)    */
/******************************************************************************/
void Node::retrieveColor() {
	//TBD
}

/******************************************************************************/
/* void Node::sol() const                                                     *
 * - Crée le fichier Node.m_name+".sol".                                      *
 * - Le remplit avec les résultats trouvés.                                   */
/******************************************************************************/
void Node::sol() const {
	//TBD
}

/******************************************************************************/
/* void dot() const                                                           *
 * - Crée le fichier Node.m_name+".dot".                                      *
 * - Le remplit avec les résultats trouvés.                                   */
/******************************************************************************/
void Node::dot() const {
	//TBD
}
