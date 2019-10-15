#ifndef __NODE_HPP
#define __NODE_HPP

#include "ca_graph.hpp"
#include <ilcplex/ilocplex.h>
#include <string>

#define RC_EPS 1.0e-6
using namespace std;

class Node {
	
	/// ATTRIBUTES ///
	private:
	string       m_name;     // Nom du fichier chargé.
	int          m_id;       // Identifiant de la node.

	gc::ca_graph m_graph;    // Graphe correspondant au fichier chargé (et modifié).
	int          m_nbUser;	 // Nombre de "vertices" dans le graphe chargé.
	int          m_nbColor;  // Nombre de couleurs à disposition.

	IloEnv       m_env;      // Variable d'environnement de la node.
	IloModel     m_master;   // Problème maître
	IloModel     m_sub;      // Sous-problème

	float        m_score;    // Score obtenu
	int          m_colors[]; // Couleur de chaque node


	/// CONSTRUCTORS ///
	public:
	Node(): m_name("unknown"),
	        m_id(0),
	        m_nbUser(0),
		m_nbColor(0),
	        m_score(0) {}
	
	/// METHODS ///
	//private:

	public:
	void load(const string fileName); /*
	* - Lit le fichier "name" et charge son contenu dans Node.m_graph.
	* - Extrait le radical du nom de fichier et le stock dans Node.m_name.
	*/
	
	void begin(); /*
	* - Crée la variable d'environnement Node.m_env.
	*/	

	void init(const int nbColor); /*
	* - Initialise le problème principale.	
	*/

	void solve(); /*
	* - Résout le problème par la méthode de génération de colonne.
	*/

	void end(); /*
	* - Détruit la variable d'environnement.
	*/

	void retrieveColor(); /*
	* - Récupère le résultat de coloration de notre problème.
	* - Ne récupère que les résultat en nombre entier (defaut:  pas colorier)
	*/	

	void sol() const; /*
	* - Crée le fichier Node.m_name+".sol".
	* - Le remplit avec les résultats trouvés.
	*/

	void dot() const; /*
	* - Crée le fichier Node.m_name+".dot".
	* - Le remplit avec les résultats trouvés.
	*/
};

#endif
