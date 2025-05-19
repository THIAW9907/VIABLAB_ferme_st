/*! \file  testZermelo.h
 *
 *
 *  \author: A. D�silles, LATSRE
 *  \brief  Fichier contant la d�finition d'un mod�le de viabilit� pour les tests
 *
 *  Ce fichier  contient les d�finitions  de tous les �l�ments qui d�crivent un
 *  un probl�me de viabilit� concret et  qui sont consid�r�s comme
 *  param�tres par le code; Ce fichier repr�sente ainsi une interface d'entr�e primitive pour le code.
 *
 *  Les �l�ments principaux d�finis ici  sont :
 *    -la  dynamique f(t,x,u)
 *    -les contraintes:
 *      -# sur l'�tat k(t,x)
 *      -# sur le controles U(x)
 *    - la cible c(t,x)
 *    -les fonctions d�finissant l'objectif  d'un probl�me d'optimisation
 *      -# fonction l(t,x,u)
 *      -# fonction m(t,x,u)
 *    - diff�rentes options d�finissant plus pr�cis�ment la nature du probl�me � r�soudre
 *    - diff�rentes options d�finissance la m�thode num�rique � utiliser
 *
 *    Voir plus loins dans les commentaires  du fichier la d�finition de chaque fonction
 *     et de chaque param�tre
 *
 *
 *
 */

#include <math.h>
#include <algorithm>
//#include <stdlib>

#ifndef MODEL_H_
#define MODEL_H_

/*
 * Some model specific parameters can be defined here
 */
#include <stdio.h>

// Définition de la structure des paramètres
typedef struct {
    double r;     // taux de croissance maximal
    double k;     // capacité maximale de biomasse
    double mf;    // taux de récolte des feuilles
    double ms;    // taux de récolte des graines
    double phi;   // proportion de biomasse en feuilles
    double alpha; // apport en nutriments en organique et en mineral
    double beta;  // absorption par la biomasse
    double rho;   // dégradation naturelle du sol
    double gamma; // taux de dépense (capital et eau)
    double R;     // apport en eau
    double E;     // évaporation
   
} Params;

// Hypothèse : bornes globales (peuvent être utilisées ailleurs)
double B_min, C_min, S_min, W_min;

// ----- Fonctions auxiliaires -----

// Efficacité de l'eau (fonction theta)
double theta(double W) {
    return W / (1.0 + W);
}

// Dérivée de theta
double dtheta(double W) {
    double denom = 1.0 + W;
    return 1.0 / (denom * denom);
}

// Fonction de revenu agricole f(B)
double f_function(double B, double k) {
    return B / (1.0 + B / k);
}

// Dérivée de f(B)
double df_dB(double B, double k) {
    double denom = 1.0 + B / k;
    return 1.0 / (denom * denom);
}

// ----- Fonction principale de dynamique -----

void dynamics(double *B, double *S, double *C, double *W, Params p, double dt) {
    double Bfeuille = p.phi * (*B);
    double Bgraines = (1.0 - p.phi) * (*B);

    double dB = p.r * (*B) * (1.0 - (*B) / p.k) * (*S) * theta(*W) - p.mf * Bfeuille - p.ms * Bgraines;
    double dS = p.alpha - p.beta * (*S) * (*B) - p.rho * (*S);
    double dC = f_function(*B, p.k) - p.gamma * (*C);
    double I = p.gamma * (*W);
    double dW = p.R - p.E - I;

    // Mise à jour des variables (Euler)
    *B += dB * dt;
    *S += dS * dt;
    *C += dC * dt;
    *W += dW * dt;
}

// ----- Fonction de calcul de la matrice jacobienne -----

void jacobian(double *x, double **jacob, Params p) {
    // x[0]=B, x[1]=S, x[2]=C, x[3]=W
    double B = x[0], S = x[1], C = x[2], W = x[3];

    double thetaW = theta(W);
    double dthetaW = dtheta(W);
    double mfphi = p.mf * p.phi;
    double ms1phi = p.ms * (1.0 - p.phi);

    // dB/dt
    jacob[0][0] = p.r * (1.0 - 2.0 * B / p.k) * S * thetaW - (mfphi + ms1phi);
    jacob[0][1] = p.r * B * (1.0 - B / p.k) * thetaW;
    jacob[0][2] = 0.0;
    jacob[0][3] = p.r * B * (1.0 - B / p.k) * S * dthetaW;

    // dS/dt
    jacob[1][0] = -p.beta * S;
    jacob[1][1] = -p.beta * B - p.rho;
    jacob[1][2] = 0.0;
    jacob[1][3] = 0.0;

    // dC/dt
    jacob[2][0] = df_dB(B, p.k);
    jacob[2][1] = 0.0;
    jacob[2][2] = -p.gamma;
    jacob[2][3] = 0.0;

    // dW/dt
    jacob[3][0] = 0.0;
    jacob[3][1] = 0.0;
    jacob[3][2] = 0.0;
    jacob[3][3] = -p.gamma;
}


#include <math.h> // pour fabs()

inline void localDynBounds(double *x, double *res, Params p) {
    // x[0] = B (biomasse), x[1] = S (fertilité du sol)
    double B = x[0];
    double S = x[1];
    double W = x[3]; // attention : si x ne contient pas W, il faut adapter

    double Bfeuilles = p.phi * B;
    double Bgraines = (1.0 - p.phi) * B;

    // Variation de B : croissance - récolte
    double dB = p.r * B * (1.0 - B / p.k) * S * theta(W) - p.mf * Bfeuilles - p.ms * Bgraines;

    // Variation de S : apport - absorption - perte
    double dS = p.alpha - p.beta * S * B - p.rho * S;

    res[0] = fabs(dB); // variation absolue de la biomasse
    res[1] = fabs(dS); // variation absolue de la fertilité
}


#include <float.h>  // pour DBL_MAX

// Définition d'une valeur flottante très grande (approx. +∞)
#define PLUS_INF DBL_MAX

// Déclaration des bornes (à définir dans un autre fichier .c ou dans Mon main)
extern double B_min;  // Biomasse minimale
extern double k;      // Capacité maximale de biomasse
extern double S_min;  // Fertilité minimale
extern double S_max;  // Fertilité maximale
extern double C_min;  // Capital économique minimal
extern double W_min;  // Eau minimale

/*!
 * Cette fonction définit l’ensemble admissible K des états (x[0]=B, x[1]=S, x[2]=C, x[3]=W).
 * Elle retourne 1.0 si les contraintes sont satisfaites, sinon une valeur "infinie".
 *
 * @param x Tableau de 4 doubles représentant l’état du système
 * @return 1.0 si admissible, PLUS_INF sinon
 */
inline double constraintsX(double *x)
{
    double B = x[0];
    double S = x[1];
    double C = x[2];
    double W = x[3];

    double res = (B >= B_min && B <= k &&
                  S >= S_min && S <= S_max &&
                  C >= C_min &&
                  W >= W_min)
                 ? 1.0 : PLUS_INF;

    return res;
}


#endif /* VIAB2D_H_ */