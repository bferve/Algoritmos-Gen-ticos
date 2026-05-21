#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
using namespace std;

//individuo
typedef vector<bool> TGenes;
typedef struct{
    TGenes genes; //cadena de bits (genotipo)
    double x;  // fenotipo
    double aptitud; //funcion de evaluacion
    double puntuacion; // puntuacion relativa: fitness/sumfitness
    double punt_acumulada; //puntuacion acumulada para sorteos
    bool elite; //indicativo de pertenecer a la elite
} TIndividuo;

//poblacion
typedef vector<TIndividuo> TPoblacion;


/****************
 * Funciones módulo
 */

TPoblacion poblacion_inicial(int tam_pob, int lcrom, double x_min, double x_max);
TIndividuo genera_individuo(int lcrom, double x_min, double x_max);
void seleccion(TPoblacion& pob, int tam_pob);
double adaptacion(TIndividuo& individuo, int lcrom, double x_min, double x_max);
double decod(TGenes genes, int lcrom, double x_min, double x_max);
int bin_dec(TGenes genes, int lcrom);

void reproduccion(TPoblacion& pob, int tam_pob, int lcrom, double prob_cruce, double x_min, double x_max);
void cruce(TIndividuo padre1, TIndividuo padre2, TIndividuo& hijo1, TIndividuo& hijo2, int lcrom, int punto_cruce, double x_min, double x_max);
void mutacion(TPoblacion& pob, int tam_pob, int lcrom, double prob_mutacion, double x_min, double x_max);
void evaluacion(TPoblacion& pob, int tam_pob, int& pos_mejor, double& sumfitness);
