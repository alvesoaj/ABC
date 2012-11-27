/*
 * main.cpp
 *
 *  Created on: Nov 06, 2012
 *      Author: zerokol
 */
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;

/* Parâmetros do algoritmo */
#define POPULATION_SIZE 40
#define FOOD_SOURCES_SIZE POPULATION_SIZE/2
#define LIMIT 100
#define MAX_NUM_CYCLES 3000
#define MAX_INTERATIONS 30
/* Parâmetros do problema */
#define PARAMS_SIZE 2
#define LOWER_BOUND -5
#define UPPER_BOUND 5

double foods_matrix[FOOD_SOURCES_SIZE][PARAMS_SIZE];
double food_sources_array[FOOD_SOURCES_SIZE]; /* f is a vector holding objective calculate_function values associated with food sources */
double fitness_array[FOOD_SOURCES_SIZE];
double trail_count_array[FOOD_SOURCES_SIZE];
double probabilities_array[FOOD_SOURCES_SIZE];
double optimum_solution;
double optimun_params_array[PARAMS_SIZE];
bool init_flag = true;

/* Funções auxiliares */
string number_to_String(double n);
double calculate_function(double solution[PARAMS_SIZE]);
double calculate_fitness(double value);
double get_random_number();
void init_bee(int index);
void get_best_source();
void send_employed_bees();
void calculate_probabilities();
void send_onlooker_bees();
void send_scout_bees();

int main(int argc, char *argv[]) {
	int interation = 0; // Inicializar o contador de interações
	srand(time(NULL));

	while (interation < MAX_INTERATIONS) {
		// Iniciar todas as abelhas
		for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
			init_bee(i);
		}
		if (init_flag == true) {
			// pegar uma solução qualquer como a melhor
			optimum_solution = food_sources_array[0];
			for (int i = 0; i < PARAMS_SIZE; i++) {
				optimun_params_array[i] = foods_matrix[0][i];
			}
			// pegar a real melhor solução
			get_best_source();
			init_flag = false;
		}
		for (int cycle = 0; cycle < MAX_NUM_CYCLES; cycle++) {
			send_employed_bees();
			calculate_probabilities();
			send_onlooker_bees();
			get_best_source();
			send_scout_bees();
		}
		string temp = "Interação(" + number_to_String(interation) + ")";
		for (int j = 0; j < PARAMS_SIZE; j++) {
			temp += "X(" + number_to_String(j + 1) + "): "
					+ number_to_String(optimun_params_array[j]) + " ";
		}
		cout << temp << endl;
		interation++;
	}
	return 0;
}

string number_to_String(double n) {
	stringstream out;
	out << n;
	return out.str();
}

double calculate_function(double solution[PARAMS_SIZE]) {
	double result = 0;
	for (int i = 0; i < PARAMS_SIZE; i++) {
		result += pow(solution[i], 2);
	}
	return result;
}

double calculate_fitness(double value) {
	if (value >= 0) {
		return 1 / (value + 1);
	} else {
		return 1 + fabs(value);
	}
}

double get_random_number() {
	return ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
}

void init_bee(int index) {
	double solution[PARAMS_SIZE];
	for (int j = 0; j < PARAMS_SIZE; j++) {
		double r = get_random_number();
		foods_matrix[index][j] = r * (UPPER_BOUND - LOWER_BOUND) + LOWER_BOUND;
		solution[j] = foods_matrix[index][j];
	}
	food_sources_array[index] = calculate_function(solution);
	fitness_array[index] = calculate_fitness(food_sources_array[index]);
	trail_count_array[index] = 0;
}

void get_best_source() {
	for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
		if (optimum_solution < food_sources_array[i]) {
			optimum_solution = food_sources_array[i];
			for (int j = 0; j < PARAMS_SIZE; j++) {
				optimun_params_array[j] = foods_matrix[i][j];
			}
		}
	}
}

void send_employed_bees() {
	double new_solution[PARAMS_SIZE];

	for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
		double r = get_random_number();
		int param_to_modify = (int) (r * PARAMS_SIZE);

		r = get_random_number();
		int neighbour = (int) (r * FOOD_SOURCES_SIZE);

		/* coso a escolhida seja a mesma */
		while (neighbour == i) {
			r = get_random_number();
			neighbour = (int) (r * FOOD_SOURCES_SIZE);
		}

		for (int j = 0; j < PARAMS_SIZE; j++) {
			new_solution[j] = foods_matrix[i][j];
		}

		/* v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
		r = get_random_number();
		new_solution[param_to_modify] = foods_matrix[i][param_to_modify]
				+ (foods_matrix[i][param_to_modify]
						- foods_matrix[neighbour][param_to_modify]) * (r - 0.5)
						* 2;
		/* se ultrapassar os limites*/
		if (new_solution[param_to_modify] < LOWER_BOUND) {
			new_solution[param_to_modify] = LOWER_BOUND;
		}
		if (new_solution[param_to_modify] > UPPER_BOUND) {
			new_solution[param_to_modify] = UPPER_BOUND;
		}
		double new_solution_value = calculate_function(new_solution);
		double new_solution_fitness = calculate_fitness(new_solution_value);

		/* verificar se a nova solução é melhor que a atual */
		if (new_solution_fitness > fitness_array[i]) {
			trail_count_array[i] = 0;
			for (int j = 0; j < PARAMS_SIZE; j++) {
				foods_matrix[i][j] = new_solution[j];
			}
			food_sources_array[i] = new_solution_value;
			fitness_array[i] = new_solution_fitness;
		} else {
			trail_count_array[i] += 1;
		}
	}
}

void calculate_probabilities() {
	double maxfit = fitness_array[0];
	for (int i = 1; i < FOOD_SOURCES_SIZE; i++) {
		if (fitness_array[i] > maxfit) {
			maxfit = fitness_array[i];
		}
	}
	for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
		probabilities_array[i] = (0.9 * (fitness_array[i] / maxfit)) + 0.1;
	}
}

void send_onlooker_bees() {
	int i, j, t;
	i = 0;
	t = 0;
	double solution[PARAMS_SIZE];

	while (t < FOOD_SOURCES_SIZE) {
		double r = get_random_number();
		/*choose a food source depending on its probability to be chosen*/
		if (r < probabilities_array[i]) {
			t++;

			/*The parameter to be changed is determined randomly*/
			r = get_random_number();
			int param_to_modify = (int) (r * PARAMS_SIZE);

			/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
			r = get_random_number();
			int neighbour = (int) (r * FOOD_SOURCES_SIZE);

			/*Randomly selected solution must be different from the solution i*/
			while (neighbour == i) {
				r = get_random_number();
				neighbour = (int) (r * FOOD_SOURCES_SIZE);
			}
			for (j = 0; j < PARAMS_SIZE; j++)
				solution[j] = foods_matrix[i][j];

			/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
			r = get_random_number();
			solution[param_to_modify] = foods_matrix[i][param_to_modify]
					+ (foods_matrix[i][param_to_modify]
							- foods_matrix[neighbour][param_to_modify])
							* (r - 0.5) * 2;

			/*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
			if (solution[param_to_modify] < LOWER_BOUND)
				solution[param_to_modify] = LOWER_BOUND;
			if (solution[param_to_modify] > UPPER_BOUND)
				solution[param_to_modify] = UPPER_BOUND;
			double new_solution_value = calculate_function(solution);
			double new_solution_fitness = calculate_fitness(new_solution_value);

			/*a greedy selection is applied between the current solution i and its mutant*/
			if (new_solution_fitness > fitness_array[i]) {
				/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trail_count_array counter of solution i*/
				trail_count_array[i] = 0;
				for (j = 0; j < PARAMS_SIZE; j++)
					foods_matrix[i][j] = solution[j];
				food_sources_array[i] = new_solution_value;
				fitness_array[i] = new_solution_fitness;
			} else { /*if the solution i can not be improved, increase its trail_count_array counter*/
				trail_count_array[i] = trail_count_array[i] + 1;
			}
		}
		i++;
		if (i == FOOD_SOURCES_SIZE)
			i = 0;
	}
}

/*determine the food sources whose trail_count_array counter exceeds the "LIMIT" value. In Basic ABC, only one scout is allowed to occur in each cycle*/
void send_scout_bees() {
	int maxtrialindex = 0;
	for (int i = 1; i < FOOD_SOURCES_SIZE; i++) {
		if (trail_count_array[i] > trail_count_array[maxtrialindex])
			maxtrialindex = i;
	}
	if (trail_count_array[maxtrialindex] >= LIMIT) {
		init_bee(maxtrialindex);
	}
}
