/*
 * ABC
 * main.cpp
 *
 *  Created on: Nov 8, 2012
 *      Author: aj.alves@zerokol.com
 */
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;

// BUG Eclipse
// #define CLOCKS_PER_SEC 1000000

/* Parâmetros do algoritmo */
#define POPULATION_SIZE 10
#define FOOD_SOURCES_SIZE POPULATION_SIZE/2
#define LIMIT (POPULATION_SIZE*PARAMS_SIZE)/2
#define MAX_NUM_CYCLES 1000
#define MAX_ITERATIONS 30
/* Parâmetros do problema */
#define PARAMS_SIZE 2
#define UPPER_BOUND 1
#define LOWER_BOUND 0

double bounds_matrix[PARAMS_SIZE][2] = { { -20, 40 }, { -30, 50 } };
double foods_matrix[FOOD_SOURCES_SIZE][PARAMS_SIZE];
double function_array[FOOD_SOURCES_SIZE];
double fitness_array[FOOD_SOURCES_SIZE];
double trail_count_array[FOOD_SOURCES_SIZE];
double probabilities_array[FOOD_SOURCES_SIZE];
double optimum_solution;
double optimun_params_array[PARAMS_SIZE];

/* Funções auxiliares */
string number_to_String(double n);
double calculate_time(clock_t start, clock_t end);
void init();
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
	clock_t time_start = clock();
	int iteration = 0; // Inicializar o contador de interações
	srand(time(NULL)); // Para um randon mais eficiente
	// comandos iniciais
	init();
	// Iniciar ciclos de busca
	while (iteration < MAX_ITERATIONS) {
		for (int cycle = 0; cycle < MAX_NUM_CYCLES; cycle++) {
			send_employed_bees();
			calculate_probabilities();
			send_onlooker_bees();
			send_scout_bees();
			get_best_source();
		}
		string temp = "Iteração(" + number_to_String(iteration) + ")-> ";
		for (int j = 0; j < PARAMS_SIZE; j++) {
			temp += "X(" + number_to_String(j + 1) + ")=" + number_to_String(
					optimun_params_array[j]) + " ";
		}
		cout << temp << endl;
		iteration++;
	}

	cout << "f(x, y) = " << optimum_solution << endl;

	cout << "Tempo de exec (ABC): " << calculate_time(time_start, clock()) << " ms"
			<< endl;
	return 0;
}

string number_to_String(double n) {
	stringstream out;
	out << n;
	return out.str();
}

double calculate_time(clock_t start, clock_t end) {
	return 1000.0 * ((double) (end - start) / (double) CLOCKS_PER_SEC);
}

void init() {
	// Iniciar todas as abelhas
	for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
		init_bee(i);
	}
	// pegar uma solução qualquer como a melhor
	optimum_solution = function_array[0];
	for (int i = 0; i < PARAMS_SIZE; i++) {
		optimun_params_array[i] = foods_matrix[0][i];
	}
	// pegar a real melhor solução
	get_best_source();
}

double calculate_function(double solution[PARAMS_SIZE]) {
	/*
	 // MIN f(x, y) = x^2 + y^2
	 return pow(solution[0], 2) + pow(solution[1], 2);
	 */
	/*
	 // MIN f(x, y) = x^2–x*y+y^2–3*y
	 return pow(solution[0], 2) - solution[0] * solution[1] + pow(solution[1], 2)
	 - 3 * solution[1];
	 */
	/*
	 // MIN f(x, y) = (x-2)^4 + (x - 2y)^2
	 return pow(solution[0] - 2, 4) + pow(solution[0] - 2 * solution[1], 2);
	 */
	// MIN f(x,y) = 100*(y-x^2)^2+(1 -x)^2
	return 100 * pow(solution[1] - pow(solution[0], 2), 2) + pow(
			1 - solution[0], 2);
}

double calculate_fitness(double value) {
	if (value >= 0) {
		return 1 / (value + 1);
	} else {
		return 1 + abs(value);
	}
}

double get_random_number() {
	return ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
}

void init_bee(int index) {
	double solution[PARAMS_SIZE];
	for (int j = 0; j < PARAMS_SIZE; j++) {
		double r = get_random_number();
		foods_matrix[index][j] = r * (bounds_matrix[j][UPPER_BOUND]
				- bounds_matrix[j][LOWER_BOUND])
				+ bounds_matrix[j][LOWER_BOUND];
		solution[j] = foods_matrix[index][j];
	}
	function_array[index] = calculate_function(solution);
	fitness_array[index] = calculate_fitness(function_array[index]);
	trail_count_array[index] = 0;
}

void get_best_source() {
	for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
		if (optimum_solution > function_array[i]) {
			optimum_solution = function_array[i];
			for (int j = 0; j < PARAMS_SIZE; j++) {
				optimun_params_array[j] = foods_matrix[i][j];
			}
		}
	}
}

void send_employed_bees() {
	double new_solution[PARAMS_SIZE];
	int neighbour = 0;

	for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
		double r = get_random_number();
		int param_to_modify = (int) (r * PARAMS_SIZE);

		/* caso a escolhida seja a mesma */
		do {
			r = get_random_number();
			neighbour = (int) (r * FOOD_SOURCES_SIZE);
		} while (neighbour == i);

		/* copiar solução atual */
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
		if (new_solution[param_to_modify]
				< bounds_matrix[param_to_modify][LOWER_BOUND]) {
			new_solution[param_to_modify]
					= bounds_matrix[param_to_modify][LOWER_BOUND];
		}
		if (new_solution[param_to_modify]
				> bounds_matrix[param_to_modify][UPPER_BOUND]) {
			new_solution[param_to_modify]
					= bounds_matrix[param_to_modify][UPPER_BOUND];
		}

		double new_solution_function = calculate_function(new_solution);
		double new_solution_fitness = calculate_fitness(new_solution_function);

		/* verificar se a nova solução é melhor que a atual */
		if (new_solution_fitness > fitness_array[i]) {
			trail_count_array[i] = 0;
			for (int j = 0; j < PARAMS_SIZE; j++) {
				foods_matrix[i][j] = new_solution[j];
			}
			function_array[i] = new_solution_function;
			fitness_array[i] = new_solution_fitness;
		} else {
			trail_count_array[i] += 1;
		}
	}
}

void calculate_probabilities() {
	// pegar o maior fitness
	double maxfit = fitness_array[0];
	for (int i = 1; i < FOOD_SOURCES_SIZE; i++) {
		if (fitness_array[i] > maxfit) {
			maxfit = fitness_array[i];
		}
	}
	// calcular a probabilidade de cada fonte de comida ser escolhida
	for (int i = 0; i < FOOD_SOURCES_SIZE; i++) {
		probabilities_array[i] = (0.9 * (fitness_array[i] / maxfit)) + 0.1;
	}
}

void send_onlooker_bees() {
	int i = 0;
	int t = 0;
	int neighbour = 0;
	double new_solution[PARAMS_SIZE];

	while (t < FOOD_SOURCES_SIZE) {
		double r = get_random_number();
		/* Escolhendo uma fonte de comida, depende da probabilidade */
		if (r < probabilities_array[i]) {
			t++;

			r = get_random_number();
			int param_to_modify = (int) (r * PARAMS_SIZE);

			/* caso a escolhida seja a mesma */
			do {
				r = get_random_number();
				neighbour = (int) (r * FOOD_SOURCES_SIZE);
			} while (neighbour == i);

			/* copiar solução atual */
			for (int j = 0; j < PARAMS_SIZE; j++) {
				new_solution[j] = foods_matrix[i][j];
			}

			/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
			r = get_random_number();
			new_solution[param_to_modify] = foods_matrix[i][param_to_modify]
					+ (foods_matrix[i][param_to_modify]
							- foods_matrix[neighbour][param_to_modify]) * (r
							- 0.5) * 2;

			/* se ultrapassar os limites*/
			if (new_solution[param_to_modify]
					< bounds_matrix[param_to_modify][LOWER_BOUND]) {
				new_solution[param_to_modify]
						= bounds_matrix[param_to_modify][LOWER_BOUND];
			}
			if (new_solution[param_to_modify]
					> bounds_matrix[param_to_modify][UPPER_BOUND]) {
				new_solution[param_to_modify]
						= bounds_matrix[param_to_modify][UPPER_BOUND];
			}

			double new_solution_function = calculate_function(new_solution);
			double new_solution_fitness = calculate_fitness(
					new_solution_function);

			/* verificar se a nova solução é melhor que a atual */
			if (new_solution_fitness > fitness_array[i]) {
				trail_count_array[i] = 0;
				for (int j = 0; j < PARAMS_SIZE; j++) {
					foods_matrix[i][j] = new_solution[j];
				}
				function_array[i] = new_solution_function;
				fitness_array[i] = new_solution_fitness;
			} else {
				trail_count_array[i] += 1;
			}
		}
		i++;
		if (i == FOOD_SOURCES_SIZE) {
			i = 0;
		}
	}
}

void send_scout_bees() {
	int max_trial_index = 0;
	for (int i = 1; i < FOOD_SOURCES_SIZE; i++) {
		if (trail_count_array[i] > trail_count_array[max_trial_index]) {
			max_trial_index = i;
		}
	}
	if (trail_count_array[max_trial_index] >= LIMIT) {
		init_bee(max_trial_index);
	}
}
