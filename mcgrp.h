#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <malloc.h>
#include <limits.h>
#define my_free(p) \
 do { \
	free(p); \
	p = NULL; \
	}while(0);
#define OPT 280
#define MAXINT INT_MAX
#define OP_NUM  3
#define ARC_NODE 3
#define ARC 2
#define NODE 1
#define NONE 0
#define INSERT 0
#define DINSERT 1
#define SWAP 2
#define TWO_OPT 3
#define INVERSE 4
#define TABU_LENGTH 10
#define STOP 2
#define RAND_STEP 2
#define RAND_STAGE 1
#define TENURE_RANGE 5
//#define alpha 0.9
#define beta 5
//#define gamma 0.2
#define TIME_OUT 1200
#define ITER_MAX 20000


typedef struct path
{
	int a;
	int serv;
	struct path *next;
} path;

typedef struct node
{
	int req;
	int serv;
} node;

typedef struct arc
{
	int req;
	int serv;
} arc;

typedef struct task
{
	int head;
	int tail;
	int arc;
	int req;
} task;

typedef struct avlMove
{
	int i;
	int target;
	int delta;
} avlMove;

typedef struct boundery
{
	int left;
	int t;
	int right;
} boundery;

int readRaw(FILE *in, int *n);
void Floyd(int n);
int greedy(int n, path **routes);
int checkEdge(int source, int sink, int n);
int countCost(int *solution);
int insertDelta(int index_i, int index_t);
int singleCost(int a, int b);
int checkFill(int *tmp);
void tabuInsert(int *solution, int *stCost, FILE *out);
void buildSinsDelta(int **moveDelta);
void updatePointDelta(int oldIndex_i, int oldIndex_t, int moveType, int upType);
void sinsSol(int *solution, int oldIndex_i, int oldIndex_t);
void searchTask(int *soluiton, int best_i, int *index_i);
void clearTabuTenure();
void setTabuTenure(int index, int iter);
int randomInt(int n);
void findBound(int a, boundery *m, int moveType);
int newPos(int i, int target, int oldPos, int moveType);
void updateRouteDelta(int left, int right, int upType);
void printSol(int *solution);
void inverse(int i);
int inverseDelta(int a);
void updateInverseDelta(int a, int upType, FILE *out);
int notTabu(int i, int iter);
void regIns(avlMove *bestMove, int iter, int *stCost, int ntCost);
void randIns(avlMove *selMove, int iter, int threshold);
void print(FILE *out, avlMove *bestMove, int ntCost, int *stCost, int moveType);
void update(avlMove *bestMove, int iter, FILE *out, int moveType);
void regInv(int *ntCost, int *stCost, FILE *out);
void regSwap(avlMove *bestMove, int iter, int ntCost, int *stCost, FILE *out);
int swapDelta(int i, int j);
void swap(int *sol, int i, int j);
void subDelta(int i, int upType);
int checkSol(int stCost, FILE *out);
int sinsAble(int i, int target, int delta);
void regDoubIns(avlMove *bestMove, int iter, int *stCost, int ntCost);
void buildDinsDelta(int **dinsDelta);
int dinsertDelta(int i, int j);
void dinsSol(int *solution, int i, int j);
int dinsAble(int i, int target, int delta);
void checkDelta(avlMove *bestMove, int moveType, int new_left_i, int new_right_i, int new_left_t, int new_right_t, FILE *out);
void buildSpDelta(int **spDelta);
void reg2Opt(avlMove *bestMove, int iter, int *stCost, int ntCost);
int twoOptDelta(int i, int j);
void twoOpt(int *solution, boundery *i, boundery *j);
void updateTwoOptDelta(int left, int right, int upType);
void update2OptPointDelta(boundery *a, boundery *b, int moveType, int upType);

static int **P;				/* previous node matrix */
static int **D;             /* shortest distance */
static int **G;             /* original graph */
static int **E;				/* label for edges */
static arc **arcs;
static node *nodes;
static int **sinsDelta;		/* 'sinsDelta[a][b] = c' means if task with index 'a' in the initial solution is inserted after 'b', then the cost of solution will increase 'c' */
							/* tips: current solution is not used as index because the elements' indexs of 'solution' change in every iteration so we'll have to update sinsDelta more deeply */
static int **dinsDelta;		/* dinsDelta[soluiton[i]][solution[j]]: the delta of insert solution[i]&solution[i + 1] after solution[j] */
static int **spDelta;		/* spDelta[solution[i]][solution[j]]: the delta of swap solution[i] & soluiton[j] */
static int *tabuTenure;		/* 'tabuTenure[a][b] = c' means task with index 'a' in the initial solution	could not be inserted at index b in c turns */
static int *load;      		/* load of each vechile */
int Q;
int depot;
int taskNum;
int rNum;
task *tasks;
int *routeLens;
int soLen;
int n;                  /* number of nodes */
static int *solution;
static int *index;
//int TENURE_RANGE = 0;
