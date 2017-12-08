#include <stdio.h>
#include <stdlib.h>
#include "mcgrp.h"

int main()
{
	int i, t, s;
	FILE *in = NULL, *out = NULL;
	path **routes;
	path *pc;
	int source, sink;

	in = fopen("./in.dat", "r");
	if (in)
		readRaw(in, &n);
	else
		printf("Warning: load graph failed!\n");


	/* calculate distance of every two nodes with Floyd-Warshall Algorithm */
	Floyd(n);

	/* find an initial solution with greedy strategy */
	routes = (path **)malloc(taskNum * sizeof(path *));
	rNum = greedy(n, routes);
	assert(rNum > 0);

	//    TENURE_RANGE = gamma * taskNum;
	soLen = taskNum + rNum - 1;
	tasks = (task *)malloc((taskNum + 1) * sizeof(task));
	solution = (int *)malloc((soLen) * sizeof(int));
	routeLens = (int *)malloc(rNum * sizeof(int));
	sinsDelta = (int **)malloc(soLen * sizeof(int *));
	dinsDelta = (int **)malloc(soLen * sizeof(int *));
	spDelta = (int **)malloc(soLen * sizeof(int *));
	tabuTenure = (int *)malloc(soLen * sizeof(int));
	for (i = 0; i < soLen; i++) {
		sinsDelta[i] = (int *)malloc(soLen * sizeof(int));
		dinsDelta[i] = (int *)malloc(soLen * sizeof(int));
		spDelta[i] = (int *)malloc(soLen * sizeof(int));
	}
	tasks[0].head = depot - 1;
	tasks[0].tail = depot - 1;
	tasks[0].arc = 0;
	tasks[0].req = 0;

	int totalLoad = 0;
	int totalCost = 0;
	int stCost = 0;
	int k = 0;
	if (rNum) {
		t = 1;
		s = 0;
		out = fopen("./route.dat", "w");
		assert(out);
		int *cost = (int *)malloc(rNum * sizeof(int));
		int serv = NONE;
		int j = 0, C = 0;
		for (i = 0; i < rNum; i++) {
			routeLens[i] = 0;
			source = depot - 1;
			sink = depot - 1;
			pc = routes[i];
			cost[i] = 0;
			C = 0;
			/* convert routes into one int array */
			while (pc != NULL) {
				sink = pc->a;
				serv = pc->serv;
				if (serv) {
					switch (serv) {
					case NODE: {
						tasks[t].head = sink;
						tasks[t].tail = sink;
						tasks[t].arc = 0;
						tasks[t].req = nodes[sink].req;
						solution[s++] = t++;
						routeLens[i]++;
						break;
					}
					case ARC: {
						assert(source != sink);
						tasks[t].head = source;
						tasks[t].tail = sink;
						tasks[t].arc = 1;
						tasks[t].req = arcs[source][sink].req;
						solution[s++] = t++;
						routeLens[i]++;
						break;
					}
					case ARC_NODE: {
						assert(source != sink);
						tasks[t].head = source;
						tasks[t].tail = sink;
						tasks[t].arc = 1;
						tasks[t].req = arcs[source][sink].req;

						solution[s++] = t++;
						routeLens[i]++;

						tasks[t].head = sink;
						tasks[t].tail = sink;
						tasks[t].arc = 0;
						tasks[t].req = nodes[sink].req;

						solution[s++] = t++;
						routeLens[i]++;
						break;
					}
					default: break;
					}
					fprintf(out, "(serv)");
				}
				fprintf(out, "%d --> ", pc->a + 1);
				cost[i] += G[source][sink];
				source = sink;
				pc = pc->next;
			}
			if (s < soLen)
				solution[s] = 0;
			s++;
			totalCost += cost[i];
			totalLoad += load[i];
			fprintf(out, "end\n");
			fprintf(out, "Cost of rout[%d] is: %d\n", i, cost[i]);
			fprintf(out, "Load of rout[%d] is: %d\n", i, load[i]);
			fprintf(out, "Corresponding solution is ");

			/* debugging code: calculate cost for each route with array solution */
			for (; j < s; j++) {
				if (j == 0)
					C += D[depot - 1][tasks[solution[j]].head];
				else if (j == soLen)
					C += D[tasks[solution[j - 1]].tail][depot - 1];
				else
					C += D[tasks[solution[j - 1]].tail][tasks[solution[j]].head];
				if (j < soLen) {
					fprintf(out, "%d(%d, %d) ", solution[j], tasks[solution[j]].head + 1, tasks[solution[j]].tail + 1);
					C += G[tasks[solution[j]].head][tasks[solution[j]].tail];
				}
			}
			fprintf(out, "\nC: %d\n\n", C);

			/* debugging code: calculate total cost and later to check if it's the same as calculated on routes solution */
			k += C;
		}
		printf("pre-calculated total cost is %d\n", k);
		fprintf(out, "pre-calculated total cost is %d\n", k);
		assert(t == taskNum + 1);
		assert(s == soLen + 1);
		printf("Total cost is %d\nTotal load is %d\n", totalCost, totalLoad);
		fprintf(out, "Total cost is %d\nTotal load is %d\n", totalCost, totalLoad);
		stCost = countCost(solution);
		printf("Shortened total cost is %d\n", stCost);
		fprintf(out, "Shortened total cost is %d\n", stCost);
		free(cost);
	}
	else
		printf("No route has been found!\n");

	int j;
	/* debugging code: check if every task is unique */
/*	for (i = 0; i < soLen; i++) {
		fprintf(out, "%d ", solution[i]);
		for (j = i; j < soLen; j++) {
			if (i != j && solution[i] == solution[j] && solution[i] != 0)
				printf("solution[%d] and solution[%d] are both %d\n", i, j, solution[i]);
		}
	}
*/
	/* reversing phase */
	//inverse(solution, stCost, out);

	/* tabu insert */
	index = (int *)malloc((soLen) * sizeof(int));
	memcpy(index, solution, soLen * sizeof(int));
	for (i = 0; i < soLen; i++)
		solution[i] = i;   // from now on, what solution contains are indexs of array 'index', which is the initial solution

/*	if (checkSol(stCost, out) == 1) {
		printf("Initial solution cost is right!\n");
		fprintf(out, "Initial solution cost is right!\n");
	}
	else {
		printf("Initial solution cost is wrong!\n");
		fprintf(out, "Initial solution cost is wrong!\n");
	}
*/
	buildSinsDelta(sinsDelta);
	buildDinsDelta(dinsDelta);
	buildSpDelta(spDelta);
	clearTabuTenure();
	srand((unsigned)time(NULL));

/*	fprintf(out, "\n\n");
	for(i = 0; i < soLen; i++) {
		for(j = 0; j < soLen; j++) {
			fprintf(out, "%d ", sinsDelta[i][j]);
		}
		fprintf(out, "\n");
	}
*/
	tabuInsert(solution, &stCost, out);
/*	printf("\nYour solution is ");
	if(checkSol(stCost, out) == 1){
		printf("right\n");
		fprintf(out, "right\n");
	}
	else {
		printf("wrong\n");
		fprintf(out, "wrong\n");
	}
	fprintf(out, "\n\n");
*/	
	for (i = 0; i < soLen; i++) {
		fprintf(out, "(%d, %d), ", tasks[index[solution[i]]].head + 1, tasks[index[solution[i]]].tail + 1);
	}

	/* exchange phase */
	//exchange(solution, stCost, out);


	for (i = 0; i < rNum; i++) {
		path *ptr;
		path *head = routes[i];
		while (head) {
			ptr = head;
			head = head->next;
			my_free(ptr);
		}
	}
	for (i = 0; i < n; i++) {
		my_free(D[i]);
		my_free(P[i]);
		my_free(E[i]);
		my_free(G[i]);
		my_free(arcs[i]);
	}
	for (i = 0; i < soLen; i++) {
		my_free(sinsDelta[i]);
		my_free(dinsDelta[i]);
		my_free(spDelta[i]);
	}
	my_free(sinsDelta);
	my_free(dinsDelta);
	my_free(spDelta);
/*	for(i = 0; i < rNum; i++) {
		path *p = routes[i];
		while(p != NULL){
			path *tp = p->next;
			free(p);
			p = tp;
		}
	}
*/    
	my_free(routes);
	my_free(tabuTenure);
	my_free(D);
	my_free(P);
	my_free(E);
	my_free(G);
	my_free(arcs);
	my_free(nodes);
	my_free(tasks);
	my_free(routeLens);
	my_free(load);
	my_free(solution);
	fclose(in);
	fclose(out);
	return 0;
}

void clearTabuTenure() {
	memset(tabuTenure, 0, soLen * sizeof(int));
	return;
}

void setTabuTenure(int index, int iter) {
	//int res = randomInt(TENURE_RANGE)+ (int)(alpha * moveNum);
	int res = randomInt(TENURE_RANGE);
	tabuTenure[solution[index]] = iter + res;
	//if(index + 1 < soLen)
	//tabuTenure[solution[index + 1]] = iter + res;
}

int randomInt(int n)
{
	//    return (int) ((rand() * n)/RAND_MAX);
	return rand() % n;
}

void buildSinsDelta(int **moveDelta) {
	int i, j; // indexes in solution

	for (i = 0; i < soLen; i++) {
		for (j = 0; j < soLen; j++) {
			moveDelta[solution[i]][solution[j]] = insertDelta(i, j);
		}
	}
	return;
}

void tabuInsert(int *solution, int *stCost, FILE *out) {
	int ntCost = *stCost;
	int iter = 0;
	int stagCount = 0, randCount = 0;
	int randInsert = 0;//, tolMoveAvl;
	int threshold;
	int moveType;
	avlMove *bestMove = NULL;
	bestMove = (avlMove *)malloc(sizeof(avlMove));
	time_t start, end;

	fprintf(out, "\n\nStart inserting...\n");
	start = time(NULL);
	end = start;
	while (difftime(end, start) < TIME_OUT && iter < ITER_MAX) {
		iter++;
		threshold = beta * ntCost;
		moveType = randomInt(OP_NUM + 1);
		if (randCount / RAND_STEP < RAND_STAGE) {
			if (stagCount < STOP) {
				//tolMoveAvl = 0;
				if (moveType == INSERT)
					regIns(bestMove, iter, stCost, ntCost);
				else if (moveType == DINSERT)
					regDoubIns(bestMove, iter, stCost, ntCost);
				else if (moveType == SWAP)
					regSwap(bestMove, iter, ntCost, stCost, out);
				else if (moveType == TWO_OPT)
					reg2Opt(bestMove, iter, stCost, ntCost);
				else
					exit(0);
/*				if(i == soLen && target == soLen && best_delta == MAXINT) {
					printf("1st, No Availabe move!\n");
					fprintf(out, "No available move!\n");
					return;
				}
*/				
				if (bestMove->delta >= 0)
					stagCount++;
				else
					stagCount = 0;
//				assert(best_delta != MAXINT);
			}
			else {
				randCount++;
				randInsert++;
				if (randInsert == RAND_STEP) {
					stagCount = 0;
					randInsert = 0;
				}
				fprintf(out, "Random insert here!\n");
				printf("Random insert here!\n");
				moveType = INSERT;
				if (moveType == INSERT)
					randIns(bestMove, iter, threshold);
				else
					continue;
//				randSwap(bestMove, iter, &tolMoveAvl, threshold);
//	            assert(best_delta != MAXINT);
			}
			ntCost += bestMove->delta;
			if (ntCost < *stCost)
				*stCost = ntCost;

			fprintf(out, "iter = %d\n", iter);
			printf("iter = %d\n", iter);
			fprintf(out, "tabu = %d\n", tabuTenure[solution[bestMove->i]] - iter);
			printf("tabu = %d\n", tabuTenure[solution[bestMove->i]] - iter);
			print(out, bestMove, ntCost, stCost, moveType);

			update(bestMove, iter, out, moveType);
		}
		else {
			fprintf(out, "iter = %d\n", iter);
			printf("iter = %d\n", iter);
			randCount = 0;
			regInv(&ntCost, stCost, out);
		}
		end = time(NULL);

/*		printf("Your solution is ");
		fprintf(out, "Your solution is ");
		if(checkSol(ntCost, out) == 1){
			printf("right\n\n");
			fprintf(out, "right\n\n");
		}
		else {
			printf("wrong\n\n");
			fprintf(out, "wrong\n\n");
			return;
		}
*/	
		printf("\n");
		fprintf(out, "\n");
		if (*stCost <= OPT)
			return;
	}
	return;
}

void reg2Opt(avlMove * bestMove, int iter, int * stCost, int ntCost)
{
	int i, j, find = 1;
	int best_i, best_j, best_delta = MAXINT, cur_delta;

	for (i = 0; i < soLen; i++) {
		for (j = i + 2; j < soLen; j++) {
			cur_delta = twoOptDelta(i, j);
			if (cur_delta < MAXINT) {
				if (cur_delta < best_delta) {
					if ((notTabu(i, iter) && notTabu(j, iter)) || ntCost + cur_delta < *stCost) { // if not tabu or shorten total cost
						best_i = i;
						best_j = j;
						best_delta = cur_delta;
						find = 1;
					}
					else
						continue;
				}
				else if (cur_delta == best_delta && notTabu(i, iter) && notTabu(j, iter)) {
					if (randomInt(++find) == 0) {
						best_i = i;
						best_j = j;
					}
				}
			}
		}
	}
	bestMove->i = best_i;
	bestMove->target = best_j;
	bestMove->delta = best_delta;
	return;
}

int twoOptDelta(int i, int j)
{
	assert(j > i + 1);

/*	In this case, it only changes the order of two routes, and the total cost will keep the same. 
/*	So we just return 0 to ignore it.*/
	if (index[solution[i]] == 0 && index[solution[j]] == 0) 
		return 0;

	int intra = 0;
	boundery *bi = NULL, *bj = NULL;
	bi = (boundery *)malloc(sizeof(boundery));
	bj = (boundery *)malloc(sizeof(boundery));
	findBound(i, bi, TWO_OPT);
	findBound(j, bj, TWO_OPT);

	if (bi->left == bj->left && bi->right == bj->right)
		intra = 1;
	else
		intra = 0;

	int *tmp = (int *)malloc(soLen * sizeof(int));
	memcpy(tmp, solution, soLen * sizeof(int));
//	printSol(tmp);
	twoOpt(tmp, bi, bj);
//	printSol(tmp);
//	printf("\n");
	if (checkFill(tmp) != 0) {
		free(tmp);
		return MAXINT;
	}

	int orig = 0, opted = 0, delta = 0;
	if (intra) {
		int m;
		for (m = i; m < j; m++) {
			orig += singleCost(index[solution[m]], index[solution[m + 1]]);
		}
		if (j < soLen - 1)
			orig += singleCost(index[solution[m]], index[solution[m + 1]]);
		else
			orig += singleCost(index[solution[m]], 0);

		for (m = i; m < j; m++) {
			opted += singleCost(index[tmp[m]], index[tmp[m + 1]]);
		}
		if (j < soLen - 1)
			opted += singleCost(index[tmp[m]], index[tmp[m + 1]]);
		else
			opted += singleCost(index[tmp[m]], 0);
	}
	else {
		/* ------- i | --------- Ri - - - - - - - - Lj -------- j | --------- Rj --------*/
		/*         a | b          c|d                           e | f          g|h       */
		int ab, cd, ef, gh;
		int af, gd, eb, ch;
		int a = i, b = i + 1, c = bi->right, d = c + 1;
		int e = j, f = e + 1, g = bj->right, h = g + 1;

//		if (c == e && e == soLen - 1)
//			return 0;
		int eqA = (c == e);
		int eqB = (e == soLen - 1);
		int eqC = (g == soLen);

		if (eqA && eqB && eqC) // In this case, it changes the soluiton in the same way as 'insert Ri after i' does. So we just ignore it and return MAXINT.
			return MAXINT;
//		else if (index[solution[i]] == 0 && bi->right == i + 1 && j == soLen - 1)
//			return MAXINT;
		else if (index[solution[i + 1]] == 0 && j == soLen - 1)	// In this case, nothing will change.
			return 0;
		else {
			ab = singleCost(index[solution[a]], index[solution[b]]);

			if (f >= soLen) {
				ef = singleCost(index[solution[e]], 0);
				af = singleCost(index[solution[a]], 0);
			}
			else {
				ef = singleCost(index[solution[e]], index[solution[f]]);
				af = singleCost(index[solution[a]], index[solution[f]]);
			}

			eb = singleCost(index[solution[e]], index[solution[b]]);

			orig = ab + ef;
			opted = af + eb;
		}

		/* as i and j are in different routes and j > i + 1, c(Ri) will always less than soLen, but d is uncertain */
/*		if (d < soLen) {
			cd = singleCost(index[solution[c]], index[solution[d]]);
			if (g < soLen)
				gd = singleCost(index[solution[g]], index[solution[d]]);
			else
				gd = singleCost(index[solution[g - 1]], index[solution[d]]);
		}
		else {
			printf("You should have returned earlier!\n");
			exit(0);// cd = singleCost(index[solution[c]], 0);			// when d == soLen, e(j) must be equal to c. This kind of move is considered to be no-change-move.
		}

		if (f < soLen) {
			ef = singleCost(index[solution[e]], index[solution[f]]);
			af = singleCost(index[solution[a]], index[solution[f]]);
		}
		else {
			ef = singleCost(index[solution[e]], 0);
			af = singleCost(index[solution[a]], 0);
		}

		if (h < soLen) {
			gh = singleCost(index[solution[g]], index[solution[h]]);
			ch = singleCost(index[solution[c]], index[solution[h]]);
		}
		else {
			gh = 0;
			ch = singleCost(index[solution[c]], 0);
		}

		if (e != c && e != soLen - 1) {
			orig = ab + cd + ef + gh;
			opted = af + gd + eb + ch;
		}
		else if (e != c && e == soLen - 1) {
			orig = ab + cd + ef;
			int ad = singleCost(index[solution[a]], index[solution[d]]);
			int cf = 0;
			opted = ad + eb + cf;
		}
		else if (e == c && e != soLen - 1) {
			orig = ab + cd + gh;
			int gb = (g < soLen) ? singleCost(index[solution[g]], index[solution[b]]) : singleCost(index[solution[g - 1]], index[solution[b]]);
			opted = af + gb + ch;
		}
*/
	}
	delta = opted - orig;

	free(bi);
	free(bj);
	free(tmp);
	return delta;
}

void twoOpt(int * solution, boundery *i, boundery *j)
{
	int m = 0, intra = 0;
	int *tmp = NULL, *tmp2 = NULL;
	int a = i->t, b = a + 1, c = j->t, d = c + 1;
	int task_a = index[solution[a]], task_c = index[solution[c]];

	if (task_a == 0 && task_c == 0)
		return;
/*	else if (task_a == 0)
		a++;
	else if (task_c == 0)
		c++;
*/
	if (i->left == j->left && i->right == j->right)
		intra = 1;
	else
		intra = 0;

	if (intra) {
	/* --------a  break  b-----c  break  d---------*/
		tmp = (int *)malloc((c - a) * sizeof(int));
		memcpy(tmp, solution + a + 1, (c - a) * sizeof(int));
		for (m = b; m <= c; m++)
			solution[m] = tmp[c - m];
		free(tmp);
	}
	else {
		/* ------- i | ---- i_right - - - - - - - - j_left ---- j | ---- j_right -------*/
		/*         a | b   len1    |             len2           c | d   len3    |       */
		/* inter-2opt: len1-->len2-->len3 ==> len3-->len2-->len1						*/

//		if (i->right == c && c == soLen - 1) // in this situation, solution won't change at all
//			return;

		int len1, len2, len3;
		int equA, equB, equC;
		equA = (i->right == j->t);
		equB = (j->t == soLen - 1);
		equC = (j->right == soLen);

		if (equA && equB && equC) // see function twoOptDelta
			return;

		if (equC) {
			len1 = i->right - i->t - 1;
		}
		else {
			len1 = i->right - i->t;
		}

		if (equA && !equB && equC) {
			len2 = 1;
		}
		else if (equA && !equB && !equC) {
			len2 = 0;
		}
		else if (!equA && equC) {
			len2 = j->t - i->right + 1;
		}
		else if(!equA && !equB && !equC)
			len2 = j->t - i->right;
		else
			exit(0);

		if (!equA && equB && equC) {
			len3 = 0;
		}
		else {
			if (!equA && !equB && equC)
				len3 = j->right - j->t - 1;
			else
				len3 = j->right - j->t;
		}

/*		len1 = i->right - a;
		if (c == i->right && j->right == soLen)
			len1 -= 1;
		len2 = (c == i->right && j->right == soLen) ? 1 : c - i->right; // i->right will always be less than soLen, or c will be nowhere. And when these two are equal, len2 will be 1, which represents Ri/c.
		len3 = (c == soLen - 1) ? 0 : ((j->right == soLen) ? (j->right - c - 1) : (j->right - c)); // when c == soLen - 1, j->right will be soLen, and this length should be revised as 0.
*/
		tmp = (int *)malloc(len1 * sizeof(int));
		tmp2 = (int *)malloc(len2 * sizeof(int));

		memcpy(tmp, solution + b, len1 * sizeof(int));

		if (!equA && !equB && !equC)
			memcpy(tmp2, solution + i->right + 1, len2 * sizeof(int));
		else
			memcpy(tmp2, solution + i->right, len2 * sizeof(int));
//		if (len2 == 1)
//			tmp2[0] = solution[i->right];
		memmove(solution + b, solution + d, len3 * sizeof(int));
		memcpy(solution + a + len3 + 1, tmp2, len2 * sizeof(int));
		memcpy(solution + a + len3 + len2 + 1, tmp, len1 * sizeof(int));
		free(tmp);
		free(tmp2);
	}
	return;
}

void regDoubIns(avlMove *bestMove, int iter, int *stCost, int ntCost) {
	int i, target, cur_delta, find = 1;
	int  best_delta, best_i, best_target;

	best_delta = MAXINT;
	cur_delta = MAXINT;
	for (i = 0; i < soLen; i++) {
		for (target = 0; target < soLen; target++) {
			cur_delta = dinsDelta[solution[i]][solution[target]];
			if (dinsAble(i, target, cur_delta)) {
				if (cur_delta < best_delta) {
					if ((notTabu(i, iter) && notTabu(i + 1, iter)) || ntCost + cur_delta < *stCost) { // if not tabu or shorten total cost
						best_i = i;
						best_target = target;
						best_delta = cur_delta;
						find = 1;
					}
					else
						continue;
				}
				else if (cur_delta == best_delta && notTabu(i, iter) && notTabu(i + 1, iter)) {
					if (randomInt(++find) == 0) {
						best_i = i;
						best_target = target;
					}
				}
			}
		}
	}
	bestMove->i = best_i;
	bestMove->target = best_target;
	bestMove->delta = best_delta;

	return;
}

void buildDinsDelta(int **dinsDelta) {
	int i, j;

	for (i = 0; i < soLen; i++) {
		for (j = 0; j < soLen; j++) {
			dinsDelta[solution[i]][solution[j]] = dinsertDelta(i, j);
		}
	}
	return;
}

int dinsertDelta(int i, int j) {
	if (j == i - 1 || j == i || j == i + 1) {
		return 0;
	}
	if (i >= soLen - 1)
		return MAXINT;
//	assert(i >= 0 && i < soLen - 1);
//	assert(j >= 0 && j <= soLen - 1);
	int *tmpSol = (int *)malloc(soLen * sizeof(int));
	memcpy(tmpSol, solution, soLen * sizeof(int));
	dinsSol(tmpSol, i, j);
	if (checkFill(tmpSol) != 0) {
		my_free(tmpSol);
		return MAXINT;
	}

	my_free(tmpSol);


	int a = i - 1, b = i + 1, c = i + 2, d = j + 1;
	int ai, bc, jd;
	int ac, ji, bd;

	if (a >= 0) {
		ai = singleCost(index[solution[a]], index[solution[i]]);
		if (c < soLen)
			ac = singleCost(index[solution[a]], index[solution[c]]);
		else
			ac = singleCost(index[solution[a]], 0);
	}
	else {
		ai = singleCost(0, index[solution[i]]);
		ac = singleCost(0, index[solution[c]]);
	}

	if (c < soLen)
		bc = singleCost(index[solution[b]], index[solution[c]]);
	else
		bc = singleCost(index[solution[b]], 0);

	if (d < soLen) {
		jd = singleCost(index[solution[j]], index[solution[d]]);
		bd = singleCost(index[solution[b]], index[solution[d]]);
	}
	else {
		jd = singleCost(index[solution[j]], 0);
		bd = singleCost(index[solution[b]], 0);
	}

	ji = singleCost(index[solution[j]], index[solution[i]]);

	int k = (ac + ji + bd) - (ai + bc + jd);
	return k;
}

void dinsSol(int *solution, int i, int j) {
	int tmp1, tmp2;
	tmp1 = solution[i];
	tmp2 = solution[i + 1];
	if (i < j) {
		memmove(solution + i, solution + i + 2, (j - i - 1) * sizeof(int));
		solution[j - 1] = tmp1;
		solution[j] = tmp2;
	}
	else if (i > j) {
		memmove(solution + j + 3, solution + j + 1, (i - j - 1) * sizeof(int));
		solution[j + 1] = tmp1;
		solution[j + 2] = tmp2;
	}
	return;
}

int dinsAble(int i, int target, int delta) {
	if (delta == MAXINT)
		return 0;
	else if (i >= soLen - 1 || target >= soLen - 1) // check source index
		return 0;
	else if (target == i || target == i - 1 || target == i + 1) // exclude meaningless insert point
		return 0;
	else if (target == i + 2) { // if three tasks[0] are adjacent, return false to end the meaningless circle
		if (!index[solution[i]] && !index[solution[i + 1]] && !index[solution[target]])
			return 0;
		else
			return 1;
	}
	else
		return 1;
}


void regSwap(avlMove *bestMove, int iter, int ntCost, int *stCost, FILE *out) {
	int i, j, find = 1;
	int best_i, best_j, best_delta = MAXINT, cur_delta;

	for (i = 0; i < soLen; i++) {
		for (j = i + 1; j < soLen; j++) {
			cur_delta = spDelta[solution[i]][solution[j]];
			if (cur_delta < MAXINT) {
				if (cur_delta < best_delta) {
					if ((notTabu(i, iter) && notTabu(j, iter)) || ntCost + cur_delta < *stCost) { // if not tabu or shorten total cost
						best_i = i;
						best_j = j;
						best_delta = cur_delta;
						find = 1;
					}
					else
						continue;
				}
				else if (cur_delta == best_delta && notTabu(i, iter) && notTabu(j, iter)) {
					if (randomInt(++find) == 0) {
						best_i = i;
						best_j = j;
					}
				}
			}
		}
	}
	bestMove->i = best_i;
	bestMove->target = best_j;
	bestMove->delta = best_delta;
	return;
}

void buildSpDelta(int **spDelta) {
	int i, j;

	for (i = 0; i < soLen; i++) {
		for (j = 0; j < soLen; j++) {
			spDelta[solution[i]][solution[j]] = swapDelta(i, j);
		}
	}
/*	for (i = 0; i < soLen; i++) {
		for (j = 0; j < soLen; j++) {
			assert(spDelta[solution[i]][solution[j]] == spDelta[solution[j]][solution[i]]);
		}
	}
*/
	return;
}


int swapDelta(int i, int j) {
	if (index[solution[i]] == 0 && index[solution[j]] == 0)
		return MAXINT;

	int *tmpSol = (int *)malloc(soLen * sizeof(int));
	memcpy(tmpSol, solution, soLen * sizeof(int));

	swap(tmpSol, i, j);
	if (checkFill(tmpSol) != 0) {
		my_free(tmpSol);
		return MAXINT;
	}

	/* ------ --- ------ --- ----- */
	/*      a  i  b    c  j  d     */
	int a = i - 1, b = i + 1;
	int c = j - 1, d = j + 1;
	int ai, ib, cj, jd;
	int aj, jb, ci, id;
	int ij, ji;
	int orig, swaped, delta;

	if (a >= 0) {
		ai = singleCost(index[solution[a]], index[solution[i]]);
		aj = singleCost(index[solution[a]], index[solution[j]]);
	}
	else {
		ai = singleCost(0, index[solution[i]]);
		aj = singleCost(0, index[solution[j]]);
	}
	if (b < soLen) {
		ib = singleCost(index[solution[i]], index[solution[b]]);
		jb = singleCost(index[solution[j]], index[solution[b]]);
	}
	else {
		ib = singleCost(index[solution[i]], 0);
		jb = singleCost(index[solution[j]], 0);
	}
	if (c >= 0) {
		cj = singleCost(index[solution[c]], index[solution[j]]);
		ci = singleCost(index[solution[c]], index[solution[i]]);
	}
	else {
		cj = singleCost(0, index[solution[j]]);
		ci = singleCost(0, index[solution[i]]);
	}
	if (d < soLen) {
		jd = singleCost(index[solution[j]], index[solution[d]]);
		id = singleCost(index[solution[i]], index[solution[d]]);
	}
	else {
		jd = singleCost(index[solution[j]], 0);
		id = singleCost(index[solution[i]], 0);
	}
	orig = ai + ib + cj + jd;
	swaped = aj + jb + ci + id;
	ij = singleCost(index[solution[i]], index[solution[j]]);
	ji = singleCost(index[solution[j]], index[solution[i]]);
	if (j == i + 1) {
		orig = ai + ij + jd;
		swaped = aj + ji + id;
	}
	else if (i == j + 1) {
		orig = cj + ji + ib;
		swaped = ci + ij + jb;
	}
	delta = swaped - orig;
	//    if(delta < 0)
	//        printf("\n\ncheck swap\n\n");
	free(tmpSol);
	return delta;
}

void swap(int *sol, int i, int j) {
	int tmp;

	tmp = sol[i];
	sol[i] = solution[j];
	sol[j] = tmp;

	return;
}
void regInv(int *ntCost, int *stCost, FILE *out) {
	int i, find = 1, best_i = 0;
	int cur_delta = MAXINT, best_delta = MAXINT;
	task taskTmp;
	int upType = 0;

	for (i = 0; i < soLen; i++) {
		taskTmp = tasks[index[solution[i]]];
		if (checkEdge(taskTmp.head, taskTmp.tail, n)) {
			cur_delta = inverseDelta(i);
			if (cur_delta < best_delta) {
				best_delta = cur_delta;
				best_i = i;
				find = 1;
			}
			else if (cur_delta == best_delta) {
				if (randomInt(find++) == 0)
					best_i = i;
			}
		}
	}
	if (best_delta < MAXINT) {
		*ntCost = *ntCost + best_delta;
		if (*ntCost < *stCost)
			*stCost = *ntCost;

		printf("Inverse tasks[%d](index[%d], solution[%d]), original: (%d, %d)\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head + 1, tasks[index[solution[best_i]]].tail + 1, best_delta, *ntCost, *stCost);
		fprintf(out, "Inverse tasks[%d](index[%d], solution[%d]), original: (%d, %d)\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head + 1, tasks[index[solution[best_i]]].tail + 1, best_delta, *ntCost, *stCost);

		inverse(best_i);
		for(upType = INSERT; upType < OP_NUM; upType++)
			updateInverseDelta(best_i, upType, out);
	}
/*	if (checkFill(solution) != 0) {
		printf("Something's wrong with function updatePointDelta!\n");
	}
	checkDelta(NULL, INVERSE, best_i, 0, 0, 0, out);
*/
	return;
}


void update(avlMove *bestMove, int iter, FILE *out, int moveType) {
	int left_i, left_t, right_i, right_t;
	int new_left_i, new_right_i, new_left_t, new_right_t;
	int best_target = bestMove->target;
	int best_i = bestMove->i;
	int best_t = bestMove->target;
	boundery *a = NULL, *b = NULL;
	a = (boundery *)malloc(sizeof(boundery));
	b = (boundery *)malloc(sizeof(boundery));
//	int new_i = newPos(best_i, best_target, best_i, moveType);
//	int new_t = newPos(best_i, best_target, best_target, moveType);
	if ((moveType == INSERT || moveType == DINSERT) && index[solution[best_target]] == 0)
		best_t = best_target + 1;

	findBound(best_i, a, moveType); // when DINSERT, right boundery is the right boundery of 'best_i + 1'
	findBound(best_t, b, (moveType != TWO_OPT) ? INSERT : TWO_OPT); // even DINSERT, right boundery will always be the right boundery of 'best_t'
	left_i = a->left;
	right_i = a->right;
	left_t = b->left;
	right_t = b->right;

	/* when insert, if target's boundery is 'i' or 'i + 1', update it */
	if (left_t == best_i || (moveType == DINSERT && left_t == best_i + 1))
		left_t = left_i;
	if (right_t == best_i || (moveType == DINSERT && right_t == best_i + 1))
		right_t = right_i;

	/* find new index of bounderies in solution */
	if (moveType == INSERT || moveType == DINSERT) {
		new_left_i = ((left_i == -1) ? 0 : newPos(best_i, best_target, left_i, moveType));
		new_right_i = newPos(best_i, best_target, right_i, moveType);
		new_left_t = ((left_t == -1) ? 0 : newPos(best_i, best_target, left_t, moveType));
		new_right_t = newPos(best_i, best_target, right_t, moveType);
	}
	else if (moveType == SWAP) {
		new_left_i = (left_t == -1) ? 0 : ((left_t == best_i) ? left_i : left_t);
		new_right_i = (right_t == best_i) ? right_i : right_t;
		new_left_t = (left_i == -1) ? 0 : ((left_i == best_t) ? left_t : left_i);
		new_right_t = (right_i == best_t) ? right_t : right_i;
	}
	else if (moveType == TWO_OPT) {
		if (left_i == left_t && right_i == right_t) {
			new_left_i = (left_i == -1) ? 0 : left_i;
			new_right_i = right_i;
			new_left_t = (left_t == -1) ? 0 : left_t;
			new_right_t = right_t;
		}
		else {
			new_left_i = (left_i == -1) ? 0 : left_i;
			new_right_i = best_i + right_t - best_t;
			new_left_t = new_right_i + left_t - right_i;
			new_right_t = right_t;
		}
	}
	else {
		new_left_i = 0;
		new_right_i = 0;
		new_left_t = 0;
		new_right_t = 0;
	}
	setTabuTenure(best_i, iter);
	if (moveType == SWAP)
		setTabuTenure(best_target, iter);
	else if (moveType == DINSERT)
		setTabuTenure(best_i + 1, iter);

	fprintf(out, "tabu set = %d\n", tabuTenure[solution[best_i]] - iter);
	printf("tabu set = %d\n", tabuTenure[solution[best_i]] - iter);
	if (moveType == INSERT)
		sinsSol(solution, best_i, best_target);
	else if (moveType == DINSERT)
		dinsSol(solution, best_i, best_target);
	else if (moveType == SWAP)
		swap(solution, best_i, best_target);
	else if (moveType == TWO_OPT)
		twoOpt(solution, a, b);
/*	if (checkFill(solution) != 0) {
		printf("Something's wrong with function updatePointDelta!\n");
	}

	if((new_left_i >= new_right_i || new_left_t >= new_right_t) && (new_left_i | new_right_i) && (new_left_t | new_right_t))
		printf("Buddy, your calculation of positions was wrong!\n");
*/
	int upType = 0;
	if (moveType != TWO_OPT) {
		/* update sinsDelta[][], dinsDelta[][] and spDelta[][] */
		for (upType = INSERT; upType <= SWAP; upType++) {
			updateRouteDelta(new_left_i, new_right_i, upType);
			updateRouteDelta(new_left_t, new_right_t, upType);
			updatePointDelta(best_i, best_target, moveType, upType);
		}
	}
	else {
		int intra = (left_i == left_t && right_i == right_t) ? 1 : 0;
		if (intra) {
			for (upType = INSERT; upType <= SWAP; upType++) {
				updateRouteDelta(new_left_i, new_right_i, upType);
				updateTwoOptDelta(new_left_i, new_right_i, upType);
			}
		}
		else {
			for (upType = INSERT; upType <= SWAP; upType++) {
				updateRouteDelta(new_left_i, new_right_i, upType);
				updateRouteDelta(new_left_t, new_right_t, upType);
				update2OptPointDelta(a, b, moveType, upType);
			}
		}
	}

	//checkDelta(bestMove, moveType, new_left_i, new_right_i, new_left_t, new_right_t, out);
	
/*	if(checkFill(solution) != 0) {
		printf("Something's wrong with function updatePointDelta!\n");
	}
*/	
	free(a);
	free(b);
	return;
}

void updateTwoOptDelta(int left, int right, int upType)
{
	int i, j;
	int **delta = NULL; // delta matrix
	int(*deltaFunc)(int, int); // function used to calculate deltas

	delta = (upType == INSERT) ? sinsDelta :
		(upType == DINSERT ? dinsDelta :
		(upType == SWAP ? spDelta : NULL));
	deltaFunc = (upType == INSERT) ? &insertDelta :
		(upType == DINSERT ? &dinsertDelta :
		(upType == SWAP ? swapDelta : NULL));

	if (delta == NULL || deltaFunc == NULL) {
		printf("In function updateRouteDelta: wrong delta matrix or deltaFunc!\n");
		exit(0);
	}

	int tmp = (right >= soLen) ? right - 1 : right;

	/* update deltas which represent "operate tasks on route left to right to all other tasks" ([Li, Ri + 1]) */
	if (upType != SWAP) { // swap has alread been updated in updateRouteDelta.
		for (i = left; i <= tmp; i++) {
			for (j = 0; j < soLen; j++)
				delta[solution[i]][solution[j]] = (*deltaFunc)(i, j);
		}
	}
	return;
}

void update2OptPointDelta(boundery * a, boundery * b, int moveType, int upType)
{
	if (index[solution[a->t]] == 0 && index[solution[b->t]] == 0)
		return;

	int len1, len2, len3;
	int equA, equB, equC;
	equA = (a->right == b->t);
	equB = (b->t == soLen - 1);
	equC = (b->right == soLen);

	if (equC) {
		len1 = a->right - a->t - 1;
	}
	else {
		len1 = a->right - a->t;
	}

	if (equA && !equB && equC) {
		len2 = 1;
	}
	else if (equA && !equB && !equC) {
		len2 = 0;
	}
	else if (!equA && equC) {
		len2 = b->t - a->right + 1;
	}
	else if (!equA && !equB && !equC)
		len2 = b->t - a->right;
	else
		exit(0);

	if (!equA && equB && equC) {
		len3 = 0;
	}
	else {
		if (!equA && !equB && equC)
			len3 = b->right - b->t - 1;
		else
			len3 = b->right - b->t;
	}
	//int len1 = a->right - a->t, len2 = b->t - a->right, len3 = b->right - b->t;
	int updateNeeded[14];

	updateNeeded[0] = a->t;									// i
	updateNeeded[1] = updateNeeded[0] + len2 + len3 + 1;	// i + 1
	updateNeeded[2] = b->right;								// i_right
	updateNeeded[3] = updateNeeded[0] + len3 + 1;			// i_right + 1
	updateNeeded[4] = updateNeeded[0] + len2 + len3;		// j
	updateNeeded[5] = updateNeeded[0] + 1;					// j + 1
	updateNeeded[6] = updateNeeded[0] + len3;				// j_right
	updateNeeded[7] = b->right + 1;							// j_right + 1
	updateNeeded[8] = updateNeeded[6] - 1;					// j_right - 1
	updateNeeded[9] = updateNeeded[2] - 1;					// i_right - 1
	updateNeeded[10] = a->left + 1;							// i_left + 1
	updateNeeded[11] = updateNeeded[4] - (b->t - b->left);	// j_left + 1
	updateNeeded[12] = updateNeeded[0] - 1;					// i - 1
	updateNeeded[13] = updateNeeded[4] - 1;					// j - 1
	
	int i;
	for (i = 0; i < 12; i++)
		subDelta(updateNeeded[i], upType);
	if (upType == DINSERT || upType == SWAP) {
		for (; i < 14; i++)
			subDelta(updateNeeded[i], upType);
	}

	return;
}

/* upType is used to decide which matrix should be updated with which function */
void updateRouteDelta(int left, int right, int upType) {
//	if(left >= right && (left | right))
//      printf("left >= right?!!\n");

/*	update those deltas(update after solution has been updated):
	1. 'left' and 'right' insert after all other tasks;
	2. all tasks insert into path between 'left' and 'right'.
*/
	int i, j;
	int **delta = NULL; // delta matrix
	int(*deltaFunc)(int, int); // function used to calculate deltas

	delta = (upType == INSERT) ? sinsDelta : 
			(upType == DINSERT ? dinsDelta : 
			(upType == SWAP ? spDelta : NULL));
	deltaFunc = (upType == INSERT) ? &insertDelta : 
				(upType == DINSERT ? &dinsertDelta : 
				(upType == SWAP ? swapDelta : NULL));

	if (delta == NULL || deltaFunc == NULL) {
		printf("In function updateRouteDelta: wrong delta matrix or deltaFunc!\n");
		exit(0);
	}

	/* update deltas which represent "insert all other tasks into route left to right" */
	int tmp = (right >= soLen) ? right - 1 : right;
	for (i = left; i <= tmp; i++) {
		for (j = 0; j < soLen; j++)
			delta[solution[j]][solution[i]] = (*deltaFunc)(j, i);
	}
	/* update deltas which represent "insert boundery of route left to right after all other tasks */
	for (i = 0; i < soLen; i++) {
		delta[solution[left]][solution[i]] = (*deltaFunc)(left, i);
	}
	if (right < soLen) {
		for (i = 0; i < soLen; i++)
			delta[solution[right]][solution[i]] = (*deltaFunc)(right, i);
	}
	/* when 'double insert', boundery should be extend to 'boundery' and 'boundery - 1' */
	if (upType == DINSERT) {
		if (left - 1 >= 0) {
			for (i = 0; i < soLen; i++)
				delta[solution[left - 1]][solution[i]] = (*deltaFunc)(left - 1, i);
		}
		if (right - 1 >= 0) {
			for (i = 0; i < soLen; i++)
				delta[solution[right - 1]][solution[i]] = (*deltaFunc)(right - 1, i);
		}
	}
	/* when SWAP, symetric updates should also be executed */
	if (upType == SWAP) {
		for (i = left; i <= tmp; i++)
			for (j = 0; j < soLen; j++)
				delta[solution[i]][solution[j]] = delta[solution[j]][solution[i]];
		for (i = 0; i < soLen; i++)
			delta[solution[i]][solution[left]] = delta[solution[left]][solution[i]];
		if (right < soLen) {
			for (i = 0; i < soLen; i++)
				delta[solution[i]][solution[right]] = delta[solution[right]][solution[i]];
		}
	}
	return;
}

/* update deltas of the tasks reltaed to best_i and best_target */
/* moveType is used to decide which deltas in the matrix should be updated */
/* upType is used to decide which matrix should be updated with which function */
void updatePointDelta(int oldIndex_i, int oldIndex_t, int moveType, int upType) {
	int newIndex_1, newIndex_2, newIndex_3;
	int newIndex_4, newIndex_5, newIndex_6;
	int newIndex_7, newIndex_8, newIndex_9;

	if (oldIndex_i < oldIndex_t || moveType == SWAP) {
		newIndex_1 = oldIndex_i - 1; 	// for moveType(INSERT, DINSERT and SWAP)&upType(INSERT, DINSERT)
		newIndex_2 = oldIndex_t;		// for moveType(INSERT, DINSERT and SWAP)&upType(INSERT, DINSERT)
		newIndex_3 = oldIndex_i;		// for moveType(INSERT, DINSERT and SWAP)&upType(INSERT, DINSERT)
		newIndex_4 = oldIndex_t - 1;	// for moveType(INSERT, DINSERT and SWAP)&upType(INSERT, DINSERT)
		newIndex_5 = oldIndex_t + 1;	// for moveType(INSERT, DINSERT and SWAP)&upType(INSERT, DINSERT)
		newIndex_6 = oldIndex_i + 1;	// for moveType(SWAP)&upType(INSERT, DINSERT)
		newIndex_7 = oldIndex_t - 2;	// for moveType(INSERT, DINSERT, SWAP)&upType(DINSERT), moveType(DINSERT)&upType(INSERT)
		newIndex_8 = oldIndex_i - 2;	// for moveType(INSERT, DINSERT, SWAP)&upType(DINSERT)
		newIndex_9 = oldIndex_t - 3;	// for moveType(DINSERT)&upType(DINSERT)
	}
	else {
		newIndex_1 = oldIndex_i;		// for moveType(INSERT)&upType(INSERT, DINSERT), moveType(DINSERT)&upType(DINSERT)
		newIndex_2 = oldIndex_t + 1;	// for moveType(INSERT, DINSERT)&upType(INSERT, DINSERT)
		newIndex_3 = oldIndex_i + 1;	// for moveType(INSERT, DINSERT)&upType(INSERT, DINSERT)
		newIndex_4 = oldIndex_t;		// for moveType(INSERT, DINSERT)&upType(INSERT, DINSERT)
		newIndex_5 = oldIndex_t + 2;	// for moveType(INSERT, DINSERT)&upType(INSERT, DINSERT)
		newIndex_6 = oldIndex_i + 2;	// for moveType(DINSERT)&upType(INSERT, DINSERT)
		newIndex_7 = oldIndex_t + 3;	// for mvoeType(DINSERT)&upType(INSERT, DINSERT)
		newIndex_8 = oldIndex_t - 1;	// for moveType(INSERT, DINSERT)&upType(DINSERT)
		newIndex_9 = oldIndex_i - 1;	// for moveType(INSERT)&upType(DINSERT)
	}

	subDelta(newIndex_2, upType);
	subDelta(newIndex_3, upType);
	subDelta(newIndex_4, upType);
	subDelta(newIndex_5, upType);
	if (oldIndex_i < oldIndex_t || moveType == SWAP) {
		subDelta(newIndex_1, upType);
		if (moveType == SWAP)
			subDelta(newIndex_6, upType);
		if (upType == DINSERT || (moveType == DINSERT && (upType == INSERT || upType == SWAP)))
			subDelta(newIndex_7, upType);
		if (upType == DINSERT)
			subDelta(newIndex_8, upType);
		if (moveType == DINSERT && upType == DINSERT)
			subDelta(newIndex_9, upType);
	}
	else {
		if (moveType == INSERT || (moveType == DINSERT && upType == DINSERT))
			subDelta(newIndex_1, upType);
		if (moveType == DINSERT) {
			subDelta(newIndex_6, upType);
			subDelta(newIndex_7, upType);
		}
		if (upType == DINSERT)
			subDelta(newIndex_8, upType);
		if (moveType == INSERT && upType == DINSERT)
			subDelta(newIndex_9, upType);
	}

	return;
}

/* update the deltas of insert 'i' after all other tasks */
void subDelta(int i, int upType) {
	int j;
	int **delta = NULL; // delta matrix
	int(*deltaFunc)(int, int); // function used to calculate deltas

	delta = (upType == INSERT) ? sinsDelta :
		(upType == DINSERT ? dinsDelta :
		(upType == SWAP ? spDelta : NULL));
	deltaFunc = (upType == INSERT) ? &insertDelta :
		(upType == DINSERT ? &dinsertDelta :
		(upType == SWAP ? swapDelta : NULL));

	if (delta == NULL || deltaFunc == NULL) {
		printf("In function subDelta: wrong delta matrix or deltaFunc!\n");
		exit(0);
	}

	if (i >= 0 && i < soLen) {
		for (j = 0; j < soLen; j++) {
			delta[solution[i]][solution[j]] = (*deltaFunc)(i, j);
		}
		if (upType == SWAP) {
			for (j = 0; j < soLen; j++)
				delta[solution[j]][solution[i]] = delta[solution[i]][solution[j]];
		}
	}
	return;
}

void checkDelta(avlMove *bestMove, int moveType, int new_left_i, int new_right_i, int new_left_t, int new_right_t, FILE *out) {
	int i, j;
	int **testDelta = (int **)malloc(soLen * sizeof(int *));
	for (i = 0; i < soLen; i++)
		testDelta[i] = (int *)malloc(soLen * sizeof(int));

	buildSinsDelta(testDelta);
	for (i = 0; i < soLen; i++) {
		for (j = 0; j < soLen; j++) {
			if (moveType != INVERSE) {
				if (testDelta[solution[i]][solution[j]] != sinsDelta[solution[i]][solution[j]]) {
					printf("moveType is %d,\nincorrect delta: sinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove: i = %d, target = %d\nnew-left-i is %d, new-right-i is %d.\nnew-left-t is %d, new-right-t is %d\n\n",
						moveType, solution[i], solution[j],
						testDelta[solution[i]][solution[j]], sinsDelta[solution[i]][solution[j]], i, j,
						bestMove->i, bestMove->target, new_left_i, new_right_i, new_left_t, new_right_t);
					fprintf(out, "moveType is %d,\nincorrect delta: sinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove: i = %d, target = %d\nnew-left-i is %d, new-right-i is %d.\nnew-left-t is %d, new-right-t is %d\n\n",
						moveType, solution[i], solution[j],
						testDelta[solution[i]][solution[j]], sinsDelta[solution[i]][solution[j]], i, j,
						bestMove->i, bestMove->target, new_left_i, new_right_i, new_left_t, new_right_t);
					printSol(solution);
					printf("\n");
					fprintf(out, "\n");
				}
			}
			else {
				if (testDelta[solution[i]][solution[j]] != sinsDelta[solution[i]][solution[j]]) {
					printf("moveType is %d,\nincorrect delta: sinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove(inserse): i = %d\n\n",
						moveType, solution[i], solution[j], testDelta[solution[i]][solution[j]], sinsDelta[solution[i]][solution[j]], i, j, new_left_i);
					fprintf(out, "moveType is %d,\nincorrect delta: sinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove(inserse): i = %d\n\n",
						moveType, solution[i], solution[j], testDelta[solution[i]][solution[j]], sinsDelta[solution[i]][solution[j]], i, j, new_left_i);
					printSol(solution);
					printf("\n");
					fprintf(out, "\n");
				}
			}
		}
	}
	buildDinsDelta(testDelta);
	for (i = 0; i < soLen; i++) {
		for (j = 0; j < soLen; j++) {
			if (moveType != INVERSE) {
				if (testDelta[solution[i]][solution[j]] != dinsDelta[solution[i]][solution[j]]) {
					printf("moveType is %d,\nincorrect delta: dinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove: i = %d, target = %d\nnew-left-i is %d, new-right-i is %d.\nnew-left-t is %d, new-right-t is %d\n\n",
						moveType, solution[i], solution[j],
						testDelta[solution[i]][solution[j]], dinsDelta[solution[i]][solution[j]], i, j,
						bestMove->i, bestMove->target, new_left_i, new_right_i, new_left_t, new_right_t);
					fprintf(out, "moveType is %d,\nincorrect delta: dinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove: i = %d, target = %d\nnew-left-i is %d, new-right-i is %d.\nnew-left-t is %d, new-right-t is %d\n\n",
						moveType, solution[i], solution[j],
						testDelta[solution[i]][solution[j]], dinsDelta[solution[i]][solution[j]], i, j,
						bestMove->i, bestMove->target, new_left_i, new_right_i, new_left_t, new_right_t);
					printSol(solution);
					printf("\n");
					fprintf(out, "\n");
				}
			}
			else {
				if (testDelta[solution[i]][solution[j]] != dinsDelta[solution[i]][solution[j]]) {
					printf("moveType is %d,\nincorrect delta: dinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove(inserse): i = %d\n\n",
						moveType, solution[i], solution[j], testDelta[solution[i]][solution[j]], dinsDelta[solution[i]][solution[j]], i, j, new_left_i);
					fprintf(out, "moveType is %d,\nincorrect delta: sinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove(inserse): i = %d\n\n",
						moveType, solution[i], solution[j], testDelta[solution[i]][solution[j]], dinsDelta[solution[i]][solution[j]], i, j, new_left_i);
					printSol(solution);
					printf("\n");
					fprintf(out, "\n");
				}
			}
		}
	}
	buildSpDelta(testDelta);
	for (i = 0; i < soLen; i++) {
		for (j = 0; j < soLen; j++) {
			if (moveType != INVERSE) {
				assert(testDelta[solution[i]][solution[j]] == testDelta[solution[j]][solution[i]]);
				assert(spDelta[solution[i]][solution[j]] == spDelta[solution[j]][solution[i]]);
				if (testDelta[solution[i]][solution[j]] != spDelta[solution[i]][solution[j]]) {
					printf("moveType is %d,\nincorrect delta: spDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove: i = %d, target = %d\nnew-left-i is %d, new-right-i is %d.\nnew-left-t is %d, new-right-t is %d\n\n",
						moveType, solution[i], solution[j],
						testDelta[solution[i]][solution[j]], spDelta[solution[i]][solution[j]], i, j,
						bestMove->i, bestMove->target, new_left_i, new_right_i, new_left_t, new_right_t);
					fprintf(out, "moveType is %d,\nincorrect delta: spDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove: i = %d, target = %d\nnew-left-i is %d, new-right-i is %d.\nnew-left-t is %d, new-right-t is %d\n\n",
						moveType, solution[i], solution[j],
						testDelta[solution[i]][solution[j]], spDelta[solution[i]][solution[j]], i, j,
						bestMove->i, bestMove->target, new_left_i, new_right_i, new_left_t, new_right_t);
					printSol(solution);
					printf("\n");
					fprintf(out, "\n");
				}
			}
			else {
				if (testDelta[solution[i]][solution[j]] != spDelta[solution[i]][solution[j]]) {
					printf("moveType is %d,\nincorrect delta: spDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove(inserse): i = %d\n\n",
						moveType, solution[i], solution[j], testDelta[solution[i]][solution[j]], spDelta[solution[i]][solution[j]], i, j, new_left_i);
					fprintf(out, "moveType is %d,\nincorrect delta: sinsDelta[%d][%d](T: %d, F: %d), i = %d, j = %d\nbestMove(inserse): i = %d\n\n",
						moveType, solution[i], solution[j], spDelta[solution[i]][solution[j]], spDelta[solution[i]][solution[j]], i, j, new_left_i);
					printSol(solution);
					printf("\n");
					fprintf(out, "\n");
				}
			}
		}
	}

	for (i = 0; i < soLen; i++)
		free(testDelta[i]);
	free(testDelta);

	return;
}

void print(FILE *out, avlMove *bestMove, int ntCost, int *stCost, int printType) {
	int best_i = bestMove->i;
	int best_target = bestMove->target;
	int best_delta = bestMove->delta;

	if (printType == INSERT) {
		fprintf(out, "Insert tasks[%d](index[%d], solution[%d])(%d, %d) after tasks[%d](index[%d], solution[%d])(%d, %d),\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head, tasks[index[solution[best_i]]].tail, index[solution[best_target]], solution[best_target],
			best_target, tasks[index[solution[best_target]]].head, tasks[index[solution[best_target]]].tail, best_delta, ntCost, *stCost);
		printf("Insert tasks[%d](index[%d], solution[%d])(%d, %d) after tasks[%d](index[%d], solution[%d])(%d, %d),\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head, tasks[index[solution[best_i]]].tail, index[solution[best_target]], solution[best_target],
			best_target, tasks[index[solution[best_target]]].head, tasks[index[solution[best_target]]].tail, best_delta, ntCost, *stCost);
	}
	else if (printType == DINSERT) {
		assert(best_i < soLen - 1);
		printf("Double Insert tasks[%d](index[%d], solution[%d])(%d, %d) & tasks[%d](index[%d], solution[%d])(%d, %d) after tasks[%d](index[%d], solution[%d])(%d, %d),\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head, tasks[index[solution[best_i]]].tail, index[solution[best_i + 1]], solution[best_i + 1], best_i + 1, 
			tasks[index[solution[best_i + 1]]].head, tasks[index[solution[best_i + 1]]].tail, index[solution[best_target]], solution[best_target],
			best_target, tasks[index[solution[best_target]]].head, tasks[index[solution[best_target]]].tail, best_delta, ntCost, *stCost);
		fprintf(out, "Double Insert tasks[%d](index[%d], solution[%d])(%d, %d) & tasks[%d](index[%d], solution[%d])(%d, %d) after tasks[%d](index[%d], solution[%d])(%d, %d),\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head, tasks[index[solution[best_i]]].tail, index[solution[best_i + 1]], solution[best_i + 1], best_i + 1,
			tasks[index[solution[best_i + 1]]].head, tasks[index[solution[best_i + 1]]].tail, index[solution[best_target]], solution[best_target],
			best_target, tasks[index[solution[best_target]]].head, tasks[index[solution[best_target]]].tail, best_delta, ntCost, *stCost);

	}
	else if (printType == SWAP) {
		printf("Swap tasks[%d](index[%d], solution[%d])(%d, %d) with tasks[%d](index[%d], solution[%d])(%d, %d)\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head + 1, tasks[index[solution[best_i]]].tail + 1,
			index[solution[best_target]], solution[best_target], best_target, tasks[index[solution[best_target]]].head + 1, tasks[index[solution[best_target]]].tail + 1,
			best_delta, ntCost, *stCost);
		fprintf(out, "Swap tasks[%d](index[%d], solution[%d])(%d, %d) with tasks[%d](index[%d], solution[%d])(%d, %d)\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head + 1, tasks[index[solution[best_i]]].tail + 1,
			index[solution[best_target]], solution[best_target], best_target, tasks[index[solution[best_target]]].head + 1, tasks[index[solution[best_target]]].tail + 1,
			best_delta, ntCost, *stCost);
	}
	else if (printType == TWO_OPT) {
		printf("2-opt: tasks[%d](index[%d], solution[%d])(%d, %d) and tasks[%d](index[%d], solution[%d])(%d, %d)\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head + 1, tasks[index[solution[best_i]]].tail + 1,
			index[solution[best_target]], solution[best_target], best_target, tasks[index[solution[best_target]]].head + 1, tasks[index[solution[best_target]]].tail + 1,
			best_delta, ntCost, *stCost);
		fprintf(out, "2-opt: tasks[%d](index[%d], solution[%d])(%d, %d) and tasks[%d](index[%d], solution[%d])(%d, %d)\ndelta = %d, ntCost = %d, stCost = %d\n",
			index[solution[best_i]], solution[best_i], best_i, tasks[index[solution[best_i]]].head + 1, tasks[index[solution[best_i]]].tail + 1,
			index[solution[best_target]], solution[best_target], best_target, tasks[index[solution[best_target]]].head + 1, tasks[index[solution[best_target]]].tail + 1,
			best_delta, ntCost, *stCost);
	}
	return;
}

void randIns(avlMove *selMove, int iter, int threshold) {
	int i, target, cur_delta = MAXINT;
	int reservior = 0;
	int tolMoveAvl = 0;
	for (i = 0; i < soLen; i++) {
		if (notTabu(i, iter)) {
			for (target = 0; target < soLen; target++) {
				cur_delta = sinsDelta[solution[i]][solution[target]];
				if (i != target && i != target + 1 && cur_delta <= threshold) {
					reservior = randomInt(++tolMoveAvl);
					if (reservior == 0) {
						selMove->i = i;
						selMove->target = target;
						selMove->delta = cur_delta;
					}
				}
			}
		}
	}
}

void regIns(avlMove *bestMove, int iter, int *stCost, int ntCost) {
	int i, target, cur_delta, find = 1;
	int  best_delta, best_i, best_target;

	best_delta = MAXINT;
	cur_delta = MAXINT;
	for (i = 0; i < soLen; i++) {
		for (target = 0; target < soLen; target++) {
			cur_delta = sinsDelta[solution[i]][solution[target]];
			if (sinsAble(i, target, cur_delta)) {
				if (cur_delta < best_delta) {
					if (notTabu(i, iter) || ntCost + cur_delta < *stCost) { // if not tabu or shorten total cost
						best_i = i;
						best_target = target;
						best_delta = cur_delta;
						find = 1;
					}
					else
						continue;
				}
				else if (cur_delta == best_delta && notTabu(i, iter)) {
					if (randomInt(++find) == 0) {
						best_i = i;
						best_target = target;
					}
				}
			}
		}
	}
	bestMove->i = best_i;
	bestMove->target = best_target;
	bestMove->delta = best_delta;

	return;
}

int sinsAble(int i, int target, int delta) {
	if (delta == MAXINT)
		return 0;
	else if (target == i || target == i - 1)
		return 0;
	else if (target == i + 1 && !index[solution[i]] && !index[solution[target]])
		return 0;
	else
		return 1;
}

int notTabu(int i, int iter) {
	if (tabuTenure[solution[i]] <= iter)
		return 1;
	else
		return 0;
}

void inverse(int i) {
	int head, tail;

	head = tasks[index[solution[i]]].head;
	tail = tasks[index[solution[i]]].tail;
	tasks[index[solution[i]]].head = tail;
	tasks[index[solution[i]]].tail = head;

	return;
}

int inverseDelta(int i) {
	int ab, cd, ac, bd;
	int a, b, c, d;

	a = (i == 0) ? tasks[0].tail : tasks[index[solution[i - 1]]].tail;
	b = tasks[index[solution[i]]].head;
	c = tasks[index[solution[i]]].tail;
	d = (i == soLen - 1) ? tasks[0].head : tasks[index[solution[i + 1]]].head;

	ab = D[a][b];
	cd = D[c][d];
	ac = D[a][c];
	bd = D[b][d];

	int delta = (ac + bd) - (ab + cd);

	return delta;
}


void updateInverseDelta(int a, int upType, FILE *out) {
	int i;
	int **delta = NULL; // delta matrix
	int(*deltaFunc)(int, int); // function used to calculate deltas

	delta = (upType == INSERT) ? sinsDelta :
		(upType == DINSERT ? dinsDelta :
		(upType == SWAP ? spDelta : NULL));
	deltaFunc = (upType == INSERT) ? &insertDelta :
		(upType == DINSERT ? &dinsertDelta :
		(upType == SWAP ? swapDelta : NULL));

	if (delta == NULL || deltaFunc == NULL) {
		printf("In function updateInverseDelta: incorrect upType!\n");
		fprintf(out, "In function updateInverseDelta: incorrect upType!\n");
		exit(0);
	}

	for (i = 0; i < soLen; i++) {
		delta[solution[a]][solution[i]] = deltaFunc(a, i);
		delta[solution[i]][solution[a]] = (upType == SWAP) ? delta[solution[a]][solution[i]] : deltaFunc(i, a);
		if (a - 1 >= 0) {
			delta[solution[i]][solution[a - 1]] = deltaFunc(i, a - 1);
			delta[solution[a - 1]][solution[i]] = (upType == SWAP) ? delta[solution[i]][solution[a - 1]] : deltaFunc(a - 1, i);			
			if (a - 2 >= 0 && upType == DINSERT) {
				delta[solution[a - 2]][solution[i]] = deltaFunc(a - 2, i);
				if (upType == SWAP)
					delta[solution[i]][solution[a - 2]] = delta[solution[a - 2]][solution[i]];
			}
		}
		if (a + 1 <= soLen - 1) {
			delta[solution[a + 1]][solution[i]] = deltaFunc(a + 1, i);
			if(upType == SWAP)
				delta[solution[i]][solution[a + 1]] = delta[solution[a + 1]][solution[i]];
		}
	}

	return;
}

void printSol(int *solution) {
	int i;
	for (i = 0; i < soLen; i++)
		printf("%d ", solution[i]);
	printf("\n");

	for (i = 0; i < soLen; i++)
		printf("%d ", index[solution[i]]);
	printf("\n");

	return;
}

void findBound(int a, boundery *m, int moveType) {
	int i = 0;
	int left = -1, right = soLen;

	m->t = a;

	i = ((moveType == TWO_OPT && index[solution[a]] == 0) ? a : a - 1);
	for (; i >= 0; --i) {
		if (index[solution[i]] == 0) {
			left = i;
			break;
		}
	}

	for (i = ((moveType == DINSERT) ? a + 2 : a + 1); i < soLen; ++i) {
		if (index[solution[i]] == 0) {
			right = i;
			break;
		}
	}
	m->left = left;
	m->right = right;
	return;
}

int newPos(int i, int target, int oldPos, int moveType) {
//	if(i == target || i == target + 1)
//		printf("i == target?!!\n");

	if (i < target) {
		if (oldPos < i || oldPos > target)
			return oldPos;
		else if (oldPos == i) {
			if (moveType == INSERT)
				return target;
			else if (moveType == DINSERT)
				return target - 1;
			else {
				printf("In function newPos: moveType error!\n");
				exit(0);
			}
		}
		else if (oldPos == i + 1 && moveType == DINSERT)
			return target;
		else {
			if (moveType == INSERT)
				return oldPos - 1;
			else if (moveType == DINSERT)
				return oldPos - 2;
		}
	}
	else {
		if (oldPos <= target || (oldPos > i && moveType == INSERT) || (oldPos > i + 1 && moveType == DINSERT))
			return oldPos;
		else if (oldPos == i)
			return target + 1;
		else if (oldPos == i + 1 && moveType == DINSERT)
			return target + 2;
		else {
			if (moveType == INSERT)
				return oldPos + 1;
			else if (moveType == DINSERT)
				return oldPos + 2;
			else {
				printf("In function newPos: moveType error!\n");
				exit(0);
			}
		}
	}
	return -1;
}


/* update solution after single insert operation */
void sinsSol(int *solution, int index_i, int index_t) {
	int tmp;
	if (index_i < index_t) {
		tmp = solution[index_i];
		memcpy(solution + index_i, solution + index_i + 1, (index_t - index_i) * sizeof(int));
		solution[index_t] = tmp;
	}
	else if (index_i > index_t) {
		tmp = solution[index_i];
		memcpy(solution + index_t + 2, solution + index_t + 1, (index_i - index_t - 1) * sizeof(int));
		solution[index_t + 1] = tmp;
	}
	return;
}

/* search index of best_i in solution[] */
void searchTask(int *solution, int best_i, int *index_i) {
	int i;
	*index_i = 0;
	for (i = 0; !*index_i && i < soLen; i++) {
		if (solution[i] == best_i)
			*index_i = i;
	}
}

void Floyd(int n) {
	int i, j, k, temp;
/*	static path *pn;
	path *ph;
	FILE *out;

	out = fopen("out.dat", "w");
	if(out == NULL) {
	printf("File error !\n");
	return;
	}
*/
	/* computing procedure */
	printf("Computing shortest distances...\n");
	for (k = 0; k < n; k++)
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				temp = D[i][k] + D[k][j];
				if (D[i][j] > temp && D[i][k] != MAXINT && D[k][j] != MAXINT)
				{
					D[i][j] = temp;
					P[i][j] = P[k][j];
				}
			}

/*	print shortest distance of every too nodes into out.txt */
/*	fprintf(out, "%d\n", n);
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++) {
			fprintf(out, "%d ", D[i][j]);
			if(j == n - 1)
			fprintf(out, "\n");
		}
*/
/*	print shortest path of every two nodes into out.txt */
/*	for(i = 0; i < n; i++)
		for(j = 0; j < n; j) {
			temp = P[i][j];
			if(temp != -1) {
//				printf("%d ", temp);
				pn = (path *)malloc(sizeof(path));
				pn->a = temp;
				pn->next = NULL;
				ph = pn;
				while(temp != i) {
					temp = P[i][temp];
					pn = (path *)malloc(sizeof(path));
					pn->a = temp;
					pn->next = ph;
					ph = pn;
				}
				while(ph != NULL) {
					fprintf(out, "%d ", ph->a);
					ph = ph->next;
				}
				fprintf(out, "%d ", j);
			}
		else if(i == j)
		fprintf(out, "%d", j);
		else
		fprintf(out, "NULL");
		if(!(i == n-1 && j == n-1))
		fprintf(out, "\n");
	}
	fclose(out);
*/
	printf("Shortest distances computing completed !\n");
}


int greedy(int n, path **routes) {
	int i, j, k, cp = depot - 1, np = depot - 1, dp = depot - 1;
	int nextReq, maxReq;
	float maxPro, currentPro;
	static path *nd;
	path *pe;
	int remain = 0;
	int minDist;
	int servType;

	load = (int *)malloc(n * sizeof(int));
	memset(load, 0, n * sizeof(int));
	for (i = 0; i < n; i++) {
		remain += nodes[i].req;
		for (j = 0; j < n; j++) {
			remain += arcs[i][j].req;
		}
	}

	for (k = 0; remain > 0; k++) {
		cp = depot - 1;
		np = depot - 1;
		load[k] = 0;
		maxReq = 0;
		currentPro = 0;
		routes[k] = (path *)malloc(sizeof(path));
		routes[k]->a = depot - 1;
		routes[k]->serv = NONE;
		routes[k]->next = NULL;
		pe = routes[k];

		do {
			/* fnd next node p with largest profit(request/cost) */
			maxPro = 0;
			minDist = MAXINT;
			nextReq = 0;
			maxReq = 0;
			servType = NONE;
			dp = depot - 1;
			for (i = 0; i < n; i++) {
				if (i != cp && G[cp][i] != MAXINT) {             /* exists unserviced request in adjacent nodes */
					if (!nodes[i].serv && !arcs[cp][i].serv)
						nextReq = nodes[i].req + arcs[cp][i].req;
					else if (!nodes[i].serv && arcs[cp][i].serv)
						nextReq = nodes[i].req;
					else if (nodes[i].serv && !arcs[cp][i].serv)
						nextReq = arcs[cp][i].req;
					else
						continue;
					if (nextReq <= Q - load[k]) {                /* make sure route feasible */
						currentPro = nextReq*1.0 / G[cp][i];
						if (currentPro > maxPro) {                   /* current choice has more profit */
							maxPro = currentPro;
							np = i;
						}
					}
				}
			}

			/* if next node hasn't been generated */
			if (np == cp) {
				dp = depot - 1;
				int empty = Q - load[k];
				if (empty > load[k]) {             /* if load is less than half of Q, find the node closest to the node/arc with largest loadabel request */
					maxReq = 0;
					for (i = 0; i < n; i++) {
						if (!nodes[i].serv && nodes[i].req > maxReq && nodes[i].req <= empty) {
							maxReq = nodes[i].req;
							dp = i;                     /* destination point */
						}
						for (j = 0; j < n; j++) {
							if (!arcs[i][j].serv && arcs[i][j].req > maxReq && arcs[i][j].req <= empty) {
								maxReq = arcs[i][j].req;
								dp = i;
							}
						}
					}
				}
				int tmpDist = MAXINT;
				for (i = 0; i < n; i++) {
					if (G[cp][i] != MAXINT && cp != i) {
						tmpDist = G[cp][i] + D[i][dp];
						if (tmpDist < minDist) {
							np = i;
							minDist = tmpDist;
						}
					}
				}
			}

			/* make sure a new node has been added */
			if (np == cp)
				printf("Warning: No move has been taken!\n");

			nd = (path *)malloc(sizeof(path));
			pe->next = nd;
			nd->a = np;
			nd->serv = NONE;
			nd->next = NULL;
			pe = nd;

			/* update load of vechile k and the total remain request */
			if (!nodes[np].serv && !arcs[cp][np].serv) {
				nextReq = nodes[np].req + arcs[cp][np].req;
				servType = ARC_NODE;
			}
			else if (!nodes[np].serv && arcs[cp][np].serv) {
				nextReq = nodes[np].req;
				servType = NODE;
			}
			else if (nodes[np].serv && !arcs[cp][np].serv) {
				nextReq = arcs[cp][np].req;
				servType = ARC;
			}
			else {
				nextReq = 0;
				servType = NONE;
			}
			if (nextReq > 0 && nextReq <= (Q - load[k])) {
				load[k] += nextReq;
				remain -= nextReq;
				nd->serv = servType;
				if (servType == NODE || servType == ARC_NODE)
					nodes[np].serv = servType;
				if (servType == ARC || servType == ARC_NODE) {
					arcs[cp][np].serv = servType;
					if (checkEdge(cp, np, n)) {       /* Required edges only need to be served for one time */
						arcs[np][cp].serv = servType;
						remain -= arcs[np][cp].req;
					}
				}
			}

			cp = np;            /* update current node */
		} while (cp != depot - 1);
		if (remain < 0)
			printf("Warning: Remain is less than zero!\n");
		assert(remain >= 0);
		printf("remain is %d\n", remain);
	}

	return k;
}

int checkEdge(int source, int sink, int n) {
	assert(source >= 0 && source < n);
	assert(sink >= 0 && sink < n);

	if (E[source][sink])
		return 1;
	else
		return 0;
}

int readRaw(FILE *in, int *n) {
	int nodeNum, edgeNum, arcNum, reN, reE, reA;
	int i, j, source, sink, curN;
	int tmp, cost, req;
	char name[50];

	fscanf(in, "Name:	%s\n", name);
	fscanf(in, "Optimal value:   %d\n", &tmp);
	fscanf(in, "#Vehicles:	%d\n", &tmp);
	fscanf(in, "Capacity:	%d\n", &Q);
	fscanf(in, "Depot Node:	%d\n", &depot);
	fscanf(in, "#Nodes:	%d\n", n);
	fscanf(in, "#Edges:	%d\n", &edgeNum);
	fscanf(in, "#Arcs:	%d\n", &arcNum);
	fscanf(in, "#Required N:	%d\n", &reN);
	fscanf(in, "#Required E:	%d\n", &reE);
	fscanf(in, "#Required A:	%d\n", &reA);
	taskNum = reE + reN + reA;

	nodeNum = *n;
	E = (int **)malloc(nodeNum * sizeof(int *));
	D = (int **)malloc(nodeNum * sizeof(int *));
	P = (int **)malloc(nodeNum * sizeof(int *));
	G = (int **)malloc(nodeNum * sizeof(int *));
	arcs = (arc **)malloc(nodeNum * sizeof(arc *));
	nodes = (node *)malloc(nodeNum * sizeof(node));
	for (i = 0; i < nodeNum; i++) {
		G[i] = (int *)malloc(nodeNum * sizeof(int));
		D[i] = (int *)malloc(nodeNum * sizeof(int));
		P[i] = (int *)malloc(nodeNum * sizeof(int));
		E[i] = (int *)malloc(nodeNum * sizeof(int));
		arcs[i] = (arc *)malloc(nodeNum * sizeof(arc));
		nodes[i].req = 0;
		nodes[i].serv = NODE;
		for (j = 0; j < nodeNum; j++) {
			if (i == j) {
				G[i][j] = 0;
				D[i][j] = 0;
			}
			else {
				G[i][j] = MAXINT;
				D[i][j] = MAXINT;
			}
			P[i][j] = -1;
			E[i][j] = 0;
			arcs[i][j].req = 0;
			arcs[i][j].serv = ARC;
		}
	}

	fscanf(in, "\n%[^\n]%*c", name);
	for (i = 0; i < reN; i++) {
		fscanf(in, "N%d	%d	%d\n", &curN, &req, &tmp);
		nodes[curN - 1].req = req;
		nodes[curN - 1].serv = NONE;
	}

	fscanf(in, "\n%[^\n]%*c", name);
	for (i = 0; i < reE; i++) {
		fscanf(in, "E%d	%d	%d	%d	%d	%d\n", &tmp, &source, &sink, &cost, &req, &tmp);
		G[source - 1][sink - 1] = cost;
		G[sink - 1][source - 1] = cost;
		D[source - 1][sink - 1] = cost;
		D[sink - 1][source - 1] = cost;
		P[source - 1][sink - 1] = source - 1;
		P[sink - 1][source - 1] = sink - 1;
		E[source - 1][sink - 1] = 1;
		E[sink - 1][source - 1] = 1;
		arcs[source - 1][sink - 1].req = req;
		arcs[sink - 1][source - 1].req = req;
		arcs[source - 1][sink - 1].serv = NONE;
		arcs[sink - 1][source - 1].serv = NONE;
	}

	fscanf(in, "\n%[^\n]%*c", name);
	int NrE = edgeNum - reE;
	for (i = 0; i < NrE; i++) {
		fscanf(in, "NrE%d	%d	%d	%d\n", &tmp, &source, &sink, &cost);
		G[source - 1][sink - 1] = cost;
		G[sink - 1][source - 1] = cost;
		D[source - 1][sink - 1] = cost;
		D[sink - 1][source - 1] = cost;
		P[source - 1][sink - 1] = source - 1;
		P[sink - 1][source - 1] = sink - 1;
		E[source - 1][sink - 1] = 1;
		E[sink - 1][source - 1] = 1;
	}

	fscanf(in, "\n%[^\n]%*c", name);
	for (i = 0; i < reA; i++) {
		fscanf(in, "A%d	%d	%d	%d	%d	%d\n", &tmp, &source, &sink, &cost, &req, &tmp);
		G[source - 1][sink - 1] = cost;
		D[source - 1][sink - 1] = cost;
		P[source - 1][sink - 1] = source - 1;
		arcs[source - 1][sink - 1].req = req;
		arcs[source - 1][sink - 1].serv = NONE;
	}

	fscanf(in, "\n%[^\n]%*c", name);
	int NrA = arcNum - reA;
	for (i = 0; i < NrA; i++) {
		fscanf(in, "NrA%d	%d	%d	%d\n", &tmp, &source, &sink, &cost);
		G[source - 1][sink - 1] = cost;
		D[source - 1][sink - 1] = cost;
		P[source - 1][sink - 1] = source - 1;
	}

	return 0;
}

int countCost(int *solution) {
	int C = 0, i;
	int t = 0;

	for (i = 0; i < soLen; i++) {
		t = abs(solution[i]);
		if (solution[i] >= 0) {
			if (i == 0)
				C += D[depot - 1][tasks[t].head];
			else {
				if (solution[i - 1] >= 0)
					C += D[tasks[solution[i - 1]].tail][tasks[t].head];
				else
					C += D[tasks[-solution[i - 1]].head][tasks[t].head];
			}
			C += G[tasks[t].head][tasks[t].tail];
		}
		else {
			if (i == 0)
				C += D[depot - 1][tasks[t].tail];
			else {
				if (solution[i - 1] >= 0)
					C += D[tasks[solution[i - 1]].tail][tasks[t].tail];
				else
					C += D[tasks[-solution[i - 1]].head][tasks[t].tail];
			}
			C += G[tasks[t].tail][tasks[t].head];
		}
	}
	if (solution[i - 1] >= 0)
		C += D[tasks[solution[i - 1]].tail][depot - 1];
	else
		C += D[tasks[-solution[i - 1]].head][depot - 1];

	return C;
}

int insertDelta(int index_i, int index_t) {
	int *tmpSol = (int *)malloc(soLen * sizeof(int));
	memcpy(tmpSol, solution, soLen * sizeof(int));
	sinsSol(tmpSol, index_i, index_t);
	if (checkFill(tmpSol) != 0) {
		my_free(tmpSol);
		return MAXINT;
	}

	my_free(tmpSol);

	if (index_i == index_t || index_i == index_t + 1) {
		return 0;
	}

	int a = index_i - 1;
	int b = index_i;
	int c = index_i + 1;
	int d = index_t;
	int e = index_t + 1;
	int ab, bc, de, ac, db, be;

	if (a >= 0) {
		ab = singleCost(index[solution[a]], index[solution[b]]);
		if (c < soLen)
			ac = singleCost(index[solution[a]], index[solution[c]]);
		else
			ac = singleCost(index[solution[a]], 0);
	}
	else {
		ab = singleCost(0, index[solution[b]]);
		ac = singleCost(0, index[solution[c]]);
	}

	if (c < soLen)
		bc = singleCost(index[solution[b]], index[solution[c]]);
	else
		bc = singleCost(index[solution[b]], 0);

	if (e < soLen) {
		de = singleCost(index[solution[d]], index[solution[e]]);
		be = singleCost(index[solution[b]], index[solution[e]]);
	}
	else {
		de = singleCost(index[solution[d]], 0);
		be = singleCost(index[solution[b]], 0);
	}

	db = singleCost(index[solution[d]], index[solution[b]]);

	int k = (ac + db + be) - (ab + bc + de);
	return k;
}

int singleCost(int a, int b) {
	assert(abs(a) >= 0 && abs(a) <= taskNum);
	assert(abs(b) >= 0 && abs(b) <= taskNum);

	if (a >= 0 && b >= 0)
		return D[tasks[a].tail][tasks[b].head];
	else if (a < 0 && b >= 0)
		return D[tasks[-a].head][tasks[b].head];
	else if (a >= 0 && b < 0)
		return D[tasks[a].tail][tasks[-b].tail];
	else
		return D[tasks[-a].head][tasks[-b].tail];
}


int checkFill(int *tmp) {
	int i, load = 0;

	//printf("Q: %d\n", Q);
	for (i = 0; i < soLen; i++) {
		if (tmp[i] >= soLen) {
			printf("solution error!\n");
			return -1;
		}
		if (index[tmp[i]] == 0) {
			if (load > Q)
				return -1;
			//assert(load <= Q && Q == 1437);
			//printf("load is %d\n", load);
			load = 0;
		}
		else
			load += tasks[abs(index[tmp[i]])].req;
	}
	if (load > Q)
		return -1;
	else
		return 0;
}

int checkSol(int stCost, FILE *out) {
	int i, cost = 0;
	if (checkFill(solution) != 0)
		return 0;
	for (i = 0; i < soLen - 1; i++) {
		cost = cost + G[tasks[index[solution[i]]].head][tasks[index[solution[i]]].tail] + singleCost(index[solution[i]], index[solution[i + 1]]);
	}
	cost = cost + G[tasks[index[solution[i]]].head][tasks[index[solution[i]]].tail] + singleCost(0, index[solution[0]]) + singleCost(index[solution[i]], 0);

	printf("\nCost just counted is %d, and the stCost is %d.\n", cost, stCost);
	fprintf(out, "\nCost just counted is %d, and the stCost is %d.\n", cost, stCost);

	if (cost == stCost)
		return 1;
	else
		return 0;
}
