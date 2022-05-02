#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <string>
#include <stack>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "DS_timer.h"
#include "DS_definitions.h"

using namespace std;

// cluster ID
#define UNCLASSIFIED -1 // default cluster ID
#define NOISE -2

double EPS = 0.0; // epsilon (radius)
int MIN_PTRS = 0; // the minimum number of points in EPS to be a CORE point

int NUM_THREADS = 0;	// the number of threads


struct Point {
public:
	double x, y;
	int clusterID;
	char dummy[40];	// dummy

	Point(int id = UNCLASSIFIED) {
		x = y = 0.0;
		clusterID = id;
	}
	
	bool isSameGroup(const Point& other) const {
		if (this->clusterID != other.clusterID) {
			return false;
		}
		return true;
	}

	void operator=(const Point& other) {
		this->x = other.x;
		this->y = other.y;
		this->clusterID = other.clusterID;
	}
};


void dbscanSerial(vector<Point>& dst);
void dbscanParallel1(vector<Point>& dst);
void dbscanParallel2(vector<Point>& dst);

void readData(const string FILENAME, vector<Point>& src, vector<Point>& serial, vector<Point>& parallel1, vector<Point>& parallel2, unsigned int& nPoints);
void writeData(const string FILENAME, vector<Point>& result, unsigned int nPoints, int flag);

int countGroup(const vector<Point>& dst);

bool compareResult(const vector<Point>& serial, const vector<Point>& parallel);
int getErr(const vector<Point>& serial, const vector<Point>& parallel);

int main(int argc, char** argv) {
	if (argc < 5) {
		printf("This Program requires four arguments\n");
		printf("Usage: Extuction_file [InputData] [EPS] [MIN_PTRS] [NUM_THREADS]\n");
		return -1;
	}

	// init
	DS_timer timer(3);
	timer.setTimerName(0, (char*)"Serial");
	timer.setTimerName(1, (char*)"Parallel1");
	timer.setTimerName(2, (char*)"Parallel2");

	EPS = stod(argv[2]);
	MIN_PTRS = atoi(argv[3]);
	NUM_THREADS = atoi(argv[4]);

	printf("EPS: %.2lf, MIN_PTRS: %d, NUM_THREADS: %d\n\n", EPS, MIN_PTRS, NUM_THREADS);
	
	// src(source) and dst(serial, parallel)
	vector<Point> src;
	vector<Point> serial;
	vector<Point> parallel1;
	vector<Point> parallel2;

	//*** 1. load data **//
	string FILENAME = argv[1];
	string FILEPATH = "./Release/data/" + FILENAME; // INPUT DATA FILE path
	unsigned int nPoints = 0; // the number of points in dataset
	readData(FILEPATH, src, serial, parallel1, parallel2, nPoints);

	printf("The number of points : %d\n", nPoints);
	
	//*** 2. serial code ***//
	timer.onTimer(0);
	dbscanSerial(serial);
	timer.offTimer(0);
	
	int groupCntSerial = countGroup(serial);
	printf("[serial] the number of Group (with NOISE) : %d\n", groupCntSerial);

	
	//*** 3-1. parallel code ***//
	timer.onTimer(1);
	dbscanParallel1(parallel1);
	timer.offTimer(1);

	int groupCntParallel1 = countGroup(parallel1);
	printf("[parallel1] the number of Group (with NOISE) : %d\n", groupCntParallel1);

	//*** 3-2. parallel code ***//
	timer.onTimer(2);
	dbscanParallel2(parallel2);
	timer.offTimer(2);

	int groupCntParallel2 = countGroup(parallel2);
	printf("[parallel2] the number of Group (with NOISE) : %d\n\n", groupCntParallel2);
	
	//*** 4. write result ***//
	writeData(FILEPATH, serial, nPoints, 0);
	writeData(FILEPATH, parallel1, nPoints, 1);
	writeData(FILEPATH, parallel2, nPoints, 2);

	//*** 5. compare result ***//
	if (groupCntSerial==groupCntParallel1) {
		cout << "[serial & parallel1] result: GREAT! \n";
	}
	else {
		cout << "[serial & parallel1] result: there are difference between serial and parallel1, check them !\n";
	}
	if (groupCntSerial == groupCntParallel2) {
		cout << "[serial & parallel2] result: GREAT! \n";
	}
	else {
		cout << "[serial & parallel2] result: there are difference between serial and parallel2, check them !\n";
	}

	cout << "\n";

	// print timer
	timer.printTimer();
	
	vector<Point>().swap(src);
	vector<Point>().swap(serial);
	vector<Point>().swap(parallel1);
	vector<Point>().swap(parallel2);

	return 0;
}


// Read data and then, assign to src and dst
void readData(const string FILENAME, vector<Point>& src, vector<Point>& serial, vector<Point>& parallel1, vector<Point>& parallel2, unsigned int& nPoints) {
	FILE* stream = NULL;
	fopen_s(&stream, FILENAME.c_str(), "r");
	if (stream == NULL) {
		cerr << "file read err\n" << endl;
		exit(1);
	}

	fscanf_s(stream, "%u\n", &nPoints);

	src.resize(nPoints);
	serial.resize(nPoints);
	parallel1.resize(nPoints);
	parallel2.resize(nPoints);


	for (int i = 0; i < nPoints; ++i) {
		fscanf_s(stream, "%lf %lf\n", &(src[i].x), &(src[i].y));
		src[i].clusterID = UNCLASSIFIED;
		serial[i] = src[i];
		parallel1[i] = src[i];
		parallel2[i] = src[i];
	}
	fclose(stream);
}

// write both serial and parallel
void writeData(const string FILENAME, vector<Point>& result, unsigned int nPoints, int flag) {
	FILE* out = NULL;
	
	if (flag == 0) {
		fopen_s(&out, (FILENAME.substr(0, FILENAME.rfind('.')) + "_serial.txt").c_str(), "w");	// ex) moon.txt ==> moon_serial.txt
		if (out == NULL) {
			cerr << "file write err\n" << endl;
			exit(1);
		}
	}
	else if (flag == 1) {
		fopen_s(&out, (FILENAME.substr(0, FILENAME.rfind('.')) + "_parallel1.txt").c_str(), "w");	// ex) moon.txt ==> moon_parallel1.txt
		if (out == NULL) {
			cerr << "file write err\n" << endl;
			exit(1);
		}
	}
	else if (flag == 2) {
		fopen_s(&out, (FILENAME.substr(0, FILENAME.rfind('.')) + "_parallel2.txt").c_str(), "w");	// ex) moon.txt ==> moon_parallel2.txt
		if (out == NULL) {
			cerr << "file write err\n" << endl;
			exit(1);
		}
	}

	fprintf(out, "%u\n", nPoints);

	Point tmp;
	for (int i = 0; i < nPoints; ++i) {
		tmp = result[i];
		fprintf(out, "%.2lf %.2lf %d\n", tmp.x, tmp.y, tmp.clusterID);
	}
	fclose(out);
}


int countGroup(const vector<Point>& dst) {
	unordered_set<int> groupCnt;
	for (int i = 0; i < dst.size(); ++i) {
		groupCnt.insert(dst[i].clusterID);
	}
	return groupCnt.size();
}

// dbscan - serial
void dbscanSerial(vector<Point>& dst) {
	int cID = 0; // cluster number

	for (int i = 0; i < dst.size(); ++i)
	{
		if (dst[i].clusterID != UNCLASSIFIED)
		{
			continue;
		}

		int cnt = 0; // ÀÌ¿ô °¹¼ö
		for (int j = 0; j < dst.size(); ++j)
		{
			if (i != j && sqrt(pow(dst[i].x - dst[j].x, 2) + pow(dst[i].y - dst[j].y, 2)) <= EPS)
			{
				cnt++;
			}
		}

		// when it is not core
		if (cnt < MIN_PTRS)
		{
			dst[i].clusterID = NOISE;
		}
		// if it is core
		else
		{
			stack<int> parallelstack;
			parallelstack.push(i);
			int cur = 0;
			while (!parallelstack.empty())
			{
				cur = parallelstack.top();
				parallelstack.pop();

				dst[cur].clusterID = cID;
				for (int j = 0; j < dst.size(); ++j)
				{
					if (cur != j && (dst[j].clusterID == UNCLASSIFIED || dst[j].clusterID == NOISE) && sqrt(pow(dst[cur].x - dst[j].x, 2) + pow(dst[cur].y - dst[j].y, 2)) <= EPS)
					{
						parallelstack.push(j);
					}
				}
			}
			++cID;
		}
	}
}

// dbscan - parallel 1
void dbscanParallel1(vector<Point>& dst) {
	vector<omp_lock_t> pointLock(dst.size()); // locks for each point

	#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < pointLock.size(); ++i) {
		omp_init_lock(&pointLock[i]);
	}
	// implicit barrier here

	//vector<unordered_set<int>> globalMap; // for merging local cluster to make global cluster

	//omp_lock_t globalMapLock;
	//omp_init_lock(&globalMapLock);

	int cID = 0; // cluster ID
	int maxCID = 0; // maximum cID;
	//#pragma omp parallel num_threads(NUM_THREADS) private(cID) shared(globalMap, globalMapLock, maxCID)
	#pragma omp parallel num_threads(NUM_THREADS) private(cID) shared(maxCID)
	{
		cID = omp_get_thread_num();

		#pragma omp for 
		for (int i = 0; i < dst.size(); ++i) {
			omp_set_lock(&pointLock[i]);
			Point curTmp = dst[i];
			omp_unset_lock(&pointLock[i]);

			omp_set_lock(&pointLock[i]);
			if (dst[i].clusterID != UNCLASSIFIED) {
				omp_unset_lock(&pointLock[i]);
				continue;
			}
			omp_unset_lock(&pointLock[i]);

			stack<int> neighbors;
			for (int j = 0; j < dst.size(); ++j) {
				if (i != j && sqrt(pow(curTmp.x - dst[j].x, 2) + pow(curTmp.y - dst[j].y, 2)) <= EPS) {
					neighbors.push(j);
				}
			}

			// when it is not core
			if (neighbors.size() < MIN_PTRS) {
				omp_set_lock(&pointLock[i]);
				dst[i].clusterID = NOISE;
				omp_unset_lock(&pointLock[i]);
			}
			// if it is core
			else {
				omp_set_lock(&pointLock[i]);
				dst[i].clusterID = cID;
				omp_unset_lock(&pointLock[i]);

				/*unordered_set<int> localMap;
				localMap.insert(cID);*/

				int cur = 0;
				while (!neighbors.empty()) {
					cur = neighbors.top();
					neighbors.pop();

					omp_set_lock(&pointLock[cur]);
					if (!(dst[cur].clusterID == UNCLASSIFIED || dst[cur].clusterID==NOISE)) {
						omp_unset_lock(&pointLock[cur]);
						continue;
					}
					dst[cur].clusterID = cID;
					omp_unset_lock(&pointLock[cur]);

					for (int j = 0; j < dst.size(); ++j) {
						omp_set_lock(&pointLock[j]);
						if (cur != j && (dst[j].clusterID == UNCLASSIFIED || dst[j].clusterID == NOISE) && sqrt(pow(dst[cur].x - dst[j].x, 2) + pow(dst[cur].y - dst[j].y, 2)) <= EPS) {
							neighbors.push(j);
						}
						//if (cur != j && !(dst[j].clusterID == cID || dst[j].clusterID == UNCLASSIFIED || dst[j].clusterID == NOISE) && sqrt(pow(dst[cur].x - dst[j].x, 2) + pow(dst[cur].y - dst[j].y, 2)) <= EPS) {
						//	localMap.insert(dst[j].clusterID);	// when meet other local group nearby
						//}
						omp_unset_lock(&pointLock[j]);
					}
				}
				//omp_set_lock(&globalMapLock);
				//if (localMap.size() > 1) { // when meet other local group nearby
				//	globalMap.push_back(localMap);
				//}
				//omp_unset_lock(&globalMapLock);
				
				cID += omp_get_num_threads();
			}
		}
		#pragma omp critical
		{
			maxCID = max(maxCID, cID);	// update Maximum Cluster ID
		}
	}

	vector<unordered_set<int>> trueMap(maxCID, unordered_set<int>());
	vector<omp_lock_t> trueMapLock(maxCID);
	
	#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < maxCID; ++i) {
		omp_init_lock(&trueMapLock[i]);
	}

	#pragma omp parallel for num_threads(NUM_THREADS) shared(trueMap, trueMapLock)
	for (int i = 0; i < dst.size(); ++i) {
		if (dst[i].clusterID == NOISE) {
			continue;
		}
		for (int j = 0; j < dst.size(); ++j) {
			if ((i != j) && (dst[j].clusterID != NOISE && dst[j].clusterID != dst[i].clusterID) && (sqrt(pow(dst[i].x - dst[j].x, 2) + pow(dst[i].y - dst[j].y, 2)) <= EPS)) {
				omp_set_lock(&trueMapLock[dst[i].clusterID]);
				trueMap[dst[i].clusterID].insert(dst[j].clusterID);
				omp_unset_lock(&trueMapLock[dst[i].clusterID]);
			}
		}
	}

	// DFS for grouping connected local cluster
	unordered_set<int> visited;
	vector<unordered_set<int>> group;	// global cluster group
	
	for (int i = 0; i < trueMap.size(); ++i) {
		if (visited.find(i) != visited.end() || trueMap[i].size()==0) {	// already visited OR there is no neighbor cluster
			
			continue;
		}
		visited.insert(i);	
		
		unordered_set<int> seen;	// local visited
		seen.insert(i);

		stack<int> stack;
		stack.push(i);
		int cur = 0;
		while (!stack.empty()) {
			cur = stack.top();
			stack.pop();
			for (unordered_set<int>::iterator it = trueMap[cur].begin(); it != trueMap[cur].end(); ++it) {
				if (seen.find(*it) == seen.end()) { // not visited
					seen.insert(*it);
					stack.push(*it);
				}
			}
		}

		for (unordered_set<int>::iterator it = seen.begin(); it != seen.end(); ++it) {
			visited.insert(*it);
		}
		group.push_back(seen);
	}

	// merge : local neighbor cluster group -> global cluster
	unordered_set<int> outGroupID; // cID whose has no neighbor

	#pragma omp parallel for num_threads(NUM_THREADS) shared(group)
	for (int i = 0; i < dst.size(); ++i) {
		if (dst[i].clusterID == NOISE) {
			continue;
		}
		bool merged = false;
		for (int j = 0; j < group.size(); ++j) {
			if (group[j].find(dst[i].clusterID) != group[j].end()) {	// when find global group
				dst[i].clusterID = j;
				merged = true;
				break;
			}
		}
		if (merged == false) {
			outGroupID.insert(dst[i].clusterID);
		}
	}

	// high local group Id to global group
	unordered_map<int, int> local2global;
	int idx = group.size();
	for (unordered_set<int>::iterator it = outGroupID.begin(); it != outGroupID.end(); ++it) {
		local2global[*it] = idx++;
	}

	
	#pragma omp parallel for num_threads(NUM_THREADS) shared(outGroupID)
	for (int i = 0; i < dst.size(); ++i) {
		if (outGroupID.find(dst[i].clusterID) != outGroupID.end()) {
			dst[i].clusterID = local2global[dst[i].clusterID];
		}
	}


	#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < maxCID; ++i) {
		omp_destroy_lock(&trueMapLock[i]);
	}

	#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < pointLock.size(); ++i) {
		omp_destroy_lock(&pointLock[i]);
	}
}


// dbscan - parallel 2
void dbscanParallel2(vector<Point>& dst)
{
	int cID = 0; // cluster number

	int* start = new int[NUM_THREADS];
	int* end = new int[NUM_THREADS];

	omp_lock_t writelock;
	omp_init_lock(&writelock);

	for (int i = 0; i < NUM_THREADS; i++)
	{
		start[i] = i * dst.size() / NUM_THREADS;
		end[i] = (i + 1) * dst.size() / NUM_THREADS;
		if (i == NUM_THREADS - 1)
		{
			end[NUM_THREADS - 1] = dst.size();
		}
	}

	for (int i = 0; i < dst.size(); ++i)
	{
		if (dst[i].clusterID != UNCLASSIFIED)
		{
			continue;
		}

		int cnt = 0; // ÀÌ¿ô °¹¼ö
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			int tID = omp_get_thread_num();
			for (int j = start[tID]; j < end[tID]; j++)
			{
				if (i != j && sqrt(pow(dst[i].x - dst[j].x, 2) + pow(dst[i].y - dst[j].y, 2)) <= EPS)
				{
					cnt++;
				}
			}
		}

		// when it is not core
		if (cnt < MIN_PTRS)
		{
			dst[i].clusterID = NOISE;
		}
		// if it is core
		else
		{
			stack<int> parallelstack;
			parallelstack.push(i);
			int cur = 0;

			while (!parallelstack.empty())
			{
				cur = parallelstack.top();
				parallelstack.pop();
				dst[cur].clusterID = cID;

				#pragma omp parallel num_threads(NUM_THREADS)
				{
					int tID = omp_get_thread_num();
					for (int j = start[tID]; j < end[tID]; j++)
					{
						if (cur != j && (dst[j].clusterID == UNCLASSIFIED || dst[j].clusterID == NOISE) && sqrt(pow(dst[cur].x - dst[j].x, 2) + pow(dst[cur].y - dst[j].y, 2)) <= EPS)
						{
							omp_set_lock(&writelock);
							parallelstack.push(j);
							omp_unset_lock(&writelock);
						}
					}
				}
			}
			++cID;
		}
	}
	omp_destroy_lock(&writelock);

	delete[] start;
	delete[] end;
}
