
#include "AKIProxie.h"
#include "CIProxie.h"
#include "ANGConProxie.h"
#include "AAPI.h"
#include <stdio.h>
#include <unordered_map>
#include <map>
#include <string>
#include <deque>
#include <queue>
#include <time.h>
#include <set>
#include <random>
#include <cmath>
#include <fstream>
// Procedures could be modified by the user

using namespace std;

double simtime = 0; //current simulation time, in second
unordered_map<int, int> optimal_lane_set; //secid-laneid
unordered_map<int, double> link_flw;
unordered_map<int, int> link_list;
unordered_map<int, int> from_turn;
unordered_map<int, int> to_turn;
unordered_map<int, double> turn_pert;
int N_lnk, N_turn, N_typ, N_step;
ifstream flink, fturn;
double** lnk_flow, ** trn_per;
const string PROJECT_DIR = "D:\\program\\Aimsun\\digital_twin_aimsun";

void printDebugLog(string s) {
	AKIPrintString(("## Debug log ##: " + s).c_str());
}

int AAPILoad()
{
	srand((uint32_t)time(NULL));

	printDebugLog("in AAPILoad");

	return 0;
}

int AAPIInit()
{

	printDebugLog("in AAPILoad");
	ANGConnEnableVehiclesInBatch(true);

	int veh_type_nb = AKIVehGetNbVehTypes();
	printDebugLog("total number of vehicle types is: " + to_string(veh_type_nb));

	int nslice = AKIStateDemandGetNumSlices(1);
	AKIPrintString(("Number of Slices: " + to_string(nslice)).c_str());
	nslice = AKIStateDemandGetNumSlices(2);
	AKIPrintString(("Number of Slices: " + to_string(nslice)).c_str());
	nslice = AKIStateDemandGetNumSlices(3);
	AKIPrintString(("Number of Slices: " + to_string(nslice)).c_str());
	nslice = AKIStateDemandGetNumSlices(4);
	AKIPrintString(("Number of Slices: " + to_string(nslice)).c_str());

	N_typ = AKIVehGetNbVehTypes();
	N_step = 24;
	lnk_flow = new double* [int(N_step)];
	trn_per = new double* [int(N_step)];


	// read the list of links and turns for traffic state update
	int lnk0, lnk1, lnk2;
	double flw, pert;
	N_lnk = 0;
	flink.open(PROJECT_DIR + "\\demand_data\\Link.txt");
	while (flink >> lnk0 >> flw) {
		link_list[N_lnk] = lnk0;
		link_flw[N_lnk] = flw;
		N_lnk++;
	}
	flink.close();
	printDebugLog("Total number links found in Link.txt is " + to_string(N_lnk - 1));

	N_turn = 0;
	fturn.open(PROJECT_DIR + "\\demand_data\\Turn.txt");
	while (fturn >> lnk1 >> lnk2 >> pert) {
		from_turn[N_turn] = lnk1;
		to_turn[N_turn] = lnk2;
		turn_pert[N_turn] = pert;
		N_turn++;
	}
	fturn.close();
	printDebugLog("Total number turns found in Turn.txt is " + to_string(N_turn - 1));

	for (int i = 0; i < N_step; i++) {
		lnk_flow[i] = new double[N_lnk];
		trn_per[i] = new double[N_turn];
	}

	return 0;
}


int AAPIManage(double time, double timeSta, double timTrans, double acicle)
{
	simtime = AKIGetCurrentSimulationTime();
	int N_hist, N_pre;


	// update traffic state at a specific time interval (5 min) with real-world traffic
	if (int(simtime * 10) % 3000 == 0) {

		// update the real-world conditions: read link flow and turn percentages in the past N_steps (5 minute interval)
		int lnk0, lnk1, lnk2;
		double flw, pert;
		N_hist = int(N_step / 2);
		//flink.open("\\demand_data\\Link_update.txt");
		flink.open("\\demand_data\\Link.txt");
		for (int j = 0; j < N_lnk; j++)
			for (int i = 0; i < N_hist; i++) {
				flink >> lnk_flow[i][j];
			}
		flink.close();

		N_turn = 0;
		//fturn.open("\\demand_data\\Turn_update.txt");
		fturn.open("\\demand_data\\Turn.txt");
		for (int j = 0; j < N_turn; j++)
			for (int i = 0; i < N_hist; i++) {
				flink >> trn_per[i][j];
			}
		fturn.close();

		// predict the road traffic conditions
		// current use moving average (will be replaced with the trained machine learning models)
		double tmpf, tmpp;
		for (int j = 0; j < N_lnk; j++) {
			for (int t = N_hist; t < N_step; t++) {
				tmpf = 0;
				for (int k = 1; k <= 5; k++) {
					tmpf += lnk_flow[t - k][j];
				}
				lnk_flow[t][j] = tmpf / 5;
			}
		}
		for (int j = 0; j < N_turn; j++) {
			for (int t = N_hist; t < N_step; t++) {
				tmpp = 0;
				for (int k = 1; k <= 5; k++) {
					tmpf += trn_per[t - k][j];
				}
				trn_per[t][j] = tmpp / 5;
			}
		}


		// update the link and turn information
		for (int t = N_hist; t < N_step; t++) {
			for (int i = 0; i < N_lnk; i++) {
				AKIStateDemandSetDemandSection(link_list[i], 1, t + 1, lnk_flow[t][i]);
				printDebugLog("Update some demand to value: " + to_string(lnk_flow[t][i]));
			}
			for (int i = 0; i < N_turn; i++) {
				AKIStateDemandSetTurningPercentage(from_turn[i], to_turn[i], 1, t + 1, trn_per[t][i]);
			}
		}
	}

	// update link-level traffic conditions at every 1 minutes, applied for dynamic routing with minimum travel time
	// can be updated wtih energy consumption by changing the cost to vehicle energy usage
	double link_tt, link_spd;
	if (int(simtime * 10) % 600 == 0) {
		int secnb = AKIInfNetNbSectionsANG();
		//printDebugLog("The number of sections from API is: " + to_string(secnb));
		for (int i = 0; i < secnb; i++) {
			int secid = AKIInfNetGetSectionANGId(i);
			A2KSectionInf secinf = AKIInfNetGetSectionANGInf(secid);
			if (secid > 0) {
				int nbveh = AKIVehStateGetNbVehiclesSection(secid, true);
				if (nbveh > 0) {
					link_spd = 0;
					for (int k = 0; k < nbveh; k++) {
						InfVeh vehinf = AKIVehStateGetVehicleInfSection(secid, k);
						link_spd += vehinf.CurrentSpeed / 3.6;						// currently use all vehicles, can change it to CV only
					}
					link_spd = link_spd / nbveh;
					if (link_spd == 0) link_spd = 0.5;
					link_tt = secinf.length / link_spd;
				}
				else {
					link_tt = secinf.length / (secinf.speedLimit / 3.6);
				}
				AKISetSectionUserDefinedCost(secid, link_tt);
			}
		}
	}

	// rewind the simulation at every one hour, i.e., only predict traffic in one hour
	if (simtime == 3600) {
		ANGSetSimulationOrder(2, 0);
	}

	return 0;
}

int AAPIPostManage(double time, double timeSta, double timTrans, double acicle)
{

	return 0;
}

int AAPIFinish()
{
	//AKIPrintString("\tFinish");

	return 0;
}

int AAPIUnLoad()
{
	//AKIPrintString("UNLOAD");
	return 0;
}

int AAPIPreRouteChoiceCalculation(double time, double timeSta)
{
	//AKIPrintString("\tPreRouteChoice Calculation");
	return 0;
}

int AAPIEnterVehicle(int idveh, int idsection)
{
	return 0;
}

int AAPIExitVehicle(int idveh, int idsection)
{
	return 0;
}

int AAPIEnterVehicleSection(int idveh, int idsection, double atime)
{
	return 0;
}

int AAPIExitVehicleSection(int idveh, int idsection, double time)
{
	return 0;
}

int AAPIEnterPedestrian(int idPedestrian, int originCentroid)
{
	AKIPrintString("A Legion Pedestrian has entered the network");
	return 0;
}

int AAPIExitPedestrian(int idPedestrian, int destinationCentroid)
{
	AKIPrintString("A Legion Pedestrian has exited the network");
	return 0;
}