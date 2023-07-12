
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
#include <chrono>
// Procedures could be modified by the user

using namespace std;

double simtime = 0; //current simulation time, in second
unordered_map<int, int> optimal_lane_set; //secid-laneid
unordered_map<int, double> link_flw;
unordered_map<int, int> link_list;
unordered_map<int, int> from_turn;
unordered_map<int, int> to_turn;
unordered_map<int, double> turn_pert;

// set variables for energy consumptions
unordered_map<int, double> ice_base;
unordered_map<int, double> pev_base;
unordered_map<int, double> phev_base1;
unordered_map<int, double> phev_base2;
unordered_map<int, double> hfcv_base;

unordered_map<int, double> veh_tt;
unordered_map<int, double> ice_eng;
unordered_map<int, double> bev_eng;
unordered_map<int, double> phev_eng1;
unordered_map<int, double> phev_eng2;
unordered_map<int, double> hfcv_eng;
unordered_map<int, double> vtyp_eng;

int N_lnk, N_turn, N_typ, N_step;
ifstream flink, fturn;
double** lnk_flow, ** trn_per;

std::chrono::system_clock::time_point start_time;
std::chrono::system_clock::time_point end_time;

void printDebugLog(string s) {
	AKIPrintString(("## Debug log ##: " + s).c_str());
}

int AAPILoad()
{
	start_time = std::chrono::system_clock::now();
	srand((uint32_t)time(NULL));

	printDebugLog("in AAPILoad");
	return 0;
}

void Emission(double spd, double grade, double acc, double& E_i, double& E_e, double& E_p1, double& E_p2, double& E_f)
{
	double m = 1928, g = 9.8006, cr = 1.75, c1 = 0.0328, c2 = 4.575, rhoa = 1.2256, etad = 0.80;
	double Af = 2.73, cd = 0.34, beta = 0.93, alpha = 0.2, ch = 0.95;
	double va = 32, Pa = 2.5, Pb = 5, P_aux = 0.7, P_idle = 1.0;
	double a0 = 0.000341, a1 = 0.0000583, a2 = 0.000001;
	double eng[5]; // ICE, PEV, HPEV1, HPEV2, HFEV
	double rst, P_i, P_w, P_b, P_f, P_h;
	double P_batt, P_fuel;
	spd = spd / 3.6;

	// VT-CPFM model for gasoline vehicles  %%% some issue about the fuel consumption calculation, shall be around 1.5 L
	rst = rhoa * cd * ch * Af * pow(spd, 2) / 25.92 + g * m * cr * (c1 * spd + c2) / 1000 + g * m * grade;
	P_i = ((rst + 1.04 * m * acc) / (3600 * etad)) * spd;
	eng[0] = a0;
	if (P_i > 0) eng[0] = a0 + a1 * P_i + a2 * pow(P_i, 2);

	// another try for VT-CPFM model
	eng[0] = a0;
	if (P_w > 0) eng[0] = a0 + a1 * P_w / 1000 + a2 * pow(P_w / 1000, 2);

	// CPEM model for EVs
	double gs, gc;
	double eta_dl = 0.92, eta_em = 0.91, eta_rb = 0;
	gs = sqrt(1 / (1 + pow(grade, 2)));
	gc = sqrt(1 - pow(gs, 2));
	P_w = (m * acc + m * g * gc * cr * (c1 * spd + c2) / 1000 + rhoa * Af * cd * pow(spd, 2) / 2 + m * g * gs) * spd;
	if (P_w < 0) P_w = 0;
	eng[1] = P_w / (eta_dl * eta_em);
	if (eng[1] < 0) {
		if (acc < 0) eta_rb = 1 / exp(0.0411 / abs(acc));
		eng[1] = eng[1] * eta_rb;
	}
	eng[1] = eng[1] / 3600;

	// PHEV energy model
	double P_max = 100 * 0.3;			//max engine motor power
	if (eng[1] > P_max) {
		P_h = (eng[1] - P_max);
		eng[2] = P_max;								// electric usage for PHEV
		eng[3] = a0 + a1 * P_h + a2 * pow(P_h, 2);  // fuel usage for PHEV
	}
	else {
		eng[2] = eng[1];
		eng[3] = 0;
	}

	// HFCV energy model
	if (P_w <= Pa) {
		P_batt = P_w + P_aux;
		P_fuel = P_idle;
	}
	else {
		if (P_w > Pa && P_w <= Pb) {
			if (spd * 3.6 <= va) {
				P_batt = P_w + P_aux;
				P_fuel = P_idle;
			}
			else {
				P_batt = P_idle + P_aux;
				P_fuel = P_w * beta;
			}
		}
		else {
			P_batt = P_w * alpha + P_aux;
			P_fuel = P_w * beta;
		}
	}
	eng[4] = P_batt + P_fuel;

	E_i = eng[0];
	E_e = eng[1];
	E_p1 = eng[2];
	E_p2 = eng[3];
	E_f = eng[4];
	//AKIPrintString(("Energy 00000000: " + to_string(grade) + ", ICE = " + to_string(E_e) + ", PHEV = " + to_string(E_p1) + ", PHEV 2 = " + to_string(E_p2) + ", HFCV = " + to_string(E_f)).c_str());
}

int AAPIInit()
{
	printDebugLog("in AAPInit");
	ANGConnEnableVehiclesInBatch(true);

	int veh_type_nb = AKIVehGetNbVehTypes();
	printDebugLog("total number of vehicle types is: " + to_string(veh_type_nb));

	N_typ = AKIVehGetNbVehTypes();
	N_step = 24;
	lnk_flow = new double* [int(N_step)];
	trn_per = new double* [int(N_step)];
	int secnb = AKIInfNetNbSectionsANG();			// obtain the number of links in the network

	// update the initial link energy consumption
	for (int i = 0; i < secnb; i++) {
		int secid = AKIInfNetGetSectionANGId(i);
		A2KSectionInf secinf = AKIInfNetGetSectionANGInf(secid);
		double spd = secinf.speedLimit / 3.6;
		double grade = secinf.slopePercentages[0] / 100;

		Emission(spd, grade, 0.0, ice_base[secid], pev_base[secid], phev_base1[secid], phev_base2[secid], hfcv_base[secid]);

		ice_base[secid] = ice_base[secid] * secinf.length / spd;
		pev_base[secid] = pev_base[secid] * secinf.length / spd;
		phev_base1[secid] = phev_base1[secid] * secinf.length / spd;
		phev_base2[secid] = phev_base2[secid] * secinf.length / spd;
		hfcv_base[secid] = hfcv_base[secid] * secinf.length / spd;
	}


	// update the initial link travel time
	double link_tt, link_spd;
	for (int i = 0; i < secnb; i++) {
		int secid = AKIInfNetGetSectionANGId(i);
		if (secid > 0) {
			A2KSectionInf& secinf = AKIInfNetGetSectionANGInf(secid);
			link_tt = secinf.length / (secinf.speedLimit / 3.6);
			AKISetSectionUserDefinedCost(secid, link_tt);
			AKISetSectionUserDefinedCost2(secid, link_tt);
			AKISetSectionUserDefinedCost3(secid, link_tt);
		}
	}


	printDebugLog("total number of sections is: " + to_string(secnb));

	return 0;
}

int AAPISimulationReady()
{
	AKIPrintString("\tAAPISimulationReady");
	return 0;
}

int AAPIManage(double time, double timeSta, double timTrans, double acicle)
{
	return 0;
}

int AAPIPostManage(double time, double timeSta, double timTrans, double acicle)
{

	return 0;
}

int AAPIFinish()
{
	//AKIPrintString("\tFinish");

	end_time = std::chrono::system_clock::now();
	printDebugLog("Total time elapsed: " + to_string(std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()) + " ms.  ###");

	return 0;
}

int AAPIUnLoad()
{
	//AKIPrintString("UNLOAD");
	return 0;
}

int AAPIPreRouteChoiceCalculation(double time, double timeSta)
{
	// update link-level traffic conditions at every 1 minutes, applied for dynamic routing with minimum travel time
	// can be updated wtih energy consumption by changing the cost to vehicle energy usage
	double link_tt, link_spd;
	int secnb = AKIInfNetNbSectionsANG();
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
			AKISetSectionUserDefinedCost2(secid, link_tt);
			AKISetSectionUserDefinedCost3(secid, link_tt);
		}
	}
	printDebugLog("AAPIPreRouteChoiceCalculation is called.");
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

int AAPIVehicleStartParking(int idveh, int idsection, double time)
{
	return 0;
}