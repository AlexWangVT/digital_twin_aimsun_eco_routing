
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

void printDebugLog(string s) {
	AKIPrintString(("## Debug log ##: " + s).c_str());
}

int AAPILoad()
{
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

	printDebugLog("in AAPILoad");
	ANGConnEnableVehiclesInBatch(true);

	int veh_type_nb = AKIVehGetNbVehTypes();
	printDebugLog("total number of vehicle types is: " + to_string(veh_type_nb));

	int nslice = AKIStateDemandGetNumSlices(1);
	AKIPrintString(("Number of Slices: " + to_string(nslice)).c_str());

	N_typ = AKIVehGetNbVehTypes();
	N_step = 24;
	lnk_flow = new double* [int(N_step)];
	trn_per = new double* [int(N_step)];


	// read the list of links and turns for traffic state update
	int lnk0, lnk1, lnk2;
	double flw, pert;
	N_lnk = 0;
	flink.open("Link.txt");
	while (flink >> lnk0 >> flw) {
		link_list[N_lnk] = lnk0;
		link_flw[N_lnk] = flw;
		N_lnk++;
	}
	flink.close();

	N_turn = 0;
	fturn.open("Turn.txt");
	while (fturn >> lnk1 >> lnk2 >> pert) {
		from_turn[N_turn] = lnk1;
		to_turn[N_turn] = lnk2;
		turn_pert[N_turn] = pert;
		N_turn++;
	}
	fturn.close();

	for (int i = 0; i < N_step; i++) {
		lnk_flow[i] = new double[N_lnk];
		trn_per[i] = new double[N_turn];
	}


	// update the initial link energy consumption
	int secnb = AKIInfNetNbSectionsANG();			// obtain the number of links in the network
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


	printDebugLog("total number of sections is: " + to_string(secnb));

	return 0;
}


int AAPIManage(double time, double timeSta, double timTrans, double acicle)
{
	simtime = AKIGetCurrentSimulationTime();

	//{ // this block will add a large user defined cost to specific sections
	//	if (time - timTrans > 600) {
	//		AKISetSectionUserDefinedCost(34745, 9999);
	//		AKISetSectionUserDefinedCost(393017, 9999);
	//		//auto re4 = AKISetSectionUserDefinedCost(393158, 99999);

	//		AKISetSectionUserDefinedCost(118843, 9999);
	//		AKISetSectionUserDefinedCost(46367, 9999);
	//		AKISetSectionUserDefinedCost(393095, 9999);
	//	}
	//}


	int N_hist, N_pre;

	//double fac, sim_step = 0.5, eng_intval = 30;
	//double E_i, E_e, E_p1, E_p2, E_f, E_cost[5];
	//double P_gas = 4.5, P_elt = 0.12;		// price of gas and electricity

	//// update the link cost with the energy consumption at the interval of eng_interval
	//if (int(simtime * 10) % int(eng_intval * 10) == 0) {
	//	int secnb = AKIInfNetNbSectionsANG();			// obtain the number of links in the network
	//	for (int i = 0; i < secnb; i++) {
	//		int secid = AKIInfNetGetSectionANGId(i);
	//		if (secid > 0) {
	//			A2KSectionInf secinf = AKIInfNetGetSectionANGInf(secid);
	//			double grade = secinf.slopePercentages[0] / 100;
	//			int nbveh = AKIVehStateGetNbVehiclesSection(secid, true);
	//			for (int k = 0; k < nbveh; k++) {
	//				InfVeh vehinf = AKIVehStateGetVehicleInfSection(secid, k);
	//				StaticInfVeh vehstat = AKIVehGetVehicleStaticInfSection(secid, k);
	//				int type_id = AKIVehTypeGetIdVehTypeANG(vehinf.type);


	//				double acc = (vehinf.CurrentSpeed - vehinf.PreviousSpeed) / (3.6 * sim_step);		// acceleration in m/s^2
	//				Emission(vehinf.CurrentSpeed / 3.6, grade, acc, E_i, E_e, E_p1, E_p2, E_f);
	//				double fac;
	//				if (vehinf.CurrentSpeed > 0) {
	//					fac = secinf.length / (vehinf.CurrentSpeed / 3.6);
	//				}
	//				else {
	//					fac = secinf.length / (5.0 / 3.6);
	//				}
	//				E_cost[0] += E_i * fac;
	//				E_cost[1] += E_e * fac;
	//				E_cost[2] += E_p1 * fac;
	//				E_cost[3] += E_p2 * fac;
	//				E_cost[4] += E_f * fac;

	//				if (ice_eng[vehinf.idVeh] == 0) {
	//					veh_tt[vehinf.idVeh] = sim_step;
	//					ice_eng[vehinf.idVeh] = E_i * sim_step;
	//					bev_eng[vehinf.idVeh] = E_e * sim_step;
	//					phev_eng1[vehinf.idVeh] = E_p1 * sim_step;
	//					phev_eng2[vehinf.idVeh] = E_p2 * sim_step;
	//					hfcv_eng[vehinf.idVeh] = E_f * sim_step;
	//				}
	//				else {
	//					veh_tt[vehinf.idVeh] += sim_step;
	//					ice_eng[vehinf.idVeh] += E_i * sim_step;
	//					bev_eng[vehinf.idVeh] += E_e * sim_step;
	//					phev_eng1[vehinf.idVeh] += E_p1 * sim_step;
	//					phev_eng2[vehinf.idVeh] += E_p2 * sim_step;
	//					hfcv_eng[vehinf.idVeh] += E_f * sim_step;
	//				}

	//				vtyp_eng[vehinf.idVeh] = vehinf.type;
	//			}

	//			if (nbveh == 0) {
	//				E_cost[0] = ice_base[secid];
	//				E_cost[1] = pev_base[secid];
	//				E_cost[2] = phev_base1[secid];
	//				E_cost[3] = phev_base2[secid];
	//				E_cost[4] = hfcv_base[secid];
	//			}
	//			else {
	//				E_cost[0] = E_cost[0] / nbveh;
	//				E_cost[1] = E_cost[1] / nbveh;
	//				E_cost[2] = E_cost[2] / nbveh;
	//				E_cost[3] = E_cost[3] / nbveh;
	//				E_cost[4] = E_cost[4] / nbveh;
	//			}

	//			AKISetSectionUserDefinedCost(secid, E_cost[0] * P_gas / 3.8);		// Fuel consumption
	//			AKISetSectionUserDefinedCost2(secid, E_cost[1] * P_elt / 3600);	// EV Energy
	//			//AKISetSectionUserDefinedCost3(secid, E_cost[2] * P_elt / 3600 + E_cost[3] * P_gas / 3.8);	// PHEV Energy in $, 0.25/kWh, $6/gallon
	//			AKISetSectionUserDefinedCost3(secid, E_cost[4] * P_elt / 3600);	// HFCV Energy
	//		}
	//	}
	//}


	// update traffic state at a specific time interval (5 min) with real-world traffic
	if (int(simtime * 10) % 3000 == 0) {

		// update the real-world conditions: read link flow and turn percentages in the past N_steps (5 minute interval)
		int lnk0, lnk1, lnk2;
		double flw, pert;
		N_hist = int(N_step / 2);
		flink.open("Link_update.txt");
		for (int j = 0; j < N_lnk; j++)
			for (int i = 0; i < N_hist; i++) {
				flink >> lnk_flow[i][j];
			}
		flink.close();

		N_turn = 0;
		fturn.open("Turn_update.txt");
		for (int j = 0; j < N_turn; j++)
			for (int i = 0; i < N_hist; i++) {
				fturn >> trn_per[i][j];
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

	//// rewind the simulation at every one hour, i.e., only predict traffic in one hour
	//if (simtime == 3600) {
	//	ANGSetSimulationOrder(2, 0);
	//}

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
	double timTrans = AKIGetDurationTransTime();
	//{ // this block will add a large user defined cost to specific sections
	//	if (time - timTrans > 600) {
	//		printDebugLog("Setting user defined cost in AAPIPreRouteChoiceCalculation");
	//		AKISetSectionUserDefinedCost(34745, 9999);
	//		printDebugLog("After setting, the cost is: " + to_string(AKIGetSectionUserDefinedCost(34745)));
	//		AKISetSectionUserDefinedCost(393017, 9999);
	//		//auto re4 = AKISetSectionUserDefinedCost(393158, 99999);

	//		AKISetSectionUserDefinedCost(118843, 9999);
	//		AKISetSectionUserDefinedCost(46367, 9999);
	//		AKISetSectionUserDefinedCost(393095, 9999);
	//	}
	//}

	double fac, sim_step = 0.5, eng_intval = 30;
	double E_i, E_e, E_p1, E_p2, E_f, E_cost[5];
	double P_gas = 4.5, P_elt = 0.12;		// price of gas and electricity

	// update the link cost with the energy consumption at the interval of eng_interval
	int secnb = AKIInfNetNbSectionsANG();			// obtain the number of links in the network
	for (int i = 0; i < secnb; i++) {
		int secid = AKIInfNetGetSectionANGId(i);
		if (secid > 0) {
			A2KSectionInf secinf = AKIInfNetGetSectionANGInf(secid);
			double grade = secinf.slopePercentages[0] / 100;
			int nbveh = AKIVehStateGetNbVehiclesSection(secid, true);
			for (int k = 0; k < nbveh; k++) {
				InfVeh vehinf = AKIVehStateGetVehicleInfSection(secid, k);
				StaticInfVeh vehstat = AKIVehGetVehicleStaticInfSection(secid, k);
				int type_id = AKIVehTypeGetIdVehTypeANG(vehinf.type);


				double acc = (vehinf.CurrentSpeed - vehinf.PreviousSpeed) / (3.6 * sim_step);		// acceleration in m/s^2
				Emission(vehinf.CurrentSpeed / 3.6, grade, acc, E_i, E_e, E_p1, E_p2, E_f);
				double fac;
				if (vehinf.CurrentSpeed > 0) {
					fac = secinf.length / (vehinf.CurrentSpeed / 3.6);
				}
				else {
					fac = secinf.length / (5.0 / 3.6);
				}
				E_cost[0] += E_i * fac;
				E_cost[1] += E_e * fac;
				E_cost[2] += E_p1 * fac;
				E_cost[3] += E_p2 * fac;
				E_cost[4] += E_f * fac;

				if (ice_eng[vehinf.idVeh] == 0) {
					veh_tt[vehinf.idVeh] = sim_step;
					ice_eng[vehinf.idVeh] = E_i * sim_step;
					bev_eng[vehinf.idVeh] = E_e * sim_step;
					phev_eng1[vehinf.idVeh] = E_p1 * sim_step;
					phev_eng2[vehinf.idVeh] = E_p2 * sim_step;
					hfcv_eng[vehinf.idVeh] = E_f * sim_step;
				}
				else {
					veh_tt[vehinf.idVeh] += sim_step;
					ice_eng[vehinf.idVeh] += E_i * sim_step;
					bev_eng[vehinf.idVeh] += E_e * sim_step;
					phev_eng1[vehinf.idVeh] += E_p1 * sim_step;
					phev_eng2[vehinf.idVeh] += E_p2 * sim_step;
					hfcv_eng[vehinf.idVeh] += E_f * sim_step;
				}

				vtyp_eng[vehinf.idVeh] = vehinf.type;
			}

			if (nbveh == 0) {
				E_cost[0] = ice_base[secid];
				E_cost[1] = pev_base[secid];
				E_cost[2] = phev_base1[secid];
				E_cost[3] = phev_base2[secid];
				E_cost[4] = hfcv_base[secid];
			}
			else {
				E_cost[0] = E_cost[0] / nbveh;
				E_cost[1] = E_cost[1] / nbveh;
				E_cost[2] = E_cost[2] / nbveh;
				E_cost[3] = E_cost[3] / nbveh;
				E_cost[4] = E_cost[4] / nbveh;
			}

			AKISetSectionUserDefinedCost(secid, E_cost[0] * P_gas / 3.8);		// Fuel consumption
			AKISetSectionUserDefinedCost2(secid, E_cost[1] * P_elt / 3600);	// EV Energy
			AKISetSectionUserDefinedCost3(secid, E_cost[2] * P_elt / 3600 + E_cost[3] * P_gas / 3.8);	// PHEV Energy in $, 0.25/kWh, $6/gallon
			//AKISetSectionUserDefinedCost3(secid, E_cost[4] * P_elt / 3600);	// HFCV Energy
		}
	}

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