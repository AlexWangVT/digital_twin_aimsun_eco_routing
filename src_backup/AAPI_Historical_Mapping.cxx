
#include "AKIProxie.h"
#include "CIProxie.h"
#include "ANGConProxie.h"
#include "AAPI.h"
#include <stdio.h>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <deque>
#include <queue>
#include <ctime>
#include <set>
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
// Procedures could be modified by the user

using namespace std;

string PROJECT_DIR = "D:\\program\\Aimsun\\digital_twin_aimsun";
const bool GENERATE_HISTORICAL_DATA = false;
const double HISTORICAL_DATA_INTERVAL = 60.0; // seconds
const double PRICE_GAS = 3.5;
const double PRICE_ELECTRICITY = 0.25;
const double SOC_MIN = 0.204; // for PHEV
const double SOC_MAX = 0.88; // for PHEV

void printDebugLog(string s);

enum class VehicleType {
	UNKNOWN = 0,
	ICE = 393772,
	BEV = 393773,
	PHEV = 393774,
	HFCV = 393895,
	ICE_NONCAV = 393941,
	BEV_NONCAV = 393942,
	PHEV_NONCAV = 393943,
	HFCV_NONCAV = 393944
};
const VehicleType VEHICLE_TYPE_ICE = VehicleType::ICE;
const VehicleType VEHICLE_TYPE_BEV = VehicleType::BEV;
const VehicleType VEHICLE_TYPE_PHEV = VehicleType::PHEV;
const VehicleType VEHICLE_TYPE_HFCV = VehicleType::HFCV;
const VehicleType VEHICLE_TYPE_ICE_NONCAV = VehicleType::ICE_NONCAV;
const VehicleType VEHICLE_TYPE_BEV_NONCAV = VehicleType::BEV_NONCAV;
const VehicleType VEHICLE_TYPE_PHEV_NONCAV = VehicleType::PHEV_NONCAV;
const VehicleType VEHICLE_TYPE_HFCV_NONCAV = VehicleType::HFCV_NONCAV;

const vector<VehicleType> valid_vehicle_type_list{VEHICLE_TYPE_ICE, VEHICLE_TYPE_ICE_NONCAV, VEHICLE_TYPE_BEV, VEHICLE_TYPE_BEV_NONCAV,
VEHICLE_TYPE_PHEV, VEHICLE_TYPE_PHEV_NONCAV, VEHICLE_TYPE_HFCV, VEHICLE_TYPE_HFCV_NONCAV};

class Vehicle {
public:
	Vehicle(int id, VehicleType vehicle_type = VehicleType::UNKNOWN) {
		id_ = id;
		vehicle_type_ = vehicle_type;
		soc_ = SOC_MAX;
		battery_capacity_ = 16.5;	// for 2013 Chevrolet Volt, in the unit of kWh
		trave_time_accumulated_ = 0;
		energy1_consumed_accumulated_ = 0;
		energy2_consumed_accumulated_ = 0;
	}
	Vehicle() : Vehicle(-1) {}

	void updateEnergy(double energy1, double energy2, double sim_step) {
		// energy1 will be in the unit of gallon
		// energy2 will be in the unit of kWh
		energy1_consumed_accumulated_ += energy1;
		energy2_consumed_accumulated_ += energy2;
		trave_time_accumulated_ += sim_step;
		if (vehicle_type_ == VEHICLE_TYPE_PHEV || vehicle_type_ == VEHICLE_TYPE_PHEV_NONCAV) {
			if (!(energy2 < 0 && soc_ > SOC_MAX)) { // no more charging when soc is greater than SOC_MAX
				soc_ -= energy2 / battery_capacity_;
			}
		}
	}

	int id_;
	VehicleType vehicle_type_;
	double trave_time_accumulated_;
	double energy1_consumed_accumulated_;
	double energy2_consumed_accumulated_;	// only used for PHEV to store fuel usage
	double soc_;
	double battery_capacity_;
};

class Link {
public:
	Link(int id) {
		id_ = id;
		resetEnergy();
		resetBaseEnergy();
		resetPredictedEnergy();
	}

	Link() :Link(-1) {}

	void resetEnergy() {
		travel_time_ = 0;
		energy_ice_ = 0;
		energy_bev_ = 0;
		energy_phev1_ = 0;
		energy_phev2_ = 0;
		energy_hfcv_ = 0;
		cnt_time_ = 0;
		cnt_ice_ = 0;
		cnt_bev_ = 0;
		cnt_phev_ = 0;
		cnt_hfcv_ = 0;
	}

	void resetBaseEnergy() {
		base_travel_time_ = 0;
		base_energy_ice_ = 0;
		base_energy_bev_ = 0;
		base_energy_phev1_ = 0;
		base_energy_phev2_ = 0;
		base_energy_hfcv_ = 0;
	}


	void resetPredictedEnergy() {
		predicted_travel_time_.clear();
		predicted_energy_ice_.clear();
		predicted_energy_bev_.clear();
		predicted_energy_phev1_.clear();
		predicted_energy_phev2_.clear();
		predicted_energy_hfcv_.clear();
	}

	void addOneVehicleData(VehicleType vt, double energy1, double energy2, double travel_time) {
		double travel_time_in_s = travel_time;
		double travel_time_in_h = travel_time / 3600.0;
		switch (vt)
		{
		case VehicleType::ICE:
		case VehicleType::ICE_NONCAV:
			energy_ice_ += energy1 * travel_time_in_s;
			cnt_ice_++;
			break;
		case VehicleType::BEV:
		case VehicleType::BEV_NONCAV:
			energy_bev_ += energy2 * travel_time_in_h;
			cnt_bev_++;
			break;
		case VehicleType::PHEV:
		case VehicleType::PHEV_NONCAV:
			energy_phev1_ += energy1 * travel_time_in_s;
			energy_phev2_ += energy2 * travel_time_in_h;
			cnt_phev_++;
			break;
		case VehicleType::HFCV:
		case VehicleType::HFCV_NONCAV:
			energy_hfcv_ += energy2 * travel_time_in_h;
			cnt_phev_++;
			break;
		default:
			break;
		}

		travel_time_ += travel_time;
		cnt_time_++;
	}

	void pushBackPredictedData() {
		if (cnt_time_ == 0)
			predicted_travel_time_.push_back(base_travel_time_);
		else
			predicted_travel_time_.push_back(travel_time_ / cnt_time_);

		if (cnt_ice_ == 0)
			predicted_energy_ice_.push_back(base_energy_ice_);
		else
			predicted_energy_ice_.push_back(energy_ice_ / cnt_ice_);

		if (cnt_bev_ == 0)
			predicted_energy_bev_.push_back(base_energy_bev_);
		else
			predicted_energy_bev_.push_back(energy_bev_ / cnt_bev_);

		if (cnt_phev_ == 0) {
			predicted_energy_phev1_.push_back(base_energy_phev1_);
			predicted_energy_phev2_.push_back(base_energy_phev2_);
		}
		else {
			predicted_energy_phev1_.push_back(energy_phev1_ / cnt_phev_);
			predicted_energy_phev2_.push_back(energy_phev2_ / cnt_phev_);
		}

		if (cnt_hfcv_ == 0)
			predicted_energy_hfcv_.push_back(base_energy_hfcv_);
		else
			predicted_energy_hfcv_.push_back(energy_ice_ / cnt_hfcv_);
	}

	void printDetail() { // only used for debug
		printDebugLog("Section id is : " + to_string(id_));
		printDebugLog("Predicted ice energy is :");
		for (auto& v : predicted_energy_ice_)
			printDebugLog(to_string(v));
	}

	int id_;
	double travel_time_;
	double energy_ice_;
	double energy_bev_;
	double energy_phev1_;
	double energy_phev2_;
	double energy_hfcv_;
	double base_travel_time_;
	double base_energy_ice_;
	double base_energy_bev_;
	double base_energy_phev1_; // for the gas usage
	double base_energy_phev2_; // for the electricity usage
	double base_energy_hfcv_;
	vector<double> predicted_travel_time_;
	vector<double> predicted_energy_ice_;
	vector<double> predicted_energy_bev_;
	vector<double> predicted_energy_phev1_; // for the gas usage
	vector<double> predicted_energy_phev2_; // for the electricity usage
	vector<double> predicted_energy_hfcv_;
	int cnt_time_, cnt_ice_, cnt_bev_, cnt_phev_, cnt_hfcv_;
};

class Network {
public:
	Network() {
		map_vehicles_.clear();
		map_links_.clear();
		overall_travel_time_ = 0;
		overall_fuel_consumed_ = 0;
		overall_electricity_used_ = 0;
		N_links_ = 0;
		N_vehicles_ = 0;
		N_vehicle_types_ = 8;
		travel_time_per_vehicle_type_.clear();
		fuel_consumed_per_vehicle_type_.clear();
		electricity_used_per_vehicle_type_.clear();
		vehicle_cnt_per_vehicle_type_.clear();

		price_gas_ = PRICE_GAS;		// define gas price here
		price_electricity_ = PRICE_ELECTRICITY; // define electric price here

		history_file_replative_path_ = "data_history\\history.txt";
		summary_log_relative_path_ = "sim_logs\\summary_historical_data.csv";

		experiment_id_ = -1;
		demand_percentage_ = 100;
		prediction_horizon_ = 5;
		cav_penetration_ = 100;
		eco_routing_with_travel_time_ = false;
	}

	void saveLogs() {

		auto t = time(nullptr);
		auto tm = *localtime(&t);
		stringstream ss1;
		stringstream ss2;
		ss1 << fixed << setprecision(0) << "DTExp_Demand_" << demand_percentage_ << "_Prediction_" << prediction_horizon_ << "_CAV_" << cav_penetration_ << "_";
		ss2 << std::put_time(&tm, "%Y%m%d%H%M%S");
		string time_now = ss2.str();
		string log_file_name = PROJECT_DIR + "\\sim_logs\\" + ss1.str() + time_now + ".log.csv";
		printDebugLog(log_file_name);

		N_vehicles_ = map_vehicles_.size();

		// write to a seperate file
		ofstream fout;
		fout.open(log_file_name);
		fout << "Demand_Percentage," << demand_percentage_ << ",%\n"
			<< "Prediction_Horizon," << prediction_horizon_ << ",min\n"
			<< "CAV_Penetration," << cav_penetration_ << ",%\n"
			<< "Eco_Routing_with_Travel_Time" << eco_routing_with_travel_time_ << "\n"
			<< "Overall_Travel_Time_Avg," << overall_travel_time_ / N_vehicles_ << ",s\n"
			<< "Overall_Fuel_Used_Avg," << overall_fuel_consumed_ / N_vehicles_ << ",gallon\n"
			<< "Overall_Electricity_Used_Avg," << overall_electricity_used_ / N_vehicles_ << ",kWh\n"
			<< "Overall_Fuel_Cost_Avg," << overall_fuel_consumed_ * price_gas_ / N_vehicles_ << ",dollars\n"
			<< "Overall_Electricity_Cost_Avg," << overall_electricity_used_ * price_electricity_ / N_vehicles_ << ",dollars\n"
			<< "Total_Number_of_Vehicles," << N_vehicles_ << "\n"
			<< "Vehicle_Type,ICE,ICE_NONCAV,BEV,BEV_NONCAV,PHEV,PHEV_NONCAN,HFCV,HFCV_NONCAV\n";

		fout << "Travel_Time_VT_Avg";
		for (auto& vt : valid_vehicle_type_list) fout << "," << ((vehicle_cnt_per_vehicle_type_[vt] == 0) ? 0 : (travel_time_per_vehicle_type_[vt] / vehicle_cnt_per_vehicle_type_[vt]));
		fout << endl;
		fout << "Fuel_Used_VT_Avg";
		for (auto& vt : valid_vehicle_type_list) fout << "," << ((vehicle_cnt_per_vehicle_type_[vt] == 0) ? 0 : (fuel_consumed_per_vehicle_type_[vt] / vehicle_cnt_per_vehicle_type_[vt]));
		fout << endl;
		fout << "Electricity_Used_VT_Avg";
		for (auto& vt : valid_vehicle_type_list) fout << "," << ((vehicle_cnt_per_vehicle_type_[vt] == 0) ? 0 : (electricity_used_per_vehicle_type_[vt] / vehicle_cnt_per_vehicle_type_[vt]));
		fout << endl;
		fout << "Vehicle_Count_VT";
		for (auto& vt : valid_vehicle_type_list) fout << "," << vehicle_cnt_per_vehicle_type_[vt];
		fout << endl;
		fout.close();

		// write to the summary file
		fout.open(PROJECT_DIR + "\\" + summary_log_relative_path_, ios::out | ios::app);
		fout << time_now
			<< "," << demand_percentage_
			<< "," << prediction_horizon_
			<< "," << cav_penetration_
			<< "," << eco_routing_with_travel_time_
			<< "," << overall_travel_time_ / N_vehicles_
			<< "," << overall_fuel_consumed_ / N_vehicles_
			<< "," << overall_electricity_used_ / N_vehicles_
			<< "," << overall_fuel_consumed_ * price_gas_ / N_vehicles_
			<< "," << overall_electricity_used_ * price_electricity_
			<< "," << N_vehicles_;
		for (auto& vt : valid_vehicle_type_list) fout << "," << ((vehicle_cnt_per_vehicle_type_[vt] == 0) ? 0 : (travel_time_per_vehicle_type_[vt] / vehicle_cnt_per_vehicle_type_[vt]));
		for (auto& vt : valid_vehicle_type_list) fout << "," << ((vehicle_cnt_per_vehicle_type_[vt] == 0) ? 0 : (fuel_consumed_per_vehicle_type_[vt] / vehicle_cnt_per_vehicle_type_[vt]));
		for (auto& vt : valid_vehicle_type_list) fout << "," << ((vehicle_cnt_per_vehicle_type_[vt] == 0) ? 0 : (electricity_used_per_vehicle_type_[vt] / vehicle_cnt_per_vehicle_type_[vt]));
		for (auto& vt : valid_vehicle_type_list) fout << "," << vehicle_cnt_per_vehicle_type_[vt];
		fout << endl;
		fout.close();
	}

	void saveHistoricalData() {
		string history_file_name = PROJECT_DIR + "\\" + history_file_replative_path_;
		printDebugLog("Will write historical data to this file: " + history_file_name);

		ofstream fout;
		fout.open(history_file_name);
		fout << "#Demand_Percentage," << demand_percentage_ << ",%\n"
			<< "#Prediction_Horizon," << prediction_horizon_ << ",min\n"
			<< "#CAV_Penetration," << cav_penetration_ << ",%\n"
			<< "#TIME is in the unit of second. ICE PHEV1 are in the unit of gallon. BEV PHEV2 and HFCV are in the unit of kWh.\n"
			<< "#Section_id,cost_type,historical_data1,historical_data2,...";
		for (auto& item : map_links_) {
			auto& link = item.second;
			fout << endl << link.id_ << ",TIME";
			for (auto& value : link.predicted_travel_time_) fout << "," << value;
			fout << endl << link.id_ << ",ICE";
			for (auto& value : link.predicted_energy_ice_) fout << "," << value;
			fout << endl << link.id_ << ",BEV";
			for (auto& value : link.predicted_energy_bev_) fout << "," << value;
			fout << endl << link.id_ << ",PHEV1";
			for (auto& value : link.predicted_energy_phev1_) fout << "," << value;
			fout << endl << link.id_ << ",PHEV2";
			for (auto& value : link.predicted_energy_phev2_) fout << "," << value;
			fout << endl << link.id_ << ",HFCV";
			for (auto& value : link.predicted_energy_hfcv_) fout << "," << value;
		}
	}

	void loadHistoricalData() {

		for (auto& item : map_links_) {
			item.second.resetPredictedEnergy();
		}

		string line, word;
		int link_id;
		string cost_type;
		string history_file_name = PROJECT_DIR + "\\" + history_file_replative_path_;
		printDebugLog("Will read historical data from this file: " + history_file_name);
		ifstream fin;
		fin.open(history_file_name);

		if (!fin.is_open())
			return;

		while (getline(fin, line)) {
			if (line[0] == '#')
				continue;
			stringstream ss(line);
			getline(ss, word, ',');
			link_id = stoi(word);
			getline(ss, cost_type, ',');
			vector<double>* target_vec = nullptr;
			if (cost_type == "TIME") {
				target_vec = &map_links_[link_id].predicted_travel_time_;
			}
			else if (cost_type == "ICE") {
				target_vec = &map_links_[link_id].predicted_energy_ice_;
			}
			else if (cost_type == "BEV") {
				target_vec = &map_links_[link_id].predicted_energy_bev_;
			}
			else if (cost_type == "PHEV1") {
				target_vec = &map_links_[link_id].predicted_energy_phev1_;
			}
			else if (cost_type == "PHEV2") {
				target_vec = &map_links_[link_id].predicted_energy_phev2_;
			}
			else if (cost_type == "HFCV") {
				target_vec = &map_links_[link_id].predicted_energy_hfcv_;
			}

			while (getline(ss, word, ',')) {
				target_vec->push_back(stod(word));
			}
		}

		fin.close();
	}

	int getTotalNumberVehicles() {
		return map_vehicles_.size();
	}

	double getMaxTravelTime() {
		double the_max = 0;
		for (auto& item : map_vehicles_) {
			if (item.second.trave_time_accumulated_ > the_max)
				the_max = item.second.trave_time_accumulated_;
		}
		return the_max;
	}

	double getAverageTravelTime() {
		double avg = 0;
		for (auto& item : map_vehicles_) {
			avg += item.second.trave_time_accumulated_;
		}
		return avg / getTotalNumberVehicles();
	}

	unordered_map<int, Vehicle> map_vehicles_;
	unordered_map<int, Link> map_links_;
	double overall_travel_time_;
	double overall_fuel_consumed_;
	double overall_electricity_used_;
	unordered_map<VehicleType, double> travel_time_per_vehicle_type_;
	unordered_map<VehicleType, double> fuel_consumed_per_vehicle_type_;
	unordered_map<VehicleType, double> electricity_used_per_vehicle_type_;
	unordered_map<VehicleType, int> vehicle_cnt_per_vehicle_type_;
	int N_links_;
	int N_vehicles_;
	int N_vehicle_types_;	// only count self-defined vehicle type
	double price_gas_;		// US dollar per gallon
	double price_electricity_;	// US dollar per kWh
	string history_file_replative_path_;
	string summary_log_relative_path_;

	int experiment_id_;
	double demand_percentage_;
	double prediction_horizon_;
	double cav_penetration_;
	bool eco_routing_with_travel_time_;
};

Network network;

double eng_sum[4];

void* link_attribute_travel_time;
void* link_attribute_ice;
void* link_attribute_bev;
void* link_attribute_phev;
void* link_attribute_hfcv;

int N_link, N_turn, N_type, N_step;
ifstream flink, fturn;
double** lnk_flow, ** trn_per;
static double _null_energy_varibale = 0;


void printDebugLog(string s) {
	AKIPrintString(("## Debug log ##: " + s).c_str());
}

// calculate energy based on vehicle type, energy1 is for gas usage and energy2 is for electricity usage
void Emission(double spd, double grade, double acc, VehicleType vehicle_type, double& energy1, double& energy2, double cur_soc = 0.88)
{
	//double spd = 100.0 / 3.6; // in the unit of m/s
	//double acc = 0;           // in the unit of m/s^2
	//double grade = 0;         // without percentage
	//double cur_soc;			// current soc of the vehicle
	energy1 = 0;
	energy2 = 0;

	double spd_in_km_per_h;
	double spd_in_m_per_s;
	spd_in_m_per_s = spd;
	spd_in_km_per_h = spd * 3.6;

	double g = 9.8066, cr = 1.75, c1 = 0.0328, c2 = 4.575, rho_air = 1.2256; // fixed value
	double m = 1453, Af = 2.32, cd = 0.30, eta_d = 0.92;                     // default for 2010 Honda Accord, eta_d is the driveline efficiency
	double a0 = 0.00059217, a1 = 0.000042709, a2 = 0.000001;                 // default for 2010 Honda Accord
	double altitude = 125;                                                   // default 125m for DC
	double ch = 1 - 0.085 * (altitude / 1000);                               // default correction factor for altitude (unitless)
	double eta_dl = 0.92, eta_em = 0.91, eta_rb = 0;                         // driveline efficiency and electric motor efficiency (fixed most of the time)
	double eta_battery = 0.9;												 // for PHEV
	double soc_min = SOC_MIN, soc_max = SOC_MAX;							 // range for PHEV

	double alpha = 0.2;                                 // for FHCV, efficiency of energy produced from the battery
	double beta = 0.93;                                 // for FHCV, overall vehicle efficiency considering a driveline efficiency and the electric motor efficiency
	double va = 32;                                     // for FHCV, in the unit of km/h
	double Pa = 2.5, Pb = 5, P_aux = 0.0, P_idle = 0.0; // for FHCV, in the unit of kW

	double grade_sin, grade_cos;
	double resistance;
	double P_ice, P_w, P_batt, P_fuel;
	double P_max = 111; // max electric motor power, in the unit of kW
	double temp_total_power;

	switch (vehicle_type)
	{
	case VehicleType::ICE:
	case VehicleType::ICE_NONCAV:
		// VT-CPFM model for gasoline vehicles
		// 2010 Honda Accord, from rakha2011virginia paper
		m = 1453, cd = 0.30, Af = 2.32, eta_d = 0.92;
		a0 = 0.00059217, a1 = 0.000042709, a2 = 0.000001;
		resistance = rho_air / 25.92 * cd * ch * Af * pow(spd_in_km_per_h, 2) + g * m * cr / 1000 * (c1 * spd_in_km_per_h + c2) + g * m * grade;
		P_ice = ((resistance + 1.04 * m * acc) / (3600 * eta_d)) * spd_in_km_per_h; // in the unit of kW
		energy1 = a0;                                                             // in the unit of liter per second
		if (P_ice > 0)
			energy1 = a0 + a1 * P_ice + a2 * pow(P_ice, 2);	// in the unit of liter per second
		energy1 /= 3.78541;												// now, energy1 is in the unit of gallon per second
		//cout << "ICE P_ice is :" << P_ice << " kW" << endl;
		//cout << "ICE Fule of ice is: " << energy1 << " gallon/s" << endl;
		break;

	case VehicleType::BEV:
	case VehicleType::BEV_NONCAV:
		// CPEM model for EVs
		// 2015 Nissan Leaf, from fiori2016power paper
		m = 1521.0, Af = 2.3316, cd = 0.28;
		grade_cos = sqrt(1.0 / (1.0 + pow(grade, 2)));
		grade_sin = sqrt(1.0 - pow(grade_cos, 2));
		if (grade < 0)
			grade_sin *= -1;
		// P_w = (m * acc + m * g * grade_cos * cr / 1000.0 * (c1 * spd_in_m_per_s + c2) + 0.5 * rho_air * Af * cd * pow(spd_in_m_per_s, 2) + m * g * grade_sin) * spd_in_m_per_s; // P_w is in the unit of Watt, could be negative
		// To be consistant on the wheel power model, use the model in HFCV model
		P_w = (m * acc + m * g * cr / 1000.0 * (c1 * spd_in_km_per_h + c2) + 1.0 / 2.0 / 3.6 / 3.6 * rho_air * Af * cd * pow(spd_in_km_per_h, 2) + m * g * grade) * spd_in_km_per_h / 3.6; // P_w is in the unit of Watt, could be negative
		P_w = P_w / 1000.0;			// now P_W is in the unit of kW
		if (P_w >= 0)
		{
			energy2 = P_w / (eta_dl * eta_em); // in the unit of kW
		}
		else
		{
			if (acc < 0)
				eta_rb = 1.0 / exp(0.0411 / abs(acc));
			else
				eta_rb = 0;
			energy2 = P_w * eta_rb; // in the unit of kW
		}
		//cout << "PEV P_w is: " << P_w << " kW" << endl;
		//cout << "PEV electric motor power is: " << energy2 << " kW" << endl;
		break;

	case VehicleType::PHEV:
	case VehicleType::PHEV_NONCAV:
		// PHEV energy model
		// 2013 Chevrolet Volt, from fiori2018microscopic paper
		// this will assume the battery state-of-charge is always greater than the minimum range
		// first calculate battery power
		m = 1860, Af = 2.1851, cd = 0.29;
		a0 = 0.00035218, a1 = 0.000032141, a2 = 0.000001;
		P_max = 111;		// max electric motor power, in the unit of kW
		grade_cos = sqrt(1 / (1 + pow(grade, 2)));
		grade_sin = sqrt(1 - pow(grade_cos, 2));
		if (grade < 0)
			grade_sin *= -1;
		// P_w = (m * acc + m * g * grade_cos * cr / 1000 * (c1 * spd_in_m_per_s + c2) + 0.5 * rho_air * Af * cd * pow(spd_in_m_per_s, 2) + m * g * grade_sin) * spd_in_m_per_s; // P_w is in the unit of Watt, could be negative
		// To be consistant on the wheel power model, use the model in HFCV model
		P_w = (m * acc + m * g * cr / 1000.0 * (c1 * spd_in_km_per_h + c2) + 1.0 / 2.0 / 3.6 / 3.6 * rho_air * Af * cd * pow(spd_in_km_per_h, 2) + m * g * grade) * spd_in_km_per_h / 3.6; // P_w is in the unit of Watt, could be negative
		P_w = P_w / 1000;				// now P_W is in the unit of kW

		if (P_w >= 0)
		{
			temp_total_power = P_w / (eta_dl * eta_em); // in the unit of kW
			if (cur_soc > soc_min) {
				if (temp_total_power > P_max)
				{
					P_ice = (temp_total_power - P_max);               // in the unit of kW
					energy2 = P_max;                                // electric usage for PHEV, in the unit of kW
					energy1 = a0 + a1 * P_ice + a2 * pow(P_ice, 2); // fuel usage for PHEV, in the unit of liter per second
					energy1 /= 3.78541;                             // now, energy1 is in the unit of gallon per second
				}
				else
				{
					energy2 = temp_total_power; // electric usage for PHEV, in the unit of kW
					energy1 = 0;                // fuel usage for PHEV, in the unit of gallon per second
				}
			}
			else {
				P_ice = temp_total_power;
				energy2 = 0;                                // electric usage for PHEV, in the unit of kW
				energy1 = a0 + a1 * P_ice + a2 * pow(P_ice, 2); // fuel usage for PHEV, in the unit of liter per second
				energy1 /= 3.78541;                             // now, energy1 is in the unit of gallon per second
			}
		}
		else
		{
			if (acc < 0)
				eta_rb = 1.0 / exp(0.0411 / abs(acc));
			else
				eta_rb = 0;
			energy2 = P_w * eta_rb * eta_battery; // in the unit of kW, will be charged back to the battery
			energy1 = 0;
		}

		//cout << "PHEV P_w is: " << P_w << " kW" << endl;
		//cout << "PHEV total needed power is: " << temp_total_power << " kW" << endl;
		//cout << "PHEV needed battery power is: " << energy2 << " kW" << endl;
		//cout << "PHEV needed fuel of ICE engine is: " << energy1 << " gallon/s" << endl;
		break;

	case VehicleType::HFCV:
	case VehicleType::HFCV_NONCAV:
		// HFCV energy model
		// 2017 Toyota Mirai, from ahn2022developing paper
		m = 1928.0, Af = 2.3316, cd = 0.28;
		grade_cos = sqrt(1.0 / (1.0 + pow(grade, 2)));
		grade_sin = sqrt(1.0 - pow(grade_cos, 2));
		if (grade < 0)
			grade_sin *= -1;
		// P_w = (m * acc + m * g * grade_cos * cr / 1000.0 * (c1 * spd_in_m_per_s + c2) + 0.5 * rho_air * Af * cd * pow(spd_in_m_per_s, 2) + m * g * grade_sin) * spd_in_m_per_s; // P_w is in the unit of Watt, could be negative
		// To be consistant on the wheel power model, use the model in HFCV model
		P_w = (m * acc + m * g * cr / 1000.0 * (c1 * spd_in_km_per_h + c2) + 1.0 / 2.0 / 3.6 / 3.6 * rho_air * Af * cd * pow(spd_in_km_per_h, 2) + m * g * grade) * spd_in_km_per_h / 3.6; // P_w is in the unit of Watt, could be negative
		P_w = P_w / 1000; // now P_W is in the unit of kW
		if (P_w >= 0) {
			if (P_w <= Pa)
			{
				P_batt = P_w + P_aux;
				P_fuel = P_idle;
			}
			else
			{
				if (P_w > Pa && P_w <= Pb)
				{
					if (spd_in_km_per_h <= va)
					{
						P_batt = P_w + P_aux;
						P_fuel = P_idle;
					}
					else
					{
						P_batt = P_idle + P_aux;
						P_fuel = P_w * beta;
					}
				}
				else
				{
					P_batt = P_w * alpha + P_aux;
					P_fuel = P_w * beta;
				}
			}
			energy2 = P_batt + P_fuel;    // in the unit of kW
		}
		else
		{
			if (acc < 0)
				eta_rb = 1.0 / exp(0.0411 / abs(acc));
			else
				eta_rb = 0;
			energy2 = P_w * eta_rb; // in the unit of kW
		}
		//cout << "HFCV P_w is: " << P_w << " kW" << endl;
		//cout << "HFCV need battery power: " << P_batt << " kW" << endl;
		//cout << "HFCV fuel cell power is : " << P_fuel << " kW" << endl;
		//cout << "HFCV total power is : " << energy2 << " kW" << endl;
		break;

	default:
		break;
	}
}

// this function will return energy cost for all vehicle type, E_phev1 if for gas usage, E_phev2 is for electricity usage
// only used to calculate base energy cost
void Emission(double spd, double grade, double acc, double& E_ice, double& E_bev, double& E_phev1, double& E_phev2, double& E_hfcv)
{
	double dummy_energy_variable;
	Emission(spd, grade, acc, VehicleType::ICE, E_ice, dummy_energy_variable);
	Emission(spd, grade, acc, VehicleType::BEV, dummy_energy_variable, E_bev);
	Emission(spd, grade, acc, VehicleType::PHEV, E_phev1, E_phev2);
	Emission(spd, grade, acc, VehicleType::HFCV, dummy_energy_variable, E_hfcv);
}

int AAPILoad()
{
	// srand((uint32_t)time(NULL));
	printDebugLog("In AAPILoad()");
	return 0;
}

int AAPIInit()
{
	//printDebugLog("In AAPIInit()");
	ANGConnEnableVehiclesInBatch(true);

	// define new attributes for link travel time cost and fuel cost (fuel costs need to be transfered into money cost before updating into the attributes)
	// Aimsun will calculate shortest path based on the values in the attributes
	link_attribute_ice = ANGConnGetAttribute(AKIConvertFromAsciiString("GKSection::ice"));
	if (link_attribute_ice == NULL) {
		link_attribute_ice = ANGConnCreateAttribute(AKIConvertFromAsciiString("GKSection"), AKIConvertFromAsciiString("GKSection::ice"), AKIConvertFromAsciiString("ice"), DOUBLE_TYPE, EXTERNAL);
	}
	link_attribute_bev = ANGConnGetAttribute(AKIConvertFromAsciiString("GKSection::bev"));
	if (link_attribute_bev == NULL) {
		link_attribute_bev = ANGConnCreateAttribute(AKIConvertFromAsciiString("GKSection"), AKIConvertFromAsciiString("GKSection::bev"), AKIConvertFromAsciiString("bev"), DOUBLE_TYPE, EXTERNAL);
	}
	link_attribute_phev = ANGConnGetAttribute(AKIConvertFromAsciiString("GKSection::phev"));
	if (link_attribute_phev == NULL) {
		link_attribute_phev = ANGConnCreateAttribute(AKIConvertFromAsciiString("GKSection"), AKIConvertFromAsciiString("GKSection::phev"), AKIConvertFromAsciiString("phev"), DOUBLE_TYPE, EXTERNAL);
	}
	link_attribute_hfcv = ANGConnGetAttribute(AKIConvertFromAsciiString("GKSection::hfcv"));
	if (link_attribute_hfcv == NULL) {
		link_attribute_hfcv = ANGConnCreateAttribute(AKIConvertFromAsciiString("GKSection"), AKIConvertFromAsciiString("GKSection::hfcv"), AKIConvertFromAsciiString("hfcv"), DOUBLE_TYPE, EXTERNAL);
	}
	link_attribute_travel_time = ANGConnGetAttribute(AKIConvertFromAsciiString("GKSection::travel_time"));
	if (link_attribute_travel_time == NULL) {
		link_attribute_travel_time = ANGConnCreateAttribute(AKIConvertFromAsciiString("GKSection"), AKIConvertFromAsciiString("GKSection::travel_time"), AKIConvertFromAsciiString("travel_time"), DOUBLE_TYPE, EXTERNAL);
	}

	// build the local link map in network
	// and update the initial link energy consumption
	void* attribute_demand = ANGConnGetAttribute(AKIConvertFromAsciiString("GKExperiment::demand_percentage"));
	void* attribute_prediction = ANGConnGetAttribute(AKIConvertFromAsciiString("GKExperiment::prediction_horizon"));
	void* attribute_cav_penetration = ANGConnGetAttribute(AKIConvertFromAsciiString("GKExperiment::cav_penetration"));
	void* attribute_eco_routing_with_travel_time = ANGConnGetAttribute(AKIConvertFromAsciiString("GKExperiment::eco_routing_with_travel_time"));
	network.experiment_id_ = ANGConnGetExperimentId();
	network.demand_percentage_ = ANGConnGetAttributeValueDouble(attribute_demand, network.experiment_id_);
	network.prediction_horizon_ = ANGConnGetAttributeValueDouble(attribute_prediction, network.experiment_id_);
	network.cav_penetration_ = ANGConnGetAttributeValueDouble(attribute_cav_penetration, network.experiment_id_);
	network.eco_routing_with_travel_time_ = ANGConnGetAttributeValueBool(attribute_eco_routing_with_travel_time, network.experiment_id_);
	network.N_links_ = AKIInfNetNbSectionsANG(); // obtain the number of links in the network
	printDebugLog("Total number of links is : " + to_string(network.N_links_));
	for (int i = 0; i < network.N_links_; i++) {
		int secid = AKIInfNetGetSectionANGId(i);

		if (network.map_links_.count(secid) < 1)
			network.map_links_[secid] = Link(secid);

		A2KSectionInf& secinf = AKIInfNetGetSectionANGInf(secid);
		double spd = secinf.speedLimit / 3.6;
		double grade = secinf.slopePercentages[0] / 100.0;

		Emission(spd, grade, 0.0, network.map_links_[secid].base_energy_ice_, network.map_links_[secid].base_energy_bev_, network.map_links_[secid].base_energy_phev1_,
			network.map_links_[secid].base_energy_phev2_, network.map_links_[secid].base_energy_hfcv_);

		double _free_flow_travel_time_in_s = secinf.length / spd;
		double _free_flow_travel_time_in_h = _free_flow_travel_time_in_s / 3600.0;
		// record base cost
		network.map_links_[secid].base_travel_time_ = _free_flow_travel_time_in_s;
		network.map_links_[secid].base_energy_ice_ *= _free_flow_travel_time_in_s;
		network.map_links_[secid].base_energy_bev_ *= _free_flow_travel_time_in_h;
		network.map_links_[secid].base_energy_phev1_ *= _free_flow_travel_time_in_s;
		network.map_links_[secid].base_energy_phev2_ *= _free_flow_travel_time_in_h;
		network.map_links_[secid].base_energy_hfcv_ *= _free_flow_travel_time_in_h;

		//ANGConnSetAttributeValueDouble(link_attribute_ice, secid, network.map_links_[secid].base_energy_ice_ * network.price_gas_);
		//ANGConnSetAttributeValueDouble(link_attribute_bev, secid, network.map_links_[secid].base_energy_bev_ * network.price_electricity_);
		//ANGConnSetAttributeValueDouble(link_attribute_phev, secid, network.map_links_[secid].base_energy_phev1_ * network.price_gas_ + network.map_links_[secid].base_energy_phev2_ * network.price_electricity_);
		//ANGConnSetAttributeValueDouble(link_attribute_hfcv, secid, network.map_links_[secid].base_energy_hfcv_ * network.price_electricity_);

		// init cost will be the travel time of the link
		ANGConnSetAttributeValueDouble(link_attribute_ice, secid, network.map_links_[secid].base_travel_time_);
		ANGConnSetAttributeValueDouble(link_attribute_bev, secid, network.map_links_[secid].base_travel_time_);
		ANGConnSetAttributeValueDouble(link_attribute_phev, secid, network.map_links_[secid].base_travel_time_);
		ANGConnSetAttributeValueDouble(link_attribute_hfcv, secid, network.map_links_[secid].base_travel_time_);
		ANGConnSetAttributeValueDouble(link_attribute_travel_time, secid, network.map_links_[secid].base_travel_time_);
	}

	if (!GENERATE_HISTORICAL_DATA)
		network.loadHistoricalData();

	return 0;
}


int AAPIManage(double time, double timeSta, double timTrans, double acicle)
{
	double sim_step = AKIGetSimulationStepTime();

	// update the total energy consumption for different types of vehicles in every simlulation step
	if (time - timTrans > 0) {
		for (auto& item : network.map_links_) {
			int secid = item.first;
			A2KSectionInf& secinf = AKIInfNetGetSectionANGInf(secid);
			double grade = secinf.slopePercentages[0] / 100.0;
			int nbveh = AKIVehStateGetNbVehiclesSection(secid, true);
			for (int k = 0; k < nbveh; k++) {
				InfVeh& vehinf = AKIVehStateGetVehicleInfSection(secid, k);
				VehicleType vehicle_type = static_cast<VehicleType>(AKIVehTypeGetIdVehTypeANG(vehinf.type));
				double spd = vehinf.CurrentSpeed / 3.6;												// speed in m/s
				double acc = (vehinf.CurrentSpeed - vehinf.PreviousSpeed) / (3.6 * sim_step);		// acceleration in m/s^2

				if (network.map_vehicles_.count(vehinf.idVeh) < 1) {
					network.map_vehicles_[vehinf.idVeh] = Vehicle(vehinf.idVeh, vehicle_type);
					network.vehicle_cnt_per_vehicle_type_[vehicle_type]++;
				}

				Vehicle& cur_vehicle = network.map_vehicles_[vehinf.idVeh];

				double energy1; // fuel usage
				double energy2; // electricity usage
				Emission(spd, grade, acc, vehicle_type, energy1, energy2, cur_vehicle.soc_);
				energy1 *= sim_step; // now will be in the unit of Gallon
				energy2 *= sim_step / 3600.0; // now will be in the unit of kWh
				cur_vehicle.updateEnergy(energy1, energy2, sim_step);

				network.overall_fuel_consumed_ += energy1;
				network.fuel_consumed_per_vehicle_type_[vehicle_type] += energy1;
				network.overall_electricity_used_ += energy2;
				network.electricity_used_per_vehicle_type_[vehicle_type] += energy2;
				network.overall_travel_time_ += sim_step;
				network.travel_time_per_vehicle_type_[vehicle_type] += sim_step;
			}
		}
	}

	if (GENERATE_HISTORICAL_DATA) {
		// only used to record historical data, recorded data will be write to file if simulation ends successfully
		// record cost for every link every 5 minutes
		static double next_record_time = 0; // seconds
		double record_interval = HISTORICAL_DATA_INTERVAL; // seconds
		double sim_time = time - timTrans;
		if (sim_time >= next_record_time - sim_step) {
			next_record_time += record_interval;
			for (auto& item : network.map_links_) {
				int secid = item.first;
				Link& current_link = item.second;
				A2KSectionInf& secinf = AKIInfNetGetSectionANGInf(secid);
				double grade = secinf.slopePercentages[0] / 100.0;
				int nbveh = AKIVehStateGetNbVehiclesSection(secid, true);
				current_link.resetEnergy();
				for (int k = 0; k < nbveh; k++) {
					InfVeh& vehinf = AKIVehStateGetVehicleInfSection(secid, k);
					VehicleType vehicle_type = static_cast<VehicleType>(AKIVehTypeGetIdVehTypeANG(vehinf.type));
					double spd = vehinf.CurrentSpeed / 3.6;												// speed in m/s
					double acc = (vehinf.CurrentSpeed - vehinf.PreviousSpeed) / (3.6 * sim_step);		// acceleration in m/s^2
					double energy1; // fuel usage
					double energy2; // electricity usage
					Emission(spd, grade, acc, vehicle_type, energy1, energy2);
					double section_travel_time = secinf.length / max(spd, 2.0);
					current_link.addOneVehicleData(vehicle_type, energy1, energy2, section_travel_time);
				}
				current_link.pushBackPredictedData();
			}
		}
	}

	return 0;
}

int AAPIPostManage(double time, double timeSta, double timTrans, double acicle)
{
	return 0;
}

int AAPIFinish()
{
	printDebugLog("Demand_Percentage is " + to_string(network.demand_percentage_) + "% , Prediction_Horizon is " + to_string(network.prediction_horizon_) + " min.");
	printDebugLog("CAV_Penetration is " + to_string(network.cav_penetration_) + "% , Eco_Routing_with_Travel_Time is " + to_string(network.eco_routing_with_travel_time_) + " .");
	printDebugLog("Network overall fuel usage is : " + to_string(network.overall_fuel_consumed_) + " Gallon.");
	printDebugLog("Network overall electricity usage is : " + to_string(network.overall_electricity_used_) + " kWh.");
	printDebugLog("Network overall travel time is : " + to_string(network.overall_travel_time_) + " seconds.");
	printDebugLog("Total simulated vehicles: " + to_string(network.getTotalNumberVehicles()));
	printDebugLog("Average travel time is : " + to_string(network.getAverageTravelTime()) + " seconds. Maimum of travel time is " + to_string(network.getMaxTravelTime()) + " seconds.");

	if (AKIGetCurrentSimulationTime() < AKIGetEndSimTime()) {
		printDebugLog("Simulation is not finished, no output files.");
		return 0;
	}

	printDebugLog("Statistic will be output to the logs directory");
	network.saveLogs();			// save simulation statistics in the \logs dir
	if (GENERATE_HISTORICAL_DATA)
		network.saveHistoricalData();	// save historical data to the file

	return 0;
}

int AAPIUnLoad()
{
	//AKIPrintString("UNLOAD");
	return 0;
}

int AAPIPreRouteChoiceCalculation(double time, double timeSta)
{
	// update link cost before rerouting
	double predition_horizon_in_s = network.prediction_horizon_ * 60;
	double timTrans = AKIGetDurationTransTime();
	if (time < timTrans) // do not update cost in the warm up period
		return 0;

	int predict_cost_idx = (int)round((time - timTrans + predition_horizon_in_s) / HISTORICAL_DATA_INTERVAL);
	if (!GENERATE_HISTORICAL_DATA) {
		// update link cost based on the predicted data
		for (auto& item : network.map_links_) {
			int secid = item.first;
			Link& current_link = item.second;

			A2KSectionInf& secinf = AKIInfNetGetSectionANGInf(secid);
			double spd = secinf.speedLimit / 3.6;
			double grade = secinf.slopePercentages[0] / 100;
			double _free_flow_travel_time_in_s = secinf.length / spd;
			double _free_flow_travel_time_in_h = _free_flow_travel_time_in_s / 3600.0;

			// this is the cost for non-cavs and will use base travel time
			ANGConnSetAttributeValueDouble(link_attribute_travel_time, secid, current_link.base_travel_time_);

			// these are the cost for all cavs
			if (!network.eco_routing_with_travel_time_) {
				ANGConnSetAttributeValueDouble(link_attribute_ice, secid, current_link.predicted_energy_ice_[predict_cost_idx] * network.price_gas_);
				ANGConnSetAttributeValueDouble(link_attribute_bev, secid, current_link.predicted_energy_bev_[predict_cost_idx] * network.price_electricity_);
				ANGConnSetAttributeValueDouble(link_attribute_phev, secid, current_link.predicted_energy_phev1_[predict_cost_idx] * network.price_gas_ + current_link.predicted_energy_phev2_[predict_cost_idx] * network.price_electricity_);
				ANGConnSetAttributeValueDouble(link_attribute_hfcv, secid, current_link.predicted_energy_hfcv_[predict_cost_idx] * network.price_electricity_);
			}
			else {
				ANGConnSetAttributeValueDouble(link_attribute_ice, secid, current_link.predicted_travel_time_[predict_cost_idx]);
				ANGConnSetAttributeValueDouble(link_attribute_bev, secid, current_link.predicted_travel_time_[predict_cost_idx]);
				ANGConnSetAttributeValueDouble(link_attribute_phev, secid, current_link.predicted_travel_time_[predict_cost_idx]);
				ANGConnSetAttributeValueDouble(link_attribute_hfcv, secid, current_link.predicted_travel_time_[predict_cost_idx]);
			}
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