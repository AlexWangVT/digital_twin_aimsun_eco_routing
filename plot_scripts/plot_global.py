import pandas as pd

df1 = pd.read_csv('../sim_logs/backups/batch_summary/20230809/summary_demand100.csv')
df2 = pd.read_csv('../sim_logs/backups/batch_summary/20230809/summary_demand150.csv')
df3 = pd.read_csv('../sim_logs/backups/batch_summary/20230809/summary_demand50.csv')

df = pd.concat([df1,df2,df3], ignore_index=True)
# print(df.shape)

statistics_df = df.copy(deep=True)
statistics_df['Overall_Fuel_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['Overall_Electricity_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['ICE_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['BEV_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['PHEV_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['HFCV_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['Money_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['ICE_TravelTime_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['BEV_TravelTime_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['PHEV_TravelTime_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
statistics_df['HFCV_TravelTime_Reduction'] = statistics_df['Fuel_Used_ICE_Avg(gallon)'].copy()
for idx in statistics_df.index:
    base_df = statistics_df[(statistics_df['Vehicle_Fleet'] == statistics_df['Vehicle_Fleet'][idx]) 
             & (statistics_df['Demand_Percentage(%)'] == statistics_df['Demand_Percentage(%)'][idx]) 
             & (statistics_df['Eco_Routing_with_Travel_Time(0/1)'] == statistics_df['Eco_Routing_with_Travel_Time(0/1)'][idx]) 
             & (statistics_df['Prediction_Horizon(min)'] == statistics_df['Prediction_Horizon(min)'][idx]) 
             & (statistics_df['CAV_Penetration(%)'] == 0)]
    base_overall_fuel = base_df['Overall_Fuel_Used(gallon)'].values[0]
    base_overall_electricity = base_df['Overall_Electricity_Used(kWh)'].values[0]
    base_ice = base_df['Fuel_Used_ICE_NONCAV_Avg(gallon)'].values[0]
    base_bev = base_df['Electricity_Used_BEV_NONCAV_Avg(kWh)'].values[0]
    base_phev = base_df['Electricity_Used_PHEV_NONCAV_Avg(kWh)'].values[0]
    base_hfcv = base_df['Electricity_Used_HFCV_NONCAV_Avg(kWh)'].values[0]
    base_money = base_df['Overall_Fuel_Cost($)'].values[0] + base_df['Overall_Electricity_Cost($)'].values[0] 
    
    base_traveltime_ice = base_df['Travel_Time_ICE_NONCAV_Avg(s)'].values[0]
    base_traveltime_bev = base_df['Travel_Time_BEV_NONCAV_Avg(s)'].values[0]
    base_traveltime_phev = base_df['Travel_Time_PHEV_NONCAV_Avg(s)'].values[0]
    base_traveltime_hfcv = base_df['Travel_Time_HFCV_NONCAV_Avg(s)'].values[0]
    
    statistics_df['Overall_Fuel_Reduction'][idx] = (statistics_df['Overall_Fuel_Used(gallon)'][idx] - base_overall_fuel) / base_overall_fuel * 100.0
    statistics_df['Overall_Electricity_Reduction'][idx] = (statistics_df['Overall_Electricity_Used(kWh)'][idx] - base_overall_electricity) / base_overall_electricity * 100.0
    if statistics_df['Fuel_Used_ICE_Avg(gallon)'][idx] > 0:
        statistics_df['ICE_Reduction'][idx] = (statistics_df['Fuel_Used_ICE_Avg(gallon)'][idx] - base_ice) / base_ice * 100.0
    else:
        statistics_df['ICE_Reduction'][idx] = 0.0
        
    if statistics_df['Electricity_Used_BEV_Avg(kWh)'][idx] > 0:
        statistics_df['BEV_Reduction'][idx] = (statistics_df['Electricity_Used_BEV_Avg(kWh)'][idx] - base_bev) / base_bev * 100.0
    else:
        statistics_df['BEV_Reduction'][idx] = 0.0
        
    if statistics_df['Electricity_Used_PHEV_Avg(kWh)'][idx] > 0:
        statistics_df['PHEV_Reduction'][idx] = (statistics_df['Electricity_Used_PHEV_Avg(kWh)'][idx] - base_phev) / base_phev * 100.0
    else:
        statistics_df['PHEV_Reduction'][idx] = 0.0
        
    if statistics_df['Electricity_Used_HFCV_Avg(kWh)'][idx] > 0:
        statistics_df['HFCV_Reduction'][idx] = (statistics_df['Electricity_Used_HFCV_Avg(kWh)'][idx] - base_hfcv) / base_hfcv * 100.0
    else:
        statistics_df['HFCV_Reduction'][idx] = 0.0
    statistics_df['Money_Reduction'][idx] = (statistics_df['Overall_Fuel_Cost($)'][idx] + statistics_df['Overall_Electricity_Cost($)'][idx] - base_money ) / base_money * 100.0
    
    if statistics_df['Travel_Time_ICE_Avg(s)'][idx] > 0:
        statistics_df['ICE_TravelTime_Reduction'][idx] = (statistics_df['Travel_Time_ICE_Avg(s)'][idx] - base_traveltime_ice) / base_traveltime_ice * 100.0
    else:
        statistics_df['ICE_TravelTime_Reduction'][idx] = 0.0
    if statistics_df['Travel_Time_BEV_Avg(s)'][idx] > 0:
        statistics_df['BEV_TravelTime_Reduction'][idx] = (statistics_df['Travel_Time_BEV_Avg(s)'][idx] - base_traveltime_bev) / base_traveltime_bev * 100.0
    else:
        statistics_df['BEV_TravelTime_Reduction'][idx] = 0.0
    if statistics_df['Travel_Time_PHEV_Avg(s)'][idx] > 0:
        statistics_df['PHEV_TravelTime_Reduction'][idx] = (statistics_df['Travel_Time_PHEV_Avg(s)'][idx] - base_traveltime_phev) / base_traveltime_phev * 100.0
    else:
        statistics_df['PHEV_TravelTime_Reduction'][idx] = 0.0
    if statistics_df['Travel_Time_HFCV_Avg(s)'][idx] > 0:
        statistics_df['HFCV_TravelTime_Reduction'][idx] = (statistics_df['Travel_Time_HFCV_Avg(s)'][idx] - base_traveltime_hfcv) / base_traveltime_hfcv * 100.0
    else:
        statistics_df['HFCV_TravelTime_Reduction'][idx] = 0.0
