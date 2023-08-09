import plot_global
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os


def plotAllTypes(vehicle_fleet='2030',
                 demand_percentage=100,
                 eco_routing_with_travel_time=0,
                 prediction_horizon=5,
                 save_fig=False):

    linewidths = [2, 2.5, 2.8, 2]
    linestyles = ['solid', 'dashed', 'dotted', 'dashdot']

    sns.set_theme(style='whitegrid', palette='bright')
    # sns.set_style("whitegrid")
    statistics_df = plot_global.statistics_df
    target_df = statistics_df[(statistics_df['Vehicle_Fleet'] == vehicle_fleet)
                              & (statistics_df['Demand_Percentage(%)'] == demand_percentage)
                              & (statistics_df['Eco_Routing_with_Travel_Time(0/1)'] == eco_routing_with_travel_time)
                              & (statistics_df['Prediction_Horizon(min)'] == prediction_horizon)]

    g = sns.lineplot(x='CAV_Penetration(%)', y='ICE_Reduction',
                     data=target_df, label='ICE_Reduction', linewidth=linewidths[0], linestyle=linestyles[0])
    if vehicle_fleet != 'ICE':
        g = sns.lineplot(x='CAV_Penetration(%)', y='BEV_Reduction',
                        data=target_df, label='BEV_Reduction', linewidth=linewidths[1], linestyle=linestyles[1])
        g = sns.lineplot(x='CAV_Penetration(%)', y='PHEV_Reduction',
                        data=target_df, label='PHEV_Reduction', linewidth=linewidths[2], linestyle=linestyles[2])
        g = sns.lineplot(x='CAV_Penetration(%)', y='HFCV_Reduction',
                        data=target_df, label='HFCV_Reduction', linewidth=linewidths[3], linestyle=linestyles[3])
    if prediction_horizon == 0:
        title = '{}: Demand {}% No Prediction '.format(
            vehicle_fleet, demand_percentage)
    else:
        title = '{}: Demand {}% Predict {}min '.format(
            vehicle_fleet, demand_percentage, prediction_horizon)
    if eco_routing_with_travel_time:
        title += 'Travel time routing'
    else:
        title += 'Eco-routing'
    plt.title(title, fontweight='bold')
    plt.ylabel("Energy Reduction (%)")
    plt.xlabel("CAV Penetration (%)")
    plt.ylim([-40, 40])
    plt.xlim([0, 100])
    plt.legend()
    if save_fig:
        file_name = 'figures/{}_Demand_{}_Predict_{}'.format(
            vehicle_fleet, demand_percentage, prediction_horizon)
        if eco_routing_with_travel_time:
            file_name += '_TravelTimeRouting.pdf'
        else:
            file_name += '_EnergyEcorouting.pdf'
        if os.path.isfile(file_name):
            os.remove(file_name)
        plt.savefig(file_name)  # , format='pdf')

    plt.show(block=False)


def plotComparisonEcoRouting(vehicle_type='ICE',
                             vehicle_fleet='2030',
                             demand_percentage=100,
                             prediction_horizon=5,
                             save_fig=False):
    eco_routing_with_travel_time = 0
    linewidths = [2, 2.5, 2.8, 2]
    linestyles = ['solid', 'dashed', 'dotted', 'dashdot']

    sns.set_theme(style='whitegrid', palette='bright')
    if vehicle_type == 'ICE':
        vehicle_type_color = sns.color_palette("bright")[0]
    elif vehicle_type == 'BEV':
        vehicle_type_color = sns.color_palette("bright")[1]
    elif vehicle_type == 'PHEV':
        vehicle_type_color = sns.color_palette("bright")[2]
    elif vehicle_type == 'HFCV':
        vehicle_type_color = sns.color_palette("bright")[3]

    statistics_df = plot_global.statistics_df
    target_df_eco_routing = statistics_df[(statistics_df['Vehicle_Fleet'] == vehicle_fleet)
                                          & (statistics_df['Demand_Percentage(%)'] == demand_percentage)
                                          & (statistics_df['Eco_Routing_with_Travel_Time(0/1)'] == 0)
                                          & (statistics_df['Prediction_Horizon(min)'] == prediction_horizon)]
    target_df_travel_time_routing = statistics_df[(statistics_df['Vehicle_Fleet'] == vehicle_fleet)
                                                  & (statistics_df['Demand_Percentage(%)'] == demand_percentage)
                                                  & (statistics_df['Eco_Routing_with_Travel_Time(0/1)'] == 1)
                                                  & (statistics_df['Prediction_Horizon(min)'] == prediction_horizon)]

    g = sns.lineplot(x='CAV_Penetration(%)', y=vehicle_type+'_Reduction',
                     data=target_df_eco_routing, label=vehicle_type+'_Energy Eco-routing', linewidth=2, linestyle='solid', color=vehicle_type_color)
    g = sns.lineplot(x='CAV_Penetration(%)', y=vehicle_type+'_Reduction',
                     data=target_df_travel_time_routing, label=vehicle_type+'_Travel time routing', linewidth=2.5, linestyle='dashed', color=vehicle_type_color)

    if prediction_horizon == 0:
        title = '{}: Demand {}% No Prediction '.format(
            vehicle_fleet, demand_percentage)
    else:
        title = '{}: Demand {}% Predict {}min '.format(
            vehicle_fleet, demand_percentage, prediction_horizon)
    title += vehicle_type
    plt.title(title, fontweight='bold')
    plt.ylabel("Energy Reduction (%)")
    plt.xlabel("CAV Penetration (%)")
    plt.ylim([-40, 40])
    plt.xlim([0, 100])
    # plt.legend()
    if save_fig:
        file_name = 'figures/CompareEcorouting_{}_{}_Demand_{}_Predict_{}.pdf'.format(
            vehicle_type, vehicle_fleet, demand_percentage, prediction_horizon)
        if os.path.isfile(file_name):
            os.remove(file_name)
        plt.savefig(file_name)  # , format='pdf')

    plt.show(block=False)


def plotOptimumCAV(
    demand_percentage = 100,
    save_fig = False):

    sns.set_theme(style='whitegrid', palette='bright')
    statistics_df = plot_global.statistics_df
    target_df = statistics_df[(statistics_df['Demand_Percentage(%)'] == demand_percentage)]
                            # & (statistics_df['Eco_Routing_with_Travel_Time(0/1)'] == 0)]
    g = sns.lineplot(data=target_df, x="CAV_Penetration(%)", y="Money_Reduction", color = 'red')
    g = sns.scatterplot(data=target_df, x="CAV_Penetration(%)", y="Money_Reduction", hue = 'Prediction_Horizon(min)')

    title = 'Demand {}%'.format(demand_percentage)
    plt.title(title, fontweight='bold')
    plt.ylabel("Total Reduction (%)")
    plt.xlabel("CAV Penetration (%)")
    plt.ylim([-30, 40])
    plt.xlim([0, 100])
    # plt.legend()
    if save_fig:
        file_name = 'figures/Optimum_CAV_Penetration_Demand_{}.pdf'.format(
            demand_percentage)
        if os.path.isfile(file_name):
            os.remove(file_name)
        plt.savefig(file_name)  # , format='pdf')

    plt.show(block=False)

if __name__ == '__main__':
    vehicle_fleet = '2030'
    demand_percentage = 100
    eco_routing_with_travel_time = 0
    prediction_horizon = 5
    plotAllTypes(vehicle_fleet, demand_percentage,
                 eco_routing_with_travel_time, prediction_horizon, save_fig=True)
