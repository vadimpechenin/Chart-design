"""Программа-аналог для построения гистограммы для диссертации Balance_virtual_KND_real_sbor"""
from constant_file import Constant_chart
from mat_data_for_chart import Mat_load_data
import Module_basic_functions as Mbf
import numpy as np
# Инициализируем константы для расчетов
ai_const = Constant_chart()
# Загружаем данные из mat файла по результатам расчетов сборочных положений
string_path_file='D:\\Задачи в MATLAB\\ИИР__2018_2019\\2019\\КНД_результаты\\Results_KND_7details_26072019_Optimize_3scheme.mat'
ai_data=Mat_load_data(string_path_file)
ai_data.import_mat_io(string_path_file)
# Приведение в соответствие систем координат, расчет дисбалансов
Table_l_dis, Center, Disbal,STD_centers,Table_l_dis_dist_bal,Center_dist_bal,Disbal_dist_bal,STD_centers_dist_bal=Mbf.function_comput_dibal(ai_data.Table_all,ai_data.RT_array,ai_data.T_array,ai_data.korteg,ai_const.Center_ost,ai_const.delt_Z,ai_const.M,ai_const.Z)
# Некоторые корректировки
STD_centers=np.zeros((ai_data.Table_all.shape[1],1))
STD_centers=Disbal_dist_bal.T
Disbal=np.zeros((ai_data.Table_all.shape[1],1))
Disbal=Disbal_dist_bal.T*1000
korteg_new=np.zeros((ai_data.korteg.shape[0],ai_data.korteg.shape[1]))
korteg_new=ai_data.korteg*180/np.pi
#Функция вычисления таблиц для построения графиков
Table1_percent, Table2_percent,Table3_percent,Table4_percent=Mbf.norm_Table_KND(ai_data.Table_all,korteg_new,Disbal,ai_const.k,ai_const.T_dopusk,ai_const.T_dopusk_Db)
#Первоначальный график
Mbf.chart_plot(Table1_percent, Table2_percent,Table3_percent,Table4_percent)
Mbf.chart_bar(Table1_percent, Table2_percent,Table3_percent,Table4_percent)
g=0