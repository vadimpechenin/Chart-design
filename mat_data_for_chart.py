"""Создание класса для загрузки данных из мат.файла"""
import numpy as np
#Библеотека импорта mat файла
import scipy.io

class Mat_load_data():
    """Класс для импорта данных в сhart"""
    def __init__(self,string_path_file):
        self.Table_all=[]
        self.korteg=[]
        self.RT_array=[]
        self.T_array=[]
        self.dist_array=[]

    def import_mat_io(self,string_path_file):
        mat = scipy.io.loadmat(string_path_file)
        self.Table_all=np.array(mat['Table_all'])
        self.korteg = np.array(mat['korteg'])
        self.RT_array = np.array(mat['RT_array'])
        self.T_array = np.array(mat['T_array'])
        self.dist_array = np.array(mat['dist_array'])