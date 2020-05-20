""" Модуль базовых функций для работы в chart"""
import numpy as np
import math

def function_of_center_with_support_Z(nm_line,O,O_opora_2,M):
    """Функция нахождения результирующих дисбалансов сборки"""
    # 2. Разворот
    # Угол вокруг оси х
    angle_x = math.pi / 2 - math.atan(nm_line[3-1] / nm_line[2-1])
    if (nm_line[2-1] == 0):
        angle_x = 0
    M_x = np.array([[1, 0, 0],[0, math.cos(angle_x), math.sin(angle_x)],[0, -math.sin(angle_x), math.cos(angle_x)]])
    nm_n = nm_line.dot(M_x)
    O_prov = O_opora_2.dot(M_x)
    # Разворот вокруг оси y
    angle_y = math.pi / 2 - math.atan(nm_n[3-1] / nm_n[1-1])
    if (nm_n[3-1] == 0):
        angle_y = 0
    M_y = np.array([[np.cos(angle_y), 0, np.sin(angle_y)],[0, 1, 0],[-np.sin(angle_y), 0, np.cos(angle_y)]])
    nm_n = nm_n.dot(M_y)
    O_prov = O_prov.dot(M_y)
    if (O_prov[2] > 0):
        an_dop = np.pi
    else:
        an_dop = 0
    M_y_ = np.array([[np.cos(an_dop), 0, np.sin(an_dop)], [0, 1, 0], [-np.sin(an_dop), 0, np.cos(an_dop)]])
    nm_n=nm_n.dot(M_y_)
    # Повороты
    O_ = O.dot(M_x).dot(M_y).dot(M_y_)
    #%1.8 Отдельные дисбалансы и суммарный дисбаланс
    #%Общий центр масс
    #idx=int(np.linspace(0, 6,num=7, endpoint=True))
    idx=np.array(range(7))
    X_sum=sum(O_[idx,1-1].T*M[idx]/math.fsum(M[idx]))
    Y_sum=sum(O_[idx,2-1].T*M[idx]/math.fsum(M[idx]))
    Z_sum=sum(O_[idx,3-1].T*M[idx]/math.fsum(M[idx]))
    Center=np.array([X_sum, Y_sum, Z_sum])
    # %Дисбалансы и расстояния (вместе с суммарным)
    Table_l_dis=np.zeros((2,8))
    Table_l_dis[0,:-1], Table_l_dis[0,-1],\
    Table_l_dis[1,:-1], Table_l_dis[1,-1]= ((O_[idx, 1-1]**2 + O_[idx, 2-1]**2)**0.5).T,\
                                           (Center[1-1]**2+Center[2-1]**2)**(0.5),\
                                           ((O_[idx,1-1]**2+O_[idx,2-1]**2)**(0.5)).T* M[idx], \
                                           ((Center[1-1]**2 + Center[2-1]**2)**(0.5)) * math.fsum(M[idx])
    Dibal = Table_l_dis[2-1, -1]
    STD_centers = np.std(Table_l_dis[1, 0:6])
    #STD_centers = np.std(O_[idx, 0],axis=0) + np.std(O_[idx, 2-1],axis=0)
    return Table_l_dis,Center,Dibal,STD_centers
# Посчитать среднее, используя плавающее окно
def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def function_comput_dibal(Table_all,RT_array,T_array,korteg,Center_ost,delt_Z,M,Z):
    """Функция пересчета дисбаланса c помощью смещений"""
    Table_l_dis=np.zeros((2,8,Table_all.shape[1]))
    Center =np.zeros((len(Table_all[0,:]),3))
    Disbal=np.zeros((len(Table_all[0,:])))
    STD_centers=np.zeros((len(Table_all[0,:])))
    Table_l_dis_dist_bal = np.zeros((2, 8, Table_all.shape[1]))
    Center_dist_bal = np.zeros((len(Table_all[0, :]), 3))
    Disbal_dist_bal = np.zeros((len(Table_all[0, :])))
    STD_centers_dist_bal = np.zeros((len(Table_all[0, :])))
    O=np.zeros((7,3))
    X=np.zeros((7))
    Y = np.zeros((7))
    Z_sdv=np.zeros((7))
    for ij in range(len(Table_all[0,:])):
        # %1.5 Смещения вдоль двух других осей (вал задний не сдвигается)
        X[0], X[1:]=0, RT_array[:,4-1,ij].T
        Y[0], Y[1:]=0, RT_array[:,5-1,ij].T
        Z_sdv[0], Z_sdv[1:]=0, (RT_array[:,6-1,ij].T+T_array[:,3-1,ij].T)
        # Сдвиг системы координат (на место первого подшипника)
        O[:,0], O[:,1],O[:,2]= X.T, Y.T, Z.T+Z_sdv.T + delt_Z
        # 1.6 Положение оси, по центрам для подшипников
        O_opora_1, O_opora_2 = np.array([X[1-1], Y[1-1], 0]), np.array([X[7-1], Y[7-1], 609.5 + delt_Z])
        nm_line = np.array((O_opora_2 - O_opora_1) / (math.fsum((O_opora_2 - O_opora_1)**2))**(0.5))
        # 1.7 Вращение системы координат. Поиск дисбалансов
        Table_l_dis[:,:, ij], Center[ij,:],\
        Disbal[ij], STD_centers[ij]=function_of_center_with_support_Z(nm_line,O,O_opora_2,M)
        # Добавление остаточного дисбаланса
        Center_ost_pov1=np.zeros(3)
        Center_ost_pov = np.zeros((7,2))
        X_obsh=np.zeros(7)
        Y_obsh = np.zeros(7)
        O_obsh=np.zeros((7,3))
        for j in range(7):
            if (j>0):
                M_z = np.array([[np.cos(korteg[j-1,ij]), -np.sin(korteg[j-1,ij]), 0],
                                [np.sin(korteg[j-1,ij]), np.cos(korteg[j-1,ij]), 0], [0, 0, 1]])
            else:
                M_z = np.array([[np.cos(0), -np.sin(0), 0],
                                [np.sin(0), np.cos(0), 0], [0, 0, 1]])
            Center_ost_pov1 = np.array([Center_ost[j,0],Center_ost[j,1], 0]).dot(M_z)
            Center_ost_pov[j,:]=Center_ost_pov1[:-1]
            X_obsh[j] = X[j] + Center_ost_pov[j, 0]
            Y_obsh[j] = Y[j] + Center_ost_pov[j, 1]
        O_obsh[:,0],O_obsh[:,1],O_obsh[:,2] = X_obsh.T, Y_obsh.T, Z.T+Z_sdv.T + delt_Z
        #%1.6 Положение оси, по центрам для подшипников
        O_opora_1, O_opora_2 = np.array([X_obsh[0], Y_obsh[0], 0]), np.array([X_obsh[-1], Y_obsh[-1], 609.5 + delt_Z])
        nm_line = np.array((O_opora_2 - O_opora_1) / (math.fsum((O_opora_2 - O_opora_1)**2))**(0.5))
        Table_l_dis_dist_bal[:, :, ij], Center_dist_bal[ij, :], \
        Disbal_dist_bal[ij], STD_centers_dist_bal[ij] = function_of_center_with_support_Z(nm_line,O,O_opora_2,M)

    return Table_l_dis, Center, Disbal,STD_centers,Table_l_dis_dist_bal,Center_dist_bal,Disbal_dist_bal,STD_centers_dist_bal

# Процедура нормирования и вычисление минимумов
def norm_Table_KND(Table_all,korteg,Disbal,k,T_dopusk,T_dopusk_Db):
    Table_KND=np.zeros((Table_all.shape[1],k.shape[0]))
    Table_KND_norm=np.zeros((Table_all.shape[1],k.shape[0]))
    Table_KND = Table_all[k,:].T
    Table_KND_norm = Table_KND/ T_dopusk
    Table_KND_norm_sum=np.sum(Table_KND_norm, axis=1)
    l_geom=np.argmin(Table_KND_norm_sum)
    Angles_result=np.zeros(korteg.shape[0])
    Angles_result = korteg[:, l_geom]
    Disbal_norm=np.zeros(Disbal.shape[0])
    Disbal_norm = Disbal/ T_dopusk_Db
    
    return Table_KND, Table_KND_norm, Table_KND_norm_sum
