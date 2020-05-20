import math
import numpy as np

class Constant_chart():
    """Класс констант для программы chart"""
    def __init__(self):
        # 1.1 Плотность, кг/мм. В виде списка
        self.ro1 = [7.831,7.831,7.831,7.831,7.831,7.831,7.831]
        self.ro=[]
        [self.ro.append(elem/1000000) for elem in self.ro1]
        # Начальные дисбалансы в деталях, гр*мм
        # Начиная с вала заднего (5-7 - придуманы)
        self.D_ish = [64, 200, 175, 4000, 371, 420, 72]
        #Углы начальных дибалансов, со стороны вала заднего если смотреть
        self.ANgle_ish1 = [153, 341, 167, 319, 126, 305, 197]
        self.ANgle_ish=[]
        [self.ANgle_ish.append(elem * math.pi/180+math.pi) for elem in self.ANgle_ish1]
        k=0
        for i in self.ANgle_ish:
            if (i>2*math.pi):
                self.ANgle_ish[k]=i-2*math.pi
            self.ANgle_ish[k]=2*math.pi-self.ANgle_ish[k]
            k=k+1
        # 1.2 Объемы тел
        #   Вал задний; Диск 3; Кольцо 2; Диск 2; Кольцо 1; Диск 1; Вал передний
        self.V=[1051947.911304562, 8596129.423485031, 816257.489545143, 8510707.920516953, 865117.808073511, 13716082.266161982, 634263.681987975]
        self.V_sum=math.fsum(self.V)
        # %1.3 Масса, кг
        self.M=np.zeros(7)
        k=0
        for i,j in zip(self.V,self.ro):
            self.M[k]=i*j
            k=k+1
        self.M_sum = math.fsum(self.M)
        # %Радиус-векторы начальных дисбалансов
        self.ro_ish=[]
        [self.ro_ish.append(i/j/1000) for i,j in zip(self.D_ish,self.M)]
        self.Center_ost = np.ones((len(self.ANgle_ish), 2))
        for i in range(len(self.ANgle_ish)):
            self.Center_ost[i-1,:]=[math.cos(self.ANgle_ish[i-1])*self.ro_ish[i-1], math.sin(self.ANgle_ish[i-1])*self.ro_ish[i-1]]
        #     %1.4 Центроид вдоль оси
        self.Z=np.array([-21.271570925, 56.581426248, 156.500111178, 253.042640978, 366.211814071, 485.713369485, 549.454065462])
        self.delt_Z=89/2
        self.Z_sum01=[]
        [self.Z_sum01.append((i*j)/self.M_sum) for i,j in zip(self.Z,self.M)]
        self.Z_sum0=math.fsum(self.Z_sum01)
        self.Z_sum0=self.Z_sum0+self.delt_Z
        #     %Центроид начальный. В виде кортежа
        self.Center0=(0, 0, self.Z_sum0)
        #Используемые параметры
        self.k =np.array([7, 6, 3, 4, 19, 17, 13, 14, 27, 25, 24, 22])
        self.T_dopusk = 0.1
        self.T_dopusk_Db = 600