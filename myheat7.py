# Have a nice day  
#  Luna


import math
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.animation as animation
#import matplotlib.image as mpimg
import matplotlib as ppl
import time as t
plt.rcParams['animation.ffmpeg_path'] = 'c:/ffmpeg/bin/ffmpeg.exe'
#plt.rcParams['animation.mencoder_path'] = 'i:/mplayer-svn-37847-x86_64/'
#fignumber=0
#import matplotlib
#matplotlib.use("Agg")
#==============================================================================
#==============================================================================
# points -- (ndarray of double, shape (npoints, ndim)) Coordinates of input points.
# vertices -- (ndarray of double, shape (nvertices, ndim)) Coordinates of the 
#             Voronoi vertices.
# ridge_points -- (ndarray of ints, shape (nridges, 2)) Indices of the points 
#                 between which each Voronoi ridge lies.
# ridge_vertices -- (list of list of ints, shape (nridges, *)) Indices of the 
#                   Voronoi vertices forming each Voronoi ridge.
# regions -- (list of list of ints, shape (nregions, *)) Indices of the Voronoi
#            vertices forming each Voronoi region. -1 indicates vertex outside 
#            the Voronoi diagram.
# point_region -- (list of ints, shape (npoints)) Index of the Voronoi region 
#                 for each input point. If qhull option “Qc” was not specified,
#                 the list will contain -1 for points that are not associated 
#                 with a Voronoi region.
#==============================================================================
# b: blue
# g: green
# r: red
# c: cyan
# m: magenta
# y: yellow
# k: black
# w: white
#==============================================================================
#==============================================================================

def euklid_distance(x,y): # расстояние между точками x и y
#    return(sqrt(sum([ (x[i]-y[i])**2 for i in range(len(x)) ])))
    return(math.sqrt(sum([(a-b)**2 for (a,b) in zip(x,y)])))
# Центр масс. На вход: множество список пар координат.
def masscenter(x):
    return list(sum(np.asarray(x))/len(x))
#en.wikipedia.org/wiki/Shoelace_formula
#stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

class Material:
    def __init__(self, hx=0, hy=0, Tmin=0, Tmax=0, x=[], y=[], phasechange=1):
#        global fignumber
        self.hx, self.hy = self.h = (hx, hy) # размеры
        self.Tmin = Tmin            # минимальная температура
        self.Tmax = Tmax            # максимальная температура
        Tplus=self.mu()/self.c()    # повышение температуры при кристаллизации
        self.phasechange = phasechange        
        #YlOrRd_r,hot,autumn
#        self.cmap1 = plt.get_cmap('flag')   # палитра для солидуса
#        self.cmap = plt.get_cmap('flag')    # палитра для ликвидуса
        self.cmap1 = plt.cm.YlOrRd_r
        self.cmap = plt.cm.YlOrRd_r
        # отображение температурного диапазона в [0,1]
        self.tcnorm=plt.Normalize(vmin=min(0,Tmin-0.50*Tplus), vmax=Tmax+Tplus) 
        # отображение [0,1] в палитру для ликвидуса и солидуса
        self.colormap=plt.cm.ScalarMappable(norm=self.tcnorm, cmap=self.cmap)      
        self.colormap1=plt.cm.ScalarMappable(norm=self.tcnorm, cmap=self.cmap1)        
        # формируем границу
        self.num=num=math.sqrt(hx*hy/len(x))    # приблизительная длина стороны квадрата
                                                # в котором в среднем есть одна точка
        d_x=list(np.arange(-num,hx+num,num))    # формирую узлы с шагом num
        d_y=list(np.arange(-num,hy+num,num))    # по осям
        granica=[d_x*2+[-num]*(len(d_y)-2)+[hx+num]*(len(d_y)-2),  # формируем границу
            [-num]*len(d_x)+[hy+num]*len(d_x)+d_y[1:-1]*2]
#        cellsnumber=cellsnumber+len(granica[0])
        # Доавляем точки границы к внутренним точкам
        #x=np.array([list(m.x)+granica[0],list(m.y)+granica[1]])
        self.x, self.y = [list(x)+granica[0],list(y)+granica[1]] # координаты точек
        self.materialcentercell = self.getcellinarea((hx/2-num,hx/2+num),(hy/2-num,hy/2+num))
        self.cellsnumber = len(self.x)
        cellsnumber = len(self.x)
        
        self.voron = Voronoi(np.array([self.x,self.y]).transpose())

        self.T=[Tmax]*cellsnumber       # температурное поле
        self.T1 = [0]*cellsnumber
#        self.TT = [self.T0,self.T1]
        self.neighbours=[]
        self.mass=[math.inf]*cellsnumber
        self.distances_to_neighbours=[]
        self.borders=[0]*cellsnumber    # для каждой точки 0-внутренняя, 1-граничная
        self.circle=[self.Tmin]*cellsnumber
        self.circle[self.materialcentercell]=self.Tmax
        self.contact_surface_area=[]
        self.complete_cell_surface_area=[]
        self.sumcontactareadevideddistance=[math.inf]*cellsnumber # для вычисления условия сходимости
        self.volume=[-1]*cellsnumber    # объем ячейки, по умолчанию -1=infinite
        self.cellphase=[0]*cellsnumber  # 0 - ликвидус, 1 - солидус
        self.complete_cell_surface_area=[math.inf]*cellsnumber # площадь поверхности клетки
        for i in range(cellsnumber):
            self.neighbours.append([])
            self.distances_to_neighbours.append([])
            self.contact_surface_area.append([])
        #  находим для каждой точкие его соседей
        for i in self.voron.ridge_points:
            self.neighbours[i[0]].append(i[1])
            self.neighbours[i[1]].append(i[0])
        # находим для каждой точки расстояния до его соседей и площади поверхностей
        # соприкосновения с ячейками соседей. Соседей перечисляем
        # согласно последовательности списка в neighbours
        for i in range(len(self.neighbours)):
            # Ищем объем ячейки, если она конечная, если бесконечная, то граница
            xx=[]
            yy=[]
            flag=0 # по умолчанию точка внутренняя
            for vert in self.voron.regions[self.voron.point_region[i]]:
                xx.append(self.voron.vertices[vert][0])
                yy.append(self.voron.vertices[vert][1])
                if vert==-1:
                    flag=1 # значит точка граничная
                    self.borders[i]=1 # запомним, что граничная
                    self.T[i]=self.Tmin
            if flag==0:
                self.volume[i]=PolyArea(xx,yy)
                self.mass[i]=self.ro()*self.volume[i]
                # клетка i неграничная -> считаем позже площадь поверхности
                # перед этим площадь math.inf переводим в 0.
                self.complete_cell_surface_area[i]=0
                self.sumcontactareadevideddistance[i]=0
        # закончили вычислять объемы (размеры) ячеек
            for ii in range(len(self.neighbours[i])):
                # ищем расстояния
                self.distances_to_neighbours[i].append(\
                    euklid_distance(\
                    (self.x[i],self.y[i]),\
                    (self.x[self.neighbours[i][ii]],self.y[self.neighbours[i][ii]])\
                    ))
                # ищем площади поверхностей (граней) касания
                # слабое место, поскольку свойство ridge_dict не описано в документации
                # но фактически есть и как раз то, что мне здесь нужно.
                # ridge_dict - словарь, ключь - пары соседних в диаграмме points,
                # значение - пары vertices вершин - концов перпендикулярного ребра
                temp1=(i,self.neighbours[i][ii])
                temp2=(self.neighbours[i][ii],i)
                tmp_vertices = []
                if temp1 in self.voron.ridge_dict:
                    tmp_vertices=self.voron.ridge_dict[temp1]
                elif temp2 in self.voron.ridge_dict:
                    tmp_vertices=self.voron.ridge_dict[temp2]
                if -1 in tmp_vertices: # бесконечно ли ребро
                    self.contact_surface_area[i].append(math.inf)
                else:
                    self.contact_surface_area[i].append(\
                        euklid_distance(\
                        self.voron.vertices[tmp_vertices[0]],\
                        self.voron.vertices[tmp_vertices[1]]\
                        ))
                self.complete_cell_surface_area[i]+=self.contact_surface_area[i][ii]
                self.sumcontactareadevideddistance[i]+=(self.contact_surface_area[i][ii]/self.distances_to_neighbours[i][ii])
            if self.sumcontactareadevideddistance[i]==self.mass[i]==math.inf:
                self.sumcontactareadevideddistance[i]=math.inf
            else:
                self.sumcontactareadevideddistance[i]/=self.mass[i]
            
        self.polygons =[] # строим полигоны, которые поместим потом на рисунок
        self.markphase=[] # точка в центре полигона, маркер ликвидуса/солидуса
        self.materialmass = sum(filter(lambda x:x<math.inf,self.mass)) # фактическая масса слитка
        self.maxcellmass = max(filter(lambda x:x<math.inf,self.mass)) # фактическая ... слитка
        self.mincellmass = min(filter(lambda x:x<math.inf,self.mass)) # фактическая ... слитка
        self.materialvolume = sum(filter(lambda x:x>0,self.volume)) # фактическая ... слитка
        self.mincellvolume = min(filter(lambda x:x>0,self.volume)) # фактическая ... слитка
        self.maxcellvolume = max(filter(lambda x:x>0,self.volume)) # фактическая ... слитка
        self.mindistance=min(list(map(min,self.distances_to_neighbours)))
        self.maxdistance=max(list(map(max,self.distances_to_neighbours)))
        self.maxcellsurfacearea=max(filter(lambda x:x<math.inf,self.complete_cell_surface_area))
        self.mincellsurfacearea=min(filter(lambda x:x<math.inf,self.complete_cell_surface_area))
        suftemp=max(filter(lambda x:x<math.inf,self.sumcontactareadevideddistance))
        self.sufficientstabilitycondition_dt=(self.c()/self.k(0))/suftemp

        # 0.97 чтобы чуточку меньше взять, чем позволено рассуждениями
        self.dt = 0.99 * self.sufficientstabilitycondition_dt # in seconds
        self.life = 0
        self.CellProtocolTime=[]   # здесь запоманаем возраст и температуру
        self.CellProtocolT =[]     # какой-нибудь клетки
        self.voron.close() # больше точек не добавляем
        self.ax = [] 
        self.pathx=[] # для рисунка полигонов
        self.pathy=[] # для рисунка полигонов
        for i in range(len(self.voron.point_region)):
            self.pathx.append([])
            self.pathy.append([])
            if self.voron.regions[self.voron.point_region[i]].count(-1)==0\
                and len(self.voron.regions[self.voron.point_region[i]])>0:
                for ii in self.voron.regions[self.voron.point_region[i]]:
                    self.pathx[i].append(self.voron.vertices[ii][0])
                    self.pathy[i].append(self.voron.vertices[ii][1])
        self.flag_fig=0
#        fignumber=fignumber+1
#        self.fig=plt.figure(fignumber, figsize=(10,10))
        self.flag_printtcircle=0
#        self.flag_nexttime=0
        self.flag_plotpolygon=0
        self.ax_printtcircle=0
        self.ax_nexttime=0
        self.count=0    # для тестирования, не нужно
        self.count1=0    # для тестирования, не нужно
#        self.alpha0 = 1 # прозрачность ликвидус, уже не нужно
#        self.alpha1 = 1 # прозрачнось солидус # непрозрачный, уже не нужно
#        self.ims=[]
#        FFMpegWriter = animation.writers['ffmpeg']
        metadata = dict(title='Movie Test', artist='Matplotlib',comment='Movie support!')
#        self.writer = animation.FFMpegWriter(fps=1, metadata=metadata, extra_args=['-vcodec', 'libx264'])
        self.writer = animation.FFMpegWriter(fps=1, metadata=metadata, bitrate=1800)
#        self.writer.setup(self.fig,'my.mp4',100)
# Вызывается только один раз
    def initpolygons(self):
        for i in range(len(self.pathx)):
            self.markphase.append([])
            self.polygons.append([])
            if len(self.pathx[i])>0:
                self.polygons[i]=\
                patches.Polygon(np.asarray([self.pathx[i],self.pathy[i]]).transpose(),\
                color=self.colormap.to_rgba(self.T[i]))
                self.markphase[i]=plt.Line2D\
                ([self.voron.points[i][0]],[self.voron.points[i][1]],marker=',',color='white')
#                self.polygons[i].set_alpha(self.alpha0)
#        self.flag_plotpolygon=1 # ставим флаг
        self.ax_nexttime=self.initfig() # инициализируем рисунок (axe)
        for i in range(len(self.polygons)):
            if self.polygons[i]!=[]:
                self.ax_nexttime.add_patch(self.polygons[i])
                self.ax_nexttime.add_line(self.markphase[i])
#        ppl.colorbar.ColorbarBase(self.ax_nexttime,cmap=self.cmap, norm=self.tcnorm)
    def initfig(self):
#        global fignumber
        if self.flag_fig==0:
            self.flag_fig=1
            temp=plt.get_fignums()
            temp.append(0)
            fignumber=max(temp)+1
            self.fig = plt.figure(fignumber, figsize=(10,10))
            self.fig.clf()
        ax = self.fig.add_axes([0.05, 0.05, 0.7,0.9],axisbg='w',aspect='equal')
        ax.set_xlim(-self.num,hx+self.num)
        ax.set_ylim(-self.num,hy+self.num)
        axbar=self.fig.add_axes([0.8, 0.05, 0.1,0.9])
#        axbar.set_xlim(0.8,0.1)
#        axbar.set_ylim(-self.num,hy+self.num)
#        mpl.colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='vertical')
        ppl.colorbar.ColorbarBase(axbar,cmap=self.cmap, norm=self.tcnorm)
        return(ax)        

    def getcellinarea(self,o1,o2):
        for i in range(len(self.x)):
            if o1[0] < self.x[i] < o1[1] and o2[0] < self.y[i] < o2[1]:
                return(i)
        return(-1)

    def Tnext(self,cell): # Meltingpoint = MP
        if self.flag_plotpolygon==0:
            self.flag_plotpolygon=1
            self.initpolygons()
        if self.borders[cell]:
            return(self.T[cell],self.cellphase[cell])
        mass=self.mass[cell] 
        phase=self.cellphase[cell]
        changephase=0
        dQ=self.dt*sum(
            list(map(
            (lambda x: self.k((self.T[self.neighbours[cell][x]]+self.T[cell])/2)*
            self.contact_surface_area[cell][x]*
            (self.T[self.neighbours[cell][x]]-self.T[cell])/
            self.distances_to_neighbours[cell][x]),
            range(len(self.neighbours[cell]))
            )))
        TnextXX=self.T[cell]+dQ/(self.c()*mass)
        # если работаем без смены фаз, то выходим
        if self.phasechange == 0: 
            return(TnextXX,0)
        # Если температура клетки ниже температуры плавления,
        # а её фаза сейчас ликвидус, вычислим next температуру в случае
        # кристаллизации
        if self.T[cell] < self.melting_point() and phase == 0:
            flag = np.random.uniform(TnextXX,TnextXX + self.mu()/self.c()) # для выбора: осуществлять ли фазовый переход
            if  flag < self.melting_point():
                changephase=1                
                phase = 1
#                self.polygons[cell].set_alpha(self.alpha1)
#                self.count=self.count+1
        # Если температура клетки выеше температуры плавления,
        # а её фаза сейчас солидус, вычислим next температуру в случае
        # плавления
        if self.T[cell] > self.melting_point() and phase == 1:
            flag = np.random.uniform(TnextXX - self.mu()/self.c(),TnextXX) # для выбора: осуществлять ли фазовый переход
            if  flag > self.melting_point():
                changephase=-1                
                phase = 0
#                self.polygons[cell].set_alpha(self.alpha0)
#                self.count1=self.count1+1
        self.cellphase[cell]=phase
        # возвращаем температуру и фазу
        return(TnextXX+changephase*self.mu()/self.c(),self.cellphase[cell])

    def nexttime(self,times):
        self.life = self.life + self.dt
        for i in range(times):
            for cell in range(len(self.T)):
                self.T1[cell]=self.Tnext(cell)[0]
            change=self.T
            self.T=self.T1
            self.T1=change
#            for cell in range(len(self.T)):
#                self.T[cell]=self.T1[cell]

    # протоколируем возраст и температуру клетки cell            
    def writeTdown(self,cell):
        self.CellProtocolTime.append(self.life)
        self.CellProtocolT.append(self.T[cell])

    # Раскрашиваем полигоны согласно температуре. Отмечаем цветом точки фазу    
    def plotpolygon(self):
        if self.flag_plotpolygon==0:
            self.flag_plotpolygon=1
            self.initpolygons()
        for cell in range(len(self.T)):
            if self.polygons[cell]!=[]:
                if self.cellphase[cell]==0:
                    self.polygons[cell].set_color(self.colormap.to_rgba(self.T[cell]))
                    self.markphase[cell].set_color('white')                    
                else:
                    self.polygons[cell].set_color(self.colormap1.to_rgba(self.T[cell]))
                    self.markphase[cell].set_color('black')                    
    # Указываем лишь только фазу
    def plotphase(self):
        if self.flag_plotpolygon==0:
            self.flag_plotpolygon=1
            self.initpolygons()
        for cell in range(len(self.T)):
            if self.polygons[cell]!=[]:
                if self.cellphase[cell]==0:
                    self.polygons[cell].set_color('white')
                else:
                    self.polygons[cell].set_color('black')

    def printtcircle(self,n):
        if self.flag_printtcircle==0:
            self.flag_printtcircle=1
            self.ax_printtcircle=self.initfig()
        temp=[]
        for i in range(n):
            for ii in range(len(self.circle)):
                for cell in self.neighbours[ii]:
                    if self.circle[cell]==self.Tmax:
                        temp.append(ii)
            for ii in temp:
                self.circle[ii]=self.Tmax
#        self.ax.clear()
#        self.ax.scatter(self.x,self.y,c=self.circle,edgecolors='None',cmap=self.cmap,marker='o',s=5)
        self.plotpolygon2(self.ax_printtcircle)

    def initT(self):
        for i in range(len(self.neighbours)):
            self.T[i]=self.Tmin
            for vert in self.voron.regions[self.voron.point_region[i]]:
                if vert==-1:
                    self.T[i]=self.Tmax

    def distance_cells(self,c1,c2):
        return(euklid_distance((self.x[c1],self.y[c1]),(self.x[c2],self.y[c2])))

    def ro(self): # [кг/м^3] плотность вещества
        return(7870)
    def c(self): # удельная теплоемкость
        # [Дж/(кг*К)] удельная теплоемкость
        return(300)
    def k(self,temp): # коэффициент теплопроводности
        # [Дж/(м*К*с)]=[Вт/(м*К)]=[кг*м^2/c^3] 50-100
        if 0<=temp<=800:
            ktemp=((temp-0)/(800-0))*(14-60)+60
        elif 800<=temp<=1520:
            ktemp=((temp-800)/(1520-800))*(34-14)+14
        elif 1520<=temp:
            ktemp=35
        return(ktemp)
    def mu(self): # теплота кристаллизации
        # [Дж/кг]
        return(220000) # сталь
    def melting_point(self):
        # [K] 1450—1520 °C для стали
        return(1500+273.15)

    def xycell(self,cell):
        return((self.x[cell],self.y[cell]))
# находит центр масс для конечного полигона  
# (ячейки, клетки и т.п. синонимов)        
    def masscenterpoints(self):
        xnew=[]
        for reg in self.voron.regions:
            tmp=[]
            if reg.count(-1)==0 and len(reg)>0:
                for i in reg:
                    tmp.append(self.voron.vertices[i])
                xnew.append(masscenter(tmp))
        return(xnew)

    def plotpolygon2(self,ax):
        ax.clear()
        for i in range(len(self.pathx)):
            if len(self.pathx[i])>0:
                ax.fill(self.pathx[i], self.pathy[i], color=self.colormap.to_rgba(self.circle[i]))        

    #[12:54:01] Anna Ivanova: поступающий расплав 1823К, граница в верхней части 1773К, в нижней части 1178К
    #[12:54:32] Anna Ivanova: высота кр-ра 0,9 м
    #[12:54:46] Anna Ivanova: вернее слитка в кристаллизаторе
    # толщина слитка 20 см, толщина расплава 1 см, высота слитка 0.9 м.

random.seed(3)

hx = 1. # Metre
hy = 1. # Metre
#Tmin = 1178 # Kelvine
Tmin=100
Tupper = 1773 # Kelvine
Tmax = 1823 # Kelvine
phasechange=1
#dm = 0.2 # толщина расплава
cellsnumber = 50**2
x = np.random.uniform(0,hx,cellsnumber)
y = np.random.uniform(0,hy,cellsnumber)
m=[]
m.append(Material(hx,hy,Tmin,Tmax,x,y,phasechange))
for i in range(1,5):
    xnew = np.array(m[len(m)-1].masscenterpoints()).transpose()
    m.append(Material(hx,hy,Tmin,Tmax,xnew[0],xnew[1],phasechange))

mat = m[len(m)-1]
del m

def savemoving():
    with mat.writer.saving(mat.fig,'my1A.mp4',100):
       for i in range(20):
            start = t.time()
            print(i)
            mat.writer.grab_frame()
            mat.nexttime(1)
            mat.plotpolygon()
            end = t.time()
            print(end-start)

# plottemperature(1,30,50) цикл от 1 до 30 по 50 единиц времени
def plottemperature(x,y,z):
    for i in range(x,y):
        print(i)
        mat.nexttime(z)
        mat.plotpolygon()
        fname='temptemp'+str(10**8+(i+1)*z)+'.jpeg'
        mat.fig.savefig(fname)
def plotphase(x,y,z):
    for i in range(x,y):
        print(i)
        mat.nexttime(z)
        mat.plotphase()
        fname='tempphase'+str(10**8+(i+1)*z)+'.jpeg'
        mat.fig.savefig(fname)

# Регулярная сетка
#==============================================================================
num=math.sqrt(hx*hy/cellsnumber)
xreg = np.linspace(0,hx,hx/num)
yreg = np.linspace(0,hy,hy/num)
x,y = np.meshgrid(xreg,yreg)
x.shape=(1,len(xreg)*len(yreg))
y.shape=(1,len(xreg)*len(yreg))
x=x[0]
y=y[0]
matreg=Material(hx,hy,Tmin,Tmax,x,y,phasechange)
matreg.dt=mat.dt
#==============================================================================

'''
for i in range(1200):
    print(i)
    mat.nexttime(10)
    matreg.nexttime(10)
    mat.writeTdown(mat.materialcentercell)
    matreg.writeTdown(matreg.materialcentercell)

temp=plt.get_fignums()
temp.append(0)
fignumber = max(temp)+1
fig = plt.figure(fignumber, figsize=(10,10))
ax = fig.add_axes([0.05, 0.05, 0.7,0.9],axisbg='w',aspect='equal')
#fig.xlabel('Time')
#fig.ylabel('Temperature')
#fig.title('Graphic of temperature')
#fig.grid(True)
x1=mat.CellProtocolTime
y1=mat.CellProtocolT
x2=matreg.CellProtocolTime
y2=matreg.CellProtocolT
#ax0 = fig.add_axes([0.05, 0.05, 0.7,0.9],axisbg='w',aspect='equal')
#ax1 = fig.add_axes([0.05, 0.05, 0.7,0.9],axisbg='w',aspect='equal')
fig.plot(x1, y1, 'r,',x2, y2, 'b,')
ax.plot(x1, y1, 'r,',x2, y2, 'b,')
#ax.set_xlim(-self.num,hx+self.num)
#ax.set_ylim(-self.num,hy+self.num)
'''



def savemovingforreg():
    global matreg
    with matreg.writer.saving(matreg.fig,'my1B.mp4',100):
       for i in range(20):
            start = t.time()
            print(i)
            matreg.writer.grab_frame()
            matreg.nexttime(1)
            matreg.plotpolygon()
            end = t.time()
            print(end-start)
def plotpolygonforreg():
    global matreg
    for i in range(30):
        matreg.nexttime(500)
        matreg.plotpolygon()
        fname='D_Reg'+str((i+1)*500)+'.jpeg'
        matreg.fig.savefig(fname)



'''
mat.nexttime(1), min(mat.T), max(mat.T), mat.cellphase.count(1), mat.plotpolygon()


for i in range(50):
    mat.nexttime(600)
    mat.plotpolygon()
    fname=str(i)+'.jpeg'
    mat.fig.savefig(fname)


'''

'''
img=[]
fig=plt.figure()
for i in range(286):
    filename='.\\b\\'+str(i)+'.jpeg'
    img.append([plt.imshow(mpimg.imread(filename))])

ani = animation.ArtistAnimation(fig, mat.ims, interval=50, blit=True,repeat_delay=1000)

FFwriter = animation.FFMpegWriter()
ani.save('teplo.mp4', writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'])
'''


# mat.nexttime(10000), (min(mat.T),max(mat.T))



#mat.initfig()
#mat.nexttime(1)

#for i in m:
#    voronoi_plot_2d(i.voron)


#
#for i in m:
#    i.printtcircle(30)    

#m.printtcircle(1)
#plt.plot(m.x,m.y,marker='o', ls='None')
#voronoi_plot_2d(m.voron)