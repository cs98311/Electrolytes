#!/usr/bin/env python
# coding: utf-8

# #### Import

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')
#%matplotlib qt
#get_ipython().run_line_magic('reload_ext', 'autoreload')
#get_ipython().run_line_magic('autoreload', '2')


# In[2]:


from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from numpy import NaN
import os



# #### Import System Info

# In[3]:


sys_info = np.loadtxt('SysInfo.txt')
box_length, numW, numT, numL, numZ = sys_info


# #### Converting float to int 

# In[4]:


numW = numW.astype(int)
numT = numT.astype(int)
numL = numL.astype(int)
numZ = numZ.astype(int)


# #### Coordinates Import 

# In[5]:


# Water coordinates
OW = np.loadtxt('Coordinates/OW.txt', dtype= 'float')
H1 = np.loadtxt('Coordinates/H1.txt', dtype= 'float')
H2 = np.loadtxt('Coordinates/H2.txt', dtype= 'float')

# Cation Coordinates
# ZN = np.loadtxt('Coordinates/ZN.txt', dtype= 'float')
LI = np.loadtxt('Coordinates/LI.txt', dtype= 'float')

# TFSI Coordinates
C1 = np.loadtxt('Coordinates/C1.txt', dtype= 'float')
C2 = np.loadtxt('Coordinates/C2.txt', dtype= 'float')

S1 = np.loadtxt('Coordinates/S1.txt', dtype= 'float')
S2 = np.loadtxt('Coordinates/S2.txt', dtype= 'float')

N1 = np.loadtxt('Coordinates/N1.txt', dtype= 'float')

O1 = np.loadtxt('Coordinates/O1.txt', dtype= 'float')
O2 = np.loadtxt('Coordinates/O2.txt', dtype= 'float')
O3 = np.loadtxt('Coordinates/O3.txt', dtype= 'float')
O4 = np.loadtxt('Coordinates/O4.txt', dtype= 'float')

F1 = np.loadtxt('Coordinates/F1.txt', dtype= 'float')
F2 = np.loadtxt('Coordinates/F2.txt', dtype= 'float')
F3 = np.loadtxt('Coordinates/F3.txt', dtype= 'float')
F4 = np.loadtxt('Coordinates/F4.txt', dtype= 'float')
F5 = np.loadtxt('Coordinates/F5.txt', dtype= 'float')
F6 = np.loadtxt('Coordinates/F6.txt', dtype= 'float')

# Combined Fluorine and Oxygen Coordinates
F = np.loadtxt('Coordinates/F.txt', dtype= 'float')
Ot = np.loadtxt('Coordinates/Ot.txt', dtype= 'float')


# In[6]:


HB_water = np.loadtxt('Plots/DataFiles/HB_with_water.txt',dtype='int')
Zero_HB = HB_water==0
One_HB = HB_water==1
Two_HB = HB_water==2
Three_HB = HB_water==3
Four_HB = HB_water==4
Five_HB = HB_water==5


# ##### Check whether file is empty

# In[7]:


#os.stat('LI.txt').st_size == 0


# # <span style="color:RoyalBlue">H2O System</span></span>

# ### Scatterplot 

# In[8]:


x1 = OW[:,0]
y1 = OW[:,1]
z1 = OW[:,2]

fig1 = plt.figure(figsize=(19.20,10.80),dpi=300)

axes1 = fig1.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes1.scatter3D(x1,y1,z1,
                #c=z1,
                #cmap='Blues',
                s = 50,
                label = r'$H_2O$',
                color = 'blue',
                edgecolor = 'k'
               )

#axes1.view_init(elev=0, azim=0)

axes1.legend(loc='upper right')

plt.savefig('Plots/Water/H2Osys.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/H2Osys.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/H2Osys.jpg',format='jpg',bbox_inches='tight')
plt.clf()


# #### No-Axes Scatterplot

# In[9]:


fig1 = plt.figure(figsize=(19.60,10.80),dpi=300)

axes1 = fig1.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes1.scatter3D(x1,y1,z1,
                #c=z1,
                #cmap='Blues',
                s = 50,
                label = r'$H_2O$',
                color = 'blue',
                edgecolor = 'k'
               )
axes1.set_axis_off()

#axes1.view_init(elev=0, azim=0)

axes1.legend(loc='upper right')

plt.savefig('Plots/Water/H2OsysNoAxes.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/H2OsysNoAxes.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/H2OsysNoAxes.jpg',format='jpg',bbox_inches='tight')
plt.clf()

#plt.savefig('H2Osys.svg',format='svg',bbox_inches='tight')


# #### Slice of Scatterplot

# In[10]:


range1 = (x1>1)*(x1<2)


# In[11]:


fig1 = plt.figure(figsize=(19.60,10.80),dpi=300)

axes1 = fig1.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes1.scatter3D(x1[range1],y1[range1],z1[range1],
                #c=z1,
                #cmap='Blues',
                s = 50,
                label = r'$H_2O$',
                color = 'blue',
                edgecolor = 'k'
               )

axes1.set_xlim3d([0,box_length])
axes1.set_ylim3d([0,box_length])
axes1.set_zlim3d([0,box_length])


axes1.view_init(elev=0, azim=0)

axes1.set_axis_off()

axes1.legend(loc='upper right')

plt.savefig('Plots/Water/H2OsysSlice.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/H2OsysSlice.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/H2OsysSlice.jpg',format='jpg',bbox_inches='tight')
plt.clf()

#plt.savefig('H2Osys.svg',format='svg',bbox_inches='tight')


# #### HB Vectors (Oxy to Oxy) 

# In[12]:


data_o2o = np.loadtxt('Plots/DataFiles/HBvectorsO2O.txt', dtype = 'float', delimiter = ',',unpack=False)


# In[13]:


a1 = data_o2o[:,0]
b1 = data_o2o[:,1]
c1 = data_o2o[:,2]
d1 = data_o2o[:,3]
e1 = data_o2o[:,4]
f1 = data_o2o[:,5]

# a1,....,f1 =zip(*data_o2o)  also same thing


# In[14]:


fig_o2o = plt.figure(figsize=(19.60,10.80),dpi=300)

axes_o2o = fig_o2o.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes_o2o.quiver(a1,b1,c1,d1,e1,f1, 
                label = r'$H_2O-H_2O$',
                color = 'blue',
                arrow_length_ratio=0.4,
                #lw = 1.5
               )

#axes_o2o.set_xlim3d([0, 1])
#axes_o2o.set_ylim3d([0, 1])
#axes_o2o.set_zlim3d([0, 1])

axes_o2o.view_init(
                   elev=10, 
                   azim=0
                  )

#plt.set_aspect('auto')
#axes_o2o.set_axis_off()

#fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)

plt.xlabel('x')
plt.ylabel('y')

axes_o2o.legend(loc='upper right')

plt.savefig('Plots/Water/HBvectors.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/HBvectors.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/HBvectors.jpg',format='jpg',bbox_inches='tight')
plt.clf()


# ### Scatter + HB Vectors 

# In[15]:


fig_o2o = plt.figure(figsize=(19.60,10.80),dpi=300)

axes_o2o = fig_o2o.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes_o2o.quiver(a1,b1,c1,d1,e1,f1, 
                label = r'$H_2O-H_2O$',
                color = 'blue',
                arrow_length_ratio=0.4
               )

axes_o2o.scatter3D(x1,y1,z1,
                #c=z1,
                #cmap='Blues',
                s = 50,
                label = r'$H_2O$',
                color = 'cyan',
                edgecolor='k'
               )

plt.savefig('Plots/Water/VectorsPlusScatter.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/VectorsPlusScatter.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/VectorsPlusScatter.jpg',format='jpg',bbox_inches='tight')
plt.clf()


# #### Slice of HB Vectors (Oxy to Oxy)

# In[16]:


range_o2o = (a1<3.5)*(a1>2.25)


# In[17]:


fig_o2o = plt.figure(figsize=(19.60,10.80),dpi=300)

axes_o2o = fig_o2o.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes_o2o.quiver(a1[range_o2o],b1[range_o2o],c1[range_o2o],d1[range_o2o],e1[range_o2o],f1[range_o2o], 
                label = r'$H_2O-H_2O$',
                color = 'blue',
                arrow_length_ratio=0.4
               )

axes_o2o.set_xlim3d([0, 3.5])
axes_o2o.set_ylim3d([0, 3.5])
axes_o2o.set_zlim3d([0, 3.5])

axes_o2o.view_init(
                   elev=10, 
                   azim=0
                  )

axes_o2o.set_aspect('auto')
axes_o2o.set_axis_off()

#fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)

plt.xlabel('x')
plt.ylabel('y')

axes_o2o.legend(loc='upper right')

plt.savefig('Plots/Water/HBvectorsSlice.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/HBvectorsSlice.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/HBvectorsSlice.jpg',format='jpg',bbox_inches='tight')
plt.clf()



# In[18]:


#plt.clf()


# In[19]:


"""
# Lengthy way of doing the same
a_t = a1
b_t = b1
c_t = c1
d_t = d1
e_t = e1
f_t = f1

# Making temp variables

for i in np.arange(len(a1)):
    if (a1[i]>3.5)+(a1[i]<2.25):
        a_t[i] = NaN
        b_t[i] = NaN
        c_t[i] = NaN
        d_t[i] = NaN
        e_t[i] = NaN
        f_t[i] = NaN

fig_o2o = plt.figure(figsize=(5,5),dpi=300)

axes_o2o = fig_o2o.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes_o2o.quiver(a_t,b_t,c_t,d_t,e_t,f_t, 
                label = r'$H_2O-H_2O$',
                color = 'blue',
                arrow_length_ratio=0.4
               )

axes_o2o.set_xlim3d([0, 3.5])
axes_o2o.set_ylim3d([0, 3.5])
axes_o2o.set_zlim3d([0, 3.5])

axes_o2o.view_init(
                   elev=10, 
                   azim=0
                  )

axes_o2o.set_aspect('auto')
#axes_o2o.set_axis_off()

#fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)

plt.xlabel('x')
plt.ylabel('y')

axes_o2o.legend(loc='upper right')
"""


# ### HB Water 

# In[20]:


x1 = OW[:,0]
y1 = OW[:,1]
z1 = OW[:,2]

fig1 = plt.figure(figsize=(19.20,10.80),dpi=300)

axes1 = fig1.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes1.scatter3D(x1[Zero_HB],y1[Zero_HB],z1[Zero_HB],s = 50,label = r'$0 HB$',color = 'Red',edgecolor = 'k')
axes1.scatter3D(x1[One_HB],y1[One_HB],z1[One_HB],s = 50,label = r'$1 HB$',color = 'blue',edgecolor = 'k')
axes1.scatter3D(x1[Two_HB],y1[Two_HB],z1[Two_HB],s = 50,label = r'$2 HB$',color = 'magenta',edgecolor = 'k')
axes1.scatter3D(x1[Three_HB],y1[Three_HB],z1[Three_HB],s = 50,label = r'$3 HB$',color = 'green',edgecolor = 'k')
axes1.scatter3D(x1[Four_HB],y1[Four_HB],z1[Four_HB],s = 50,label = r'$4 HB$',color = 'cyan',edgecolor = 'k')
axes1.scatter3D(x1[Five_HB],y1[Five_HB],z1[Five_HB],s = 500,label = r'$5 HB$',color = 'lime',edgecolor = 'k')

#axes1.view_init(elev=0, azim=0)
axes1.set_axis_off()

axes1.view_init(elev=8, azim=-140)

axes1.legend(loc='upper right')

plt.savefig('Plots/Water/HBcolored.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/HBcolored.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/HBcolored.jpg',format='jpg',bbox_inches='tight')
plt.clf()


### Colored + Vectors

# In[46]:

x1 = OW[:,0]
y1 = OW[:,1]
z1 = OW[:,2]

fig1 = plt.figure(figsize=(19.20,10.80),dpi=300)

axes1 = fig1.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes1.scatter3D(x1[Zero_HB],y1[Zero_HB],z1[Zero_HB],s = 50,label = r'$0 HB$',color = 'Red',edgecolor = 'k')
axes1.scatter3D(x1[One_HB],y1[One_HB],z1[One_HB],s = 50,label = r'$1 HB$',color = 'blue',edgecolor = 'k')
axes1.scatter3D(x1[Two_HB],y1[Two_HB],z1[Two_HB],s = 50,label = r'$2 HB$',color = 'magenta',edgecolor = 'k')
axes1.scatter3D(x1[Three_HB],y1[Three_HB],z1[Three_HB],s = 50,label = r'$3 HB$',color = 'green',edgecolor = 'k')
axes1.scatter3D(x1[Four_HB],y1[Four_HB],z1[Four_HB],s = 50,label = r'$4 HB$',color = 'cyan',edgecolor = 'k')
axes1.scatter3D(x1[Five_HB],y1[Five_HB],z1[Five_HB],s = 500,label = r'$5 HB$',color = 'lime',edgecolor = 'k')

axes1.quiver(a1,b1,c1,d1,e1,f1, 
                label = r'$H_2O-H_2O$',
                color = 'navy',
                arrow_length_ratio=0.4
               )

axes1.view_init(elev=8, azim=-140)

axes1.set_axis_off()

axes1.legend(loc='upper right')

plt.savefig('Plots/Water/HBcolored+Vectors.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Water/HBcolored+Vectors.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Water/HBcolored+Vectors.jpg',format='jpg',bbox_inches='tight')
plt.clf()



#plt.savefig('H2Osys.svg',format='svg',bbox_inches='tight')


# # <span style="color:MediumSeaGreen">TFSI </span></span>

# ### COM Scatterplot

# In[21]:


data_com = np.loadtxt('Plots/DataFiles/TFSIcom.txt')


# In[22]:


x_com = data_com[:,0]
y_com = data_com[:,1]
z_com = data_com[:,2]


# In[23]:


fig_com = plt.figure(figsize=(19.60,10.80),dpi=300)

axes_com = fig_com.add_axes([0.1,0.1,0.9,0.9], projection='3d')

axes_com.scatter3D(x_com,y_com,z_com, s=50, label = 'TFSI Centre of Mass', color = 'green',edgecolor = 'k')

#axes1.view_init(elev=0, azim=0)
#axes_o2o.set_axis_off()

plt.legend(loc='upper right')

plt.savefig('Plots/TFSI/TFSIcom.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSIcom.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSIcom.jpg',format='jpg',bbox_inches='tight')
plt.clf()


#plt.savefig('TFSIcom.svg',format='svg',bbox_inches='tight')


# ### Vector TFSI Skeleton

# ##### Function to take care of long cross-box vectors due to PBC

# In[24]:


def pbc_issue(alk):
    for i in np.arange(len(alk)):
        for j in range(0,3):
            if(np.absolute(alk[i,j])>1):
                alk[i,j]=box_length-np.absolute(alk[i,j])


# ##### Direction vectors for skeleton N,S,C

# In[25]:


v1 = S1 - N1
v2 = S2 - N1
v3 = C1 - S1
v4 = C2 - S2


# ##### Direction vectors for F

# In[26]:


v5 = F1 - C1
v6 = F2 - C1
v7 = F3 - C1
v8 = F4 - C2
v9 = F5 - C2
v10 = F6 - C2


# ###### Direction vectors for O

# In[27]:


v11 = O1 - S1
v12 = O2 - S1
v13 = O3 - S2
v14 = O4 - S2


# In[28]:


pbc_issue(v1)
pbc_issue(v2)
pbc_issue(v3)
pbc_issue(v4)
pbc_issue(v5)
pbc_issue(v6)
pbc_issue(v7)
pbc_issue(v8)
pbc_issue(v9)
pbc_issue(v10)
pbc_issue(v11)
pbc_issue(v12)
pbc_issue(v13)
pbc_issue(v14)


# In[29]:


# Function to print vectors by giving axes, start point, direction vector, color, label, arrow length ratio as input

def plot_vectors(axes,rng,start,dirn,clr,lbl,ratio,lwd):
    axes.quiver(
                start[:,0][rng],start[:,1][rng],start[:,2][rng],dirn[:,0][rng],dirn[:,1][rng],dirn[:,2][rng], 
                color = clr,
                label = lbl,
                arrow_length_ratio=ratio,
                lw =lwd
               )


# In[30]:


fig_sk = plt.figure(figsize = (19.60,10.80), dpi=300)

axes_sk = fig_sk.add_axes([0.1,0.1,0.9,0.9], projection ='3d')

rng = np.ones(numT)==1

plot_vectors(axes_sk,rng,N1,v1,'lime','',0,1)
plot_vectors(axes_sk,rng,N1,v2,'lime','',0,1)
plot_vectors(axes_sk,rng,S1,v3,'lime','',0,1)
plot_vectors(axes_sk,rng,S2,v4,'lime','',0,1)

plot_vectors(axes_sk,rng,C1,v5,'red','',0,1)
plot_vectors(axes_sk,rng,C1,v6,'red','',0,1)
plot_vectors(axes_sk,rng,C1,v7,'red','',0,1)

plot_vectors(axes_sk,rng,C2,v8,'red','',0,1)
plot_vectors(axes_sk,rng,C2,v9,'red','',0,1)
plot_vectors(axes_sk,rng,C2,v10,'red','',0,1)

plot_vectors(axes_sk,rng,S1,v11,'blue','',0,1)
plot_vectors(axes_sk,rng,S1,v12,'blue','',0,1)
plot_vectors(axes_sk,rng,S2,v13,'blue','',0,1)
plot_vectors(axes_sk,rng,S2,v14,'blue','',0,1)

axes_sk.view_init(
                   elev=10, 
                   azim=0
                  )

#axes_sk.set_axis_off()
#plt.legend(loc='upper right')

plt.savefig('Plots/TFSI/TFSIskeleton.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSIskeleton.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSIskeleton.jpg',format='jpg',bbox_inches='tight')
plt.clf()


# ### Slice of TFSI Skeleton 

# ##### Set range for slice 

# In[31]:


rng = (N1[:,0]<2)*(N1[:,0]>1.5)


# In[32]:


##### To Print a single TFSI 

"""
single = np.ones(144)
single[18]=0
rng = single==0
"""


# In[33]:


fig_sk = plt.figure(figsize = (19.60,10.80), dpi=300)

axes_sk = fig_sk.add_axes([0.1,0.1,0.9,0.9], projection ='3d')

plot_vectors(axes_sk,rng,N1,v1,'lime','N-S',0,1)
plot_vectors(axes_sk,rng,N1,v2,'lime','',0,1)
plot_vectors(axes_sk,rng,S1,v3,'lime','S-C',0,1)
plot_vectors(axes_sk,rng,S2,v4,'lime','',0,1)

plot_vectors(axes_sk,rng,C1,v5,'red','C-F',0,1)
plot_vectors(axes_sk,rng,C1,v6,'red','',0,1)
plot_vectors(axes_sk,rng,C1,v7,'red','',0,1)

plot_vectors(axes_sk,rng,C2,v8,'red','',0,1)
plot_vectors(axes_sk,rng,C2,v9,'red','',0,1)
plot_vectors(axes_sk,rng,C2,v10,'red','',0,1)

plot_vectors(axes_sk,rng,S1,v11,'gold','S-O',0,4)
plot_vectors(axes_sk,rng,S1,v12,'gold','',0,4)
plot_vectors(axes_sk,rng,S2,v13,'gold','',0,4)
plot_vectors(axes_sk,rng,S2,v14,'gold','',0,4)

axes_sk.set_xlim3d([0, box_length])
axes_sk.set_ylim3d([0, box_length])
axes_sk.set_zlim3d([0, box_length])

axes_sk.view_init(
                   elev=10, 
                   azim=0
                  )

axes_sk.set_aspect('auto')
#axes_sk.set_axis_off()

plt.legend(loc = 'right')

plt.savefig('Plots/TFSI/TFSIskeletonSlice.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSIskeletonSlice.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSIskeletonSlice.jpg',format='jpg',bbox_inches='tight')
plt.clf()


# ### With isolated water and cation

# In[36]:


data_cat = np.loadtxt('Coordinates/LI.txt', dtype='float')
x_cat = data_cat[:,0]
y_cat = data_cat[:,1]
z_cat = data_cat[:,2]


# In[37]:


range_iso = (OW[Zero_HB][:,0]>1.5)*(OW[Zero_HB][:,0]<2) 
range_iso2 = (LI[:,0]>1.5)*(LI[:,0]<2)
rng = (N1[:,0]<2)*(N1[:,0]>1.5)


# In[38]:


fig_sk = plt.figure(figsize = (19.60,10.80), dpi=300)

axes_sk = fig_sk.add_axes([0.1,0.1,0.9,0.9], projection ='3d')

plot_vectors(axes_sk,rng,N1,v1,'green','',0,1)
plot_vectors(axes_sk,rng,N1,v2,'green','',0,1)
plot_vectors(axes_sk,rng,S1,v3,'green','',0,1)
plot_vectors(axes_sk,rng,S2,v4,'green','',0,1)

plot_vectors(axes_sk,rng,C1,v5,'red','C-F',0,1)
plot_vectors(axes_sk,rng,C1,v6,'red','',0,1)
plot_vectors(axes_sk,rng,C1,v7,'red','',0,1)

plot_vectors(axes_sk,rng,C2,v8,'red','',0,1)
plot_vectors(axes_sk,rng,C2,v9,'red','',0,1)
plot_vectors(axes_sk,rng,C2,v10,'red','',0,1)

plot_vectors(axes_sk,rng,S1,v11,'green','',0,1)
plot_vectors(axes_sk,rng,S1,v12,'green','',0,1)
plot_vectors(axes_sk,rng,S2,v13,'green','',0,1)
plot_vectors(axes_sk,rng,S2,v14,'green','',0,1)

axes_sk.scatter3D(OW[:,0][Zero_HB][range_iso],OW[:,1][Zero_HB][range_iso],OW[:,2][Zero_HB][range_iso], s = 200, color = 'cyan', label = r'$H_2O(iso)$', edgecolor = 'k')

axes_sk.scatter3D(x_cat[range_iso2],y_cat[range_iso2],z_cat[range_iso2], s= 50, color = 'lime', label = r'$Li^{+2}$', edgecolor = 'k')


axes_sk.set_xlim3d([0, box_length])
axes_sk.set_ylim3d([0, box_length])
axes_sk.set_zlim3d([0, box_length])

axes_sk.view_init(
                   elev=20, 
                   azim=-20
                  )

axes_sk.set_aspect('auto')
axes_sk.set_axis_off()

plt.legend(loc = 'right')

plt.savefig('Plots/TFSI/TFSI+H2O+Cation.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSI+H2O+Cation.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/TFSI/TFSI+H2O+Cation.jpg',format='jpg',bbox_inches='tight')
plt.clf()



# #  <span style="color:OrangeRed">Cation</span>

# In[39]:


data_cat = np.loadtxt('Coordinates/LI.txt', dtype='float')
x_cat = data_cat[:,0]
y_cat = data_cat[:,1]
z_cat = data_cat[:,2]


# In[40]:


fig_cat = plt.figure(figsize=(4,4),dpi=300)
axes_cat = fig_cat.add_axes([0.1,0.1,0.9,0.9], projection = '3d')

axes_cat.scatter3D(x_cat, y_cat, z_cat, s=50, color = 'orange', label = r'$Li^{+2}$', edgecolor = 'k')

axes_cat.legend(loc='upper right')

plt.savefig('Plots/Cation/Cation.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Cation/Cation.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Cation/Cation.jpg',format='jpg',bbox_inches='tight')
plt.clf()


# #### Slice 

# In[41]:


range_cat = (x_cat>1)*(x_cat<2)


# In[42]:


fig_cat = plt.figure(figsize=(4,4),dpi=300)
axes_cat = fig_cat.add_axes([0.1,0.1,0.9,0.9], projection = '3d')

axes_cat.scatter3D(x_cat[range_cat], y_cat[range_cat], z_cat[range_cat], s=50, color = 'orange', label = r'$Li^{+2}$', edgecolor = 'k')

axes_cat.view_init(elev=0, azim=0)

axes_cat.set_xlim3d([0,box_length])
axes_cat.set_ylim3d([0,box_length])
axes_cat.set_zlim3d([0,box_length])

axes_cat.set_axis_off()

axes_cat.legend(loc='upper right')

plt.savefig('Plots/Cation/CationSlice.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Cation/CationSlice.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Cation/CationSlice.jpg',format='jpg',bbox_inches='tight')
plt.clf()


# ### With isolated water 

# In[43]:


fig_iso = plt.figure(figsize=(19.60,10.80),dpi=300)
axes_iso = fig_iso.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes_iso.scatter3D(x_cat,y_cat,z_cat, s= 2000, color = 'orangered', label = r'$Li^{+2}$', edgecolor = 'k')
axes_iso.scatter3D(OW[:,0][Zero_HB],OW[:,1][Zero_HB],OW[:,2][Zero_HB], s = 50, color = 'cyan', label = r'$H_2O(iso)$', edgecolor = 'k')

axes_iso.set_axis_off()
axes_iso.legend(loc='upper right')

plt.savefig('Plots/Cation/Cation+H2O.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Cation/Cation+H2O.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Cation/Cation+H2O.jpg',format='jpg',bbox_inches='tight')
plt.clf()

# ### Slice 

# In[44]:


range_iso = (OW[Zero_HB][:,0]>1)*(OW[Zero_HB][:,0]<2) 
range_iso2 = (LI[:,0]>1)*(LI[:,0]<2)


# In[45]:


fig_iso = plt.figure(figsize=(19.60,10.80),dpi=300)
axes_iso = fig_iso.add_axes([0.1,0.1,0.9,0.9],projection='3d')

axes_iso.scatter3D(x_cat[range_iso2],y_cat[range_iso2],z_cat[range_iso2], s= 50, color = 'orangered', label = r'$Li^{+2}$', edgecolor = 'k')
axes_iso.scatter3D(OW[:,0][Zero_HB][range_iso],OW[:,1][Zero_HB][range_iso],OW[:,2][Zero_HB][range_iso], s = 200, color = 'cyan', label = r'$H_2O(iso)$', edgecolor = 'k')

axes_iso.set_xlim3d([0,box_length])
axes_iso.set_ylim3d([0,box_length])
axes_iso.set_zlim3d([0,box_length])

axes_iso.set_axis_off()
axes_iso.legend(loc='upper right')

plt.savefig('Plots/Cation/CationH2Oslice.svg',format='svg',bbox_inches='tight')
plt.savefig('Plots/Cation/CationH2Oslice.eps',format='eps',bbox_inches='tight')
plt.savefig('Plots/Cation/CationH2Oslice.jpg',format='jpg',bbox_inches='tight')
plt.clf()