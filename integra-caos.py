### Esse programa le dados de estrelas e faz um mapa de densidade
### Dados do GAIA DR2, assim como um mapa dinâmico com o uso dos
### Expoentes de Lyapunov
###
### Autor: Lucas de Amorim (IAG/USP) - Data: 18/06/2020

### Para contar o tempo que o programa levou para rodar
import time
start_time = time.time()
print(time.strftime("\n%H:%M:%S", time.localtime()), "Início do programa.\n\n")

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
# from astropy.coordinates import Galactic
# from astropy.coordinates import ICRS
from astropy.coordinates import Galactocentric
# import winsound

### Pacotes que usaremos depois, para a integração
from galpy.potential import LogarithmicHaloPotential
from galpy.potential import SpiralArmsPotential
# from galpy.potential.JunqueiraPotential import JunqueiraPotential
from galpy.orbit import Orbit
from galpy.potential import plotRotcurve

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd
import seaborn as sns

### Define o potencial do halo logarítimico
lhp = LogarithmicHaloPotential()
sap = SpiralArmsPotential(N = 4, omega = 25)


### Faz o gráfico das curvas de rotação
# plotRotcurve(lhp, Rrange=[0.01,10.],grid=1001,yrange=[0.,1.2])
# plotRotcurve(sap, Rrange=[0.01,10.],grid=1001,yrange=[0.,1.2], phi = 0.)
# plotRotcurve(lhp + sap, Rrange=[0.01,10.],grid=1001,yrange=[0.,1.2], phi = 0.)
# plt.show()


### Leitura do arquivo com os dados
lines = [line.rstrip("\n") for line in open("./gaia.csv")]
### Os dados são, em ordem, ra, dec, pmra, pmdec, parallax e radial_velocity
ra, dec, pmra, pmdec, plx, rv, R, distance = ([] for i in range(8))


N = 200

divisoes = 1

numcores = 4    ### Max = 64

for line in lines[1:N]:
    table = line.split(",")
    ra_i = float(table[0])
    dec_i = float(table[1])
    pmra_i = float(table[2])
    pmdec_i = float(table[3])
    plx_i = float(table[4])
    rv_i = float(table[5])

    ### Queremos apenas as estrelas com distância menor que 200 pc do Sol
    distance_check = 1/plx_i * 10 ** 3

    if distance_check < 200:

        distance.append(1/plx_i * 10 ** 3)
        ra.append(ra_i)
        dec.append(dec_i)
        pmra.append(pmra_i)
        pmdec.append(pmdec_i)
        plx.append(plx_i)
        rv.append(rv_i)


### Definimos as coordenadas com unidades com o SkyCoord
c = SkyCoord(ra*u.degree, dec*u.degree, distance = distance*u.pc, \
            pm_ra_cosdec = pmra*u.mas/u.yr, pm_dec = pmdec*u.mas/u.yr, \
            radial_velocity = rv*u.km/u.s, frame = 'icrs')

### Transformamos as coordenadas de ICRS para Galactocentrica
gc = c.transform_to(coord.Galactocentric)

##################### Sem limites nas velocidades ###########################
### Pegamos as componentes x e y da velocidade cartesiana
# vx = np.array(gc.v_x)
# vy = np.array(gc.v_y)
# # vz = np.array(gc.v_z)

#############################################################################


##################### Com limites nas velocidades ###########################
### Pegamos as componentes x e y da velocidade cartesiana
U = np.array(gc.v_x)
V = np.array(gc.v_y)
# W = np.array(gc.v_z)


### Limitamos os limites aceitados para as velocidades
# cond = np.where((U > -80) & (U < 80) & (V > 170) & (V < 260))
cond = np.where(U < 10000000)
vx = U[cond]
vy = V[cond]

#############################################################################



### Definimos os passos de integração
ts = np.linspace(0., 10 ** 3., 1000)*u.Myr



### Inicializamos as órbitas, com a posição do Sol, e a velocidade
### encontrada da estrela
### Vetores com as coordenadas de posição do Sol
x = [8.3]*(len(vx))
y = [0]*(len(vx))
# z = [0]*(len(vx))
z = [27]*(len(vx))


coord_R = SkyCoord(x = x*u.kpc, y = y*u.kpc, z = z*u.pc, \
              v_x = vx*u.km/u.s, v_y = vy*u.km/u.s, v_z = 0*u.km/u.s, \
              frame = 'galactocentric', representation_type='cartesian')
o = Orbit(coord_R)

### Órbita com pequena diferença de raio para o cálculo do
### expoente de Lyapunov
### Vetor com as coordenadas de posição do Sol dR
x_dR = [8.2]*(len(vx))

coord_dR = SkyCoord(x = x_dR*u.kpc, y = y*u.pc, z = z*u.pc, \
              v_x = vx*u.km/u.s, v_y = vy*u.km/u.s, v_z = 0*u.km/u.s, \
              frame = 'galactocentric', representation_type='cartesian')
odR = Orbit(coord_dR)


###########################    Jeito Demorado    ##############################
# orb = o[0]
# orb_dR = odR[0]
# orb.integrate(ts, lhp+sap, method='rk6_c')
# orb_dR.integrate(ts, lhp+sap, method='rk6_c')
# MLCE_matrix = np.array(np.log(abs(orb.R(ts)-orb_dR.R(ts))/(x[0] - x_dR[0])))
# # print(MLCE_matrix, MLCE_matrix.size)
# for i in range(1, len(vx)):
#     orb= o[i]
#     orb_dR = odR[i]
#     orb.integrate(ts, lhp+sap, method='rk6_c')
#     orb_dR.integrate(ts, lhp+sap, method='rk6_c')
#
#     MLCE_matrix = np.array(np.append(MLCE_matrix, np.array(np.log(abs(orb.R(ts)-orb_dR.R(ts))/(x[0] - x_dR[0]))), axis = 0))
#     # print(MLCE_matrix)
#
#     # orb.plot()
#
# MLCE_matrix = np.reshape(MLCE_matrix, (len(vx), -1))
# MLCE = np.nanmax(MLCE_matrix, 1)

# print(time.strftime("%H:%M:%S", time.localtime()),
#     "Fim da integração. Preparando os gráficos. \n")

###############################################################################



#########################    Jeito Menos Demorado    ##########################

divis = [round(len(vx)*i/divisoes) for i in range(1, divisoes+1)]

print(f"Integrando as {len(vx)} órbitas em {divisoes} partes.\n")
print(time.strftime("%H:%M:%S", time.localtime()), "Intervalo 1:", "0", divis[0])

orb = o[0:divis[0]+1]
orb_dR = odR[0:divis[0]+1]
orb.integrate(ts, lhp+sap, method='odeint', numcores = numcores)

orb_dR.integrate(ts, lhp+sap, method='odeint', numcores = numcores)
MLCE_matrix = (np.log(abs(orb.R(ts)-orb_dR.R(ts))/(x[0] - x_dR[0])))/len(ts)
# print(orb.size)
# print(MLCE_matrix, MLCE_matrix.size)
for i in range(1, len(divis)):
    # print(o[divis[i-1]].vxvv)
    # print(o[divis[i-1]+1].vxvv)
    print(time.strftime("%H:%M:%S", time.localtime()), f"Intervalo {i+1}:", (divis[i-1]+1), divis[i])
    orb = o[(divis[i-1]+1):(divis[i]+1)]
    orb_dR = odR[(divis[i-1]+1):(divis[i]+1)]
    orb.integrate(ts, lhp+sap, method='rk6_c', numcores = numcores)
    orb_dR.integrate(ts, lhp+sap, method='rk6_c', numcores = numcores)

    MLCE_matrix = np.array(np.append(MLCE_matrix, (np.array(np.log(abs(orb.R(ts)-orb_dR.R(ts))/(x[0] - x_dR[0]))))/len(ts), axis = 0))
    # print(o[divis[i]].vxvv)

MLCE_matrix = np.reshape(MLCE_matrix, (len(vx), -1))
MLCE = np.nanmax(MLCE_matrix, 1)


print(time.strftime("%H:%M:%S", time.localtime()),
    "Fim da integração. Preparando os gráficos. \n")
###############################################################################




##############################  Jeito Rápido  #################################

# ### Fazemos a integração
# o.integrate(ts, lhp+sap, method='rk6_c')
# odR.integrate(ts, lhp+sap, method='rk6_c')
#
# ### Calculando o Maximum Lyapunov Characteristc Exponent
# MLCE_matrix = np.log(abs(o.R(ts)-odR.R(ts))/(x[0] - x_dR[0]))
# MLCE = np.nanmax(MLCE_matrix, 1)
#
#
# print(time.strftime("%H:%M:%S", time.localtime()),
#     "Fim da integração. Preparando os gráficos. \n")
###############################################################################


# import matplotlib.image as mpimg
# img = mpimg.imread('your_image.png')
# imgplot = plt.imshow(img)
# from numpy import diff

for i in range(len(vx)):

    # # print(MLCE_matrix[i][:])
    # plt.figure(figsize=(20,7))
    # plt.subplot(2, 1, 1)

    # plt.plot(ts, orb[i].R(ts))
    # # lyap[i] = (1 / Niter) * sum( np.log( abs( dfdx(x_lyap, p[i]) ) ) )
    # lyap = 1/len(ts) * sum(np.log( abs (1)))
    #
    #
    # plt.subplot(1, 2, 1)
    plt.plot(orb[i].x(ts), orb[i].y(ts))
    #
    #
    # plt.subplot(1, 2, 2)
    # plt.plot(MLCE_matrix[i][:], 'r')
    # # plt.plot(t, np.real(square_serie(99, frequencias[i])))
    plt.suptitle(f'i = {i}')
    # # plt.xlim(-0.05, 0.05)
    # # print(frequencias[i])
    #
    #
    # plt.pause(0.5)
    # plt.close()

plt.show()




####### Construção do histograma 2D ########
print("Número de estrelas: ", len(vx))

x_inf, x_sup = -80, 80
y_inf, y_sup = 170, 260

# print(x_inf, x_sup)

plt.figure(1)
## Configurações do histograma 2D UV (bins de 1 km/s)
plt.hist2d(vx, vy, bins=((x_sup-x_inf), (y_sup-y_inf)), cmap=plt.cm.jet)
plt.colorbar(matplotlib.cm.ScalarMappable(cmap=plt.cm.jet))
# plt.colorbar()
## Caso não queria estabelecer os limites:
# plt.hist2d(vx, vy, bins=(300, 300), cmap=plt.cm.jet)

plt.xlabel('U [km/s]', fontsize=14)
plt.ylabel('V [km/s]', fontsize=14)

plt.xlim(x_inf, x_sup)
plt.ylim(y_inf, y_sup)

# print("nbins = ", x_sup-x_inf, (y_sup-y_inf)*2)
plt.savefig('MapaVel.png')



plt.figure(2)

data = np.array([vx, vy])
data = np.reshape(data, (2, len(vx))).T
data = pd.DataFrame(data)


data['z'] = list(MLCE)
# print(data.info())

nCut = [(x_sup-x_inf)*2 - 1, (y_sup-y_inf)*3 - 1]

cuts = pd.DataFrame({str(feature) + 'Bin' : pd.cut(data[feature], nCut[feature]) for feature in [0, 1]})

cuts.columns = ["U [km/s]", "V [km/s]"]

# print(cuts)
# print('at first cuts are pandas intervalindex.')
# print(cuts.head())
# print(cuts.info())


means = data.join(cuts).groupby( list(cuts) ).mean()
means = means.unstack(level = 0) # Use level 0 to put 0Bin as columns.

# Reverse the order of the rows as the heatmap will print from top to bottom.
means = means.iloc[::-1]

# arr_vx = means['z'].columns.map(lambda x : x.left).tolist()
# means['z'] = means['z'].astype('float')

# means[means.select_dtypes(['category']).columns] = means[means.select_dtypes(['category']).columns].apply(lambda x : x.left)

# print(means.head())
# print(means['z'])
# print("***************************")
# print(means['z'].columns.map(lambda x : x.left).tolist()[::20])
# print(means['z'].columns.codes)



ax = sns.heatmap(means['z'], xticklabels = 23,
                        yticklabels = 26,
                        cmap = "Purples")
# ax = sns.heatmap(means['z'], xticklabels = means['z'].columns.map(lambda x : x.left),
#                         yticklabels = means['z'].index.map(lambda x : x.left),
#                         cmap = "Purples")
# plt.title('Means of z vs Features 0 and 1')



# print(ax.get_xlim(), ax.get_ylim())

# xlim = [x_inf, x_sup]
# ylim = [y_inf, y_sup]
#
# ax.set_xlim((110.8+x_inf)*1.222, (110.8+x_sup)*1.222)
# ax.set_ylim((273-y_inf), (273-y_sup))

ax.set_xlabel('U [km/s]', fontsize=14)
ax.set_ylabel('V [km/s]', fontsize=14)

plt.tight_layout()

labels_x = means['z'].columns.map(lambda x : x.left).tolist()[::23]

labels_y = means['z'].index.map(lambda x : x.left).tolist()[::26]

labels_x = [int(round(num, -1)) for num in labels_x]
labels_y = [int(round(num)) for num in labels_y]

ax.set_xticklabels(labels_x)

ax.set_yticklabels(labels_y)

plt.savefig('MapaDin.png')


### Sol: plt.plot(11.1, 232.24, marker= '*', color = 'white')

### Mostra quanto tempo o programa demorou para rodar
if ((time.time() - start_time) > 60.):
    if ((time.time() - start_time) > 3600.):
        print("--- %d hours ---" % (float(time.time() - start_time) // 3600))
        print("--- %d minutes ---" % (float(time.time() - start_time) % 3600 // 60))
        print("--- %s seconds ---" % (float(time.time() - start_time) % 3600 % 60 ))
    else:
        print("--- %d minutes ---" % (float(time.time() - start_time) // 60))
        print("--- %s seconds ---" % (float(time.time() - start_time) % 60))
else:
    print("--- %s seconds ---" % (time.time() - start_time))
# winsound.Beep(2300, 400)
# winsound.Beep(2300, 400)
# winsound.Beep(2300, 400)

# plt.show()
