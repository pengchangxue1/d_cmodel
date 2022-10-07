import itasca as it
it.command("python-reset-state false")#如果没有给出此命令，则在通过model new或model restore命令重置FLAC3D模型状态时，将删除Python变量和函数。
it.command("""
program load module 'zone'
program load guimodule 'zone'
program load module 'body'
program load guimodule 'body'
program load module 'extruder'
program load guimodule 'extruder'
program load module 'sel'
program load module 'zonesel'
program load module 'selpython'
program load module 'wallsel'
program load guimodule 'sel'
program load module 'pfcsel'
program load module 'rblocksel'
program load module 'dfnzone'
program load module 'ballzone'
program load module 'zonepython'
program load module 'wallzone'
""")
import math


Pa=101325.0
c=110e3
fi0=48.5
delt_fi=0
Rf=0.79
Ki=704
n=0.38
Kb=303
m=0.18
Kur=844.8

segm1=-it.zone.near((0,0,0.6)).stress_prin_x()
segm3=-it.zone.near((0,0,0.6)).stress_prin_z()

fi=fi0-delt_fi*math.log((segm3+1)/Pa,10)
fi=fi0

fi=fi*math.pi/180

s=(1-math.sin(fi))*(segm1-segm3)/(2*c*math.cos(fi)+2*segm3*math.sin(fi))
Ei=Ki*Pa*(segm3/Pa)**n  
Emin=0.25*Ki*Pa*(0.02)**n
Et=(1-Rf*s)**2*Ei

Eur=Kur*Pa*(segm3/Pa)**n
Kt=Kb*Pa*(segm3/Pa)**m 
if Et<=Emin :
    Et=Emin


if Kt<Et/3 :
    Kt=Et/3 
    
if Kt>17*Et :
    Kt=17*Et

Gt=3*Kt*Et/(9*Kt-Et)


#print(str(Kt)+','+str(Ei)+','+str(it.zone.near((0,0,0.6)).stress_prin_x())+','+str(it.zone.near((0,0,0.6)).stress_prin_z()))
#
print(Et)
print(Eur)
#print(it.zone.near((1,1,1)).props())
print(it.zone.near((0,0,0.6)).prop('young'))
print('_____________')
#print(it.zone.near((1,1,1)).prop('poisson'))

print(Kt)
print(it.zone.near((0,0,0.6)).prop('bulk'))
print('_____________')
print(Gt)
print(it.zone.near((0,0,0.6)).prop('shear'))

print('_____________')

print(it.zone.near((0,0,0.6)).stress_prin_x())
print(it.zone.near((0,0,0.6)).stress_prin_y())
print(it.zone.near((0,0,0.6)).stress_prin_z())
print('#############')