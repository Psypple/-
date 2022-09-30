# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 10:14:50 2022

@author: psypple
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xlwt

'''定义初始参数'''
table=pd.read_excel("structure.xlsx")
λ=table['波长（nm）'].values.tolist()[0]*1E-9
k0=2*np.pi/λ+0j
b=1E-6 #迭代用参数，决定初始步长，太大会导致跳出波导解
δ=1E-10 #迭代用参数
ε=1E-16 #迭代用参数
iter_max=1000 #设置一个最大迭代步数
active_index_list=table['是否为有源区（是为1）'].values.tolist() #记录是否为有源区
active_index=-1
for i in range(len(active_index_list)):
    if active_index_list[i]==1:
        active_index=i

'''输入折射率和厚度'''
n_real=table['折射率实部'].values.tolist()
n_imag=table['吸收（/cm）'].values.tolist()
for i in range(len(n_imag)):
    n_imag[i]=n_imag[i]*100*λ/(4*np.pi)
global n,d_m
n=[1+0j] #折射率，注意这个数组index=0时为空气，index最大时为衬底
for i in range(len(n_imag)):
    n.append(n_real[i]+n_imag[i]*1j)
d=table['层厚（nm）'].values.tolist() #每层厚度（从上往下）
del d[-1] #删除多余的nan项，这一项是因为无限衬底的厚度没有输入值
d_m=d #约化为国际标准单位m
for i in range (len(d_m)):
    d_m[i]=d_m[i]/1E9

'''三层波导试探解'''
d_waveguide=table['波导厚度（nm）'].values.tolist()[0]*1E-9 #波导总厚度
n1=table['cladding折射率实部'].values.tolist()[0] #cladding折射率，注意只要实部
n2=table['波导折射率实部'].values.tolist()[0] #波导折射率，注意只要实部
R_square=k0**2*d_waveguide**2/4*(n2**2-n1**2)
x_0=[]
δy_even=[]
δy_odd=[]
for i in range(10000):
    x_0.append(i*np.sqrt(R_square)/10000)
    δy_even.append(abs(np.sqrt(R_square-x_0[i]**2)-x_0[i]*np.tan(x_0[i])))
    δy_odd.append(abs(np.sqrt(R_square-x_0[i]**2)+x_0[i]/np.tan(x_0[i])))
n0=[]
for i in range(9998):
    if δy_even[i+1]<δy_even[i] and δy_even[i+1]<δy_even[i+2]:
        n0.append(np.sqrt(n2**2*k0**2-4*x_0[i+1]**2/d_waveguide**2))
    if δy_odd[i+1]<δy_odd[i] and δy_odd[i+1]<δy_odd[i+2]:
        n0.append(np.sqrt(n2**2*k0**2-4*x_0[i+1]**2/d_waveguide**2))
for i in range(len(n0)):
    n0[i]=n0[i]/k0 #这里n0记录了所有的有效折射率试探解的解析解
print('Trial solution(s) for neff:')
for i in range(len(n0)):
    print('%f'%n0[i].real)
print ('substrate refractive index: %f' %n_real[-1])

'''传递矩阵'''
def mode_profile(β):
    layer_num=len(d)
    γ=[] #γ即生长方向的各层中光场的复波矢
    for i in range(layer_num+2):
        γ.append(np.sqrt(β**2-k0**2*n[i]**2))
    AB=[np.array([1,0])] #A和B为光场通解的系数，AB列表记录了每一层（从0层空气开始）的每一组系数[A,B]
    for i in range(layer_num): #迭代layer_num次后，AB的最后一个组[Az-1,Bz-1]其实是倒数第二层，还没有得到无限衬底的系数
        Ti=np.array([[(1+γ[i]/γ[i+1])*np.exp(γ[i+1]*d_m[i])/2 , (1-γ[i]/γ[i+1])*np.exp(γ[i+1]*d_m[i])/2],
                     [(1-γ[i]/γ[i+1])*np.exp(-γ[i+1]*d_m[i])/2 , (1+γ[i]/γ[i+1])*np.exp(-γ[i+1]*d_m[i])/2]])
        AB_i1=np.dot(Ti,AB[i])
        AB.append(AB_i1)
    Tzmin1=np.array([[(1+γ[layer_num]/γ[layer_num+1])/2 , (1-γ[layer_num]/γ[layer_num+1])/2],
                     [(1-γ[layer_num]/γ[layer_num+1])/2 , (1+γ[layer_num]/γ[layer_num+1])/2]])
    AB.append(np.dot(Tzmin1,AB[-1])) #这里加入的是无限衬底层的系数，但是参考点与前面不同，计算光场时要单独处理
    return AB[layer_num+1][0],AB,γ #输出的第一项是Az，由于设置了A0=0，因此Az=t11用于迭代求解；后两项用于计算光场

'''迭代'''
def downhill(β_init):
    β_iter=[β_init] #记录迭代中的β，最后一项是最终符合收敛条件的β，即解
    t11=abs(mode_profile(β_init)[0])
    t11_0=[t11] #记录迭代中的t11，最后一项是最终符合收敛条件的t11<δ
    δβ=b*t11
    crimdef=0 #这个参数叫做前科值，即上次的δβ过大造成t11迭代没有减小留下记录，下次迭代若成果缩小则消除记录，但当此不减半δβ
    count=0
    while t11>=δ and count<iter_max:
        count=count+1
        β_temp=[] #方便对照四个t11'中哪个是下降最快的解
        t11_1=[] #记录四个试探解的t11的模
        t11_11=mode_profile(β_iter[-1]+δβ)[0]
        t11_12=mode_profile(β_iter[-1]-δβ)[0]
        t11_13=mode_profile(β_iter[-1]+δβ*1j)[0]
        t11_14=mode_profile(β_iter[-1]-δβ*1j)[0]
        t11_1.append(abs(t11_11))
        t11_1.append(abs(t11_12))
        t11_1.append(abs(t11_13))
        t11_1.append(abs(t11_14))
        β_temp.append(β_iter[-1]+δβ)
        β_temp.append(β_iter[-1]-δβ)
        β_temp.append(β_iter[-1]+1j*δβ)
        β_temp.append(β_iter[-1]-1j*δβ)
        if min(t11_1)<t11:
            β_iter.append(β_temp[t11_1.index(min(t11_1))])
            t11=t11_1[t11_1.index(min(t11_1))]
            t11_0.append(t11)
            if crimdef==0:
                δβ=1.1*δβ
            crimdef=0
        else:
            crimdef=1
            if not δβ<ε:
                δβ=0.5*δβ
    if t11>=δ:
        print ('Failed to convert.')
    else:
        print ('Converting succeeds.')
    for i in range (len(β_iter)):
        β_iter[i]=β_iter[i]/k0
    return β_iter,t11_0 #此时输出的是有效折射率

'''主控制脚本'''
def main(neff_0,q): #主程序要求输入迭代初始的有效折射率，注意初始值非常关键！特别是存在高阶模时，输入的第二项q是为了记录解的数量保存图片用
    q=q+1
    matrix=downhill(neff_0*k0) #导出迭代步骤数据
    β_final=matrix[0] #迭代过程中有效折射率中间数值
    t11=matrix[1] #迭代过程中t11中间数值
    main_ABγ=mode_profile(β_final[-1]*k0)
    iter_num=[]
    for i in range(len(β_final)):
        iter_num.append(i)
    # plt.figure()
    # plt.plot(iter_num,β_final)
    # plt.xlabel('iter number')
    # plt.ylabel('neff_iter')
    # plt.show()
    plt.figure()
    plt.title('Solution %d' %q)
    plt.plot(iter_num,t11)
    plt.xlabel('iter number')
    plt.ylabel('t11_iter')
    plt.yscale('log')
    # plt.savefig('Convergence for Solution %d.png' %q)
    plt.show()
    return main_ABγ,β_final[-1]

'''主求解程序'''
for q in range(len(n0)):
    q1=q+1
    ABγ=main(n0[q],q) #把解系数导出
    AB=ABγ[0][1]
    A=[]
    B=[]
    γ=ABγ[0][2]
    for i in range(len(n)):
        A.append(AB[i][0])
        B.append(AB[i][1])
    # print ('Az=%e + %e i' %(A[-1].real,A[-1].imag)) #输出t11的值判断精确度
    # print ('Bz=%e + %e i' %(B[-1].real,B[-1].imag)) #无限衬底开始处电场强度
    print ('neff=%f + %f i' %(ABγ[1].real,ABγ[1].imag))
    
    '''绘制光场'''
    t=[0] #设置边界点，计算光场
    for i in range(len(d)):
        t.append(t[i]+d[i])
    x=[] #x是坐标，用于绘制光场，单位nm
    E=[] #E是光场
    n_plot=[] #画光场轮廓
    integral_active=0 #计算限制因子分子的积分
    integral=0 #计算限制因子分母的积分
    for i in range(200): #上空气层光场
        x.append(1E-9*(-199+i))
        E.append(np.exp(γ[0]*(i-199)*1E-9))
        n_plot.append(n[0].real)
    for i in range(len(d)): #计算有限厚度层光场
        for k in range(500):
            x.append(t[i]+k*d[i]/500)
            Elec_field=A[i+1]*np.exp(γ[i+1]*(k*d[i]/500-d[i]))+B[i+1]*np.exp(-γ[i+1]*(k*d[i]/500-d[i]))
            E.append(Elec_field)
            n_plot.append(n[i+1].real)
            if active_index>0: #设置存在有源区时计算限制因子
                if i==active_index:
                    integral_active=integral_active+d[i]/500*abs(Elec_field)**2
    for i in range(1000): #计算无限衬底光场
        x.append(t[-1]+i*1E-9)
        E.append(B[-1]*np.exp(-γ[-1]*i*1E-9)) #从1.1版本开始不再计算设计Az的项，因为当衬底折射率远小于材料时，γ值很大，即使Az收敛至很小，随着距离增加也会导致光场强度发散，而继续降低Az对材料内模式影响很小，意义不大
        n_plot.append(n[-1].real)
    E_real=[]
    E_imag=[]
    I_modal=[] #复光场的模
    for i in range(len(E)):
        E_real.append(E[i].real)
        E_imag.append(E[i].imag)
        I_modal.append(abs(E[i])**2)
    for i in range(len(E)-1):
        integral=integral+I_modal[i]*(x[i+1]-x[i])
    for i in range(len(x)):
        x[i]=x[i]*1E6
    Γ=integral_active/integral
    print ('Γ=%f' %Γ)
    fig, ax1=plt.subplots(1,1)
    ax2=ax1.twinx()
    plt.title('Solution %d, Γ=%f' %(q1,Γ))
    ax1.plot(x,I_modal,color='red')
    ax1.set_xlabel('Distance (μm)')
    ax1.set_ylabel('Intensity (arb.unit)',color='red')
    ax1.tick_params(axis='y',colors='red')
    ax2.plot(x,n_plot,color='blue')
    ax2.set_ylabel('refractive index',color='blue')
    ax2.set_ylim(0,max(n_plot)*1.5)
    ax2.tick_params(axis='y',colors='blue')
    ax2.spines['left'].set_color('red')
    ax2.spines['right'].set_color('blue')
    plt.savefig('Modal Intensity for Solution %d.png' %q1)
    plt.show()
    
    '''打印结果'''
    f = xlwt.Workbook() #创建工作薄
    sheet1 = f.add_sheet(u'sheet1',cell_overwrite_ok=True) #创建sheet
    list = x
    list.insert(0,'distance (μm)')
    list.insert(0,'neff')
    list.insert(0,'Γ')
    j = 0
    for i in list:
        sheet1.write(j,0,i) #循环写入
        j=j+1
    list = I_modal
    list.insert(0,'Intensity (arb.unit)')
    list.insert(0,'%f + %f i' %(ABγ[1].real,ABγ[1].imag))
    list.insert(0,'%f' %Γ)
    j = 0
    for i in list:
        sheet1.write(j,1,i) #循环写入
        j=j+1
    f.save('solution%d.xls' %q1)#保存文件
input('Press any button to exit.')